#!/usr/bin/env python3
import argparse
import os
import time
import random
import subprocess
import pickle
import bisect
import gzip
import shutil
import hashlib
from collections import defaultdict, Counter
import pysam
from tqdm import tqdm


def verify_tools():
    """
    Determine if the needed tools are on path
    If any of the rquired tools is not on PATH kill siRmap
    
    """
    samtools_path   = shutil.which("samtools")
    bowtie_path     = shutil.which("bowtie")
    bowtieb_path    = shutil.which("bowtie-build")
    pigz_path       = shutil.which("pigz")
    
    if bowtieb_path is None:
        raise RuntimeError("Bowtie1 indexer ('bowtie-build') not found in PATH.")
    if bowtie_path is None:
        raise RuntimeError("Bowtie1 is not found in PATH, please verify srnanf_env is active or call apptainer properly")
    if samtools_path is None:
        raise RuntimeError("Samtools is not found in PATH, please verify srnanf_env is active or call apptainer properly")
    if pigz_path is None:
        raise RuntimeError("pigz is not found in PATH, please verify srnanf_env is active or call apptainer properly")
        

# helper to open files
def _open_any(
    path,
    mode="rt"
    ):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# helper function to call pigz to speedup ziping steps
def _pigz_compress(
    in_file,
    out_file,
    threads=0
    ):
    """
    Removes inFile on success
    """
    # Check outs are gz just in case
    if not out_file.endswith(".gz"):
        raise ValueError("out_file must end with .gz")
    # is pigz avalable on path?
    if shutil.which("pigz") is None:
        raise RuntimeError("pigz not found in PATH")
    # Open in and out files
    with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
        cmd = ["pigz", "-c"] # base pigz line
        # if to check the numebr of threads
        if threads and threads > 0:
            cmd = ["pigz", "-p", str(threads), "-c"]
        # Call pigz
        proc = subprocess.Popen(cmd, stdin=fin, stdout=fout)
        # wait till pigz ends and capture rc as a control
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"pigz failed compressing {in_file} to {out_file}")
    os.remove(in_file)


# collapse identical reads only
def collapse_fastq_to_fastq(
    fastq_path,
    out_fastq,
    rc_map_pickle,
    pigz_threads=0,
    summary_log_file=None
    ):
    """
    Collapse identical FASTQ reads to one FASTQ seq. This is helpfull for hd libraries.
    - Quality: taken from the first instance of each distinct sequence
      (keeps meaninful qualities and avoids inventing them)
    - Saves {read_id: RC} in rc_map_pickle.
    Returns (num_unique, num_total, col_gz)
    """
    print(f"\nCollapsing FASTQ: {fastq_path}")

    # Read counter
    seq_counts = {}
    total = 0
    # parse the fastq file in chunks of 4 lines
    with _open_any(fastq_path, "rt") as fh, tqdm(desc="Collapsing FASTQ", unit="reads", total=None) as pbar:
        while True:
            # header; if empty, break
            h = fh.readline()
            if not h:
                break
            # sequence
            seq = fh.readline().rstrip("\n\r")
            # plus line
            fh.readline()   # '+'
            # quality
            qual = fh.readline().rstrip("\n\r")
            # increment count of total reads seen
            total += 1
            if seq in seq_counts:
                seq_counts[seq][0] += 1
            else:
                # store first quality seen for that sequence
                seq_counts[seq] = [1, qual]
            pbar.update(1)

    # Now write a collapsed fastq: one record per unique sequence,
    # assigning "rand" read IDs
    readid_to_rc = {}

    # Wirte tmp and gz collapsed reads
    col_gz = out_fastq if out_fastq.endswith(".gz") else (out_fastq + ".gz")
    tmp_plain  = col_gz[:-3]  # strip .gz for the temporary plain file


    # Add ID
    with open(tmp_plain, "wt") as fq:
        for i, (seq, (rc, qual)) in enumerate(seq_counts.items(), start=1):
            # I will add a hash to 100% avoid possible duplicated read ids
            hs = hashlib.blake2b(seq.encode(), digest_size=5).hexdigest()
            rid = f"rid{hs}_{i:06d}"  # new id
            readid_to_rc[rid] = rc  # multiplicity for later RC tag
            fq.write(f"@{rid}\n{seq}\n+\n{qual}\n")  # write new fastq

    # Save IDs
    with open(rc_map_pickle, "wb") as pf:
        pickle.dump(readid_to_rc, pf)

    # Mandatory pigz compression
    _pigz_compress(tmp_plain, col_gz, threads=pigz_threads)

    # append a short summary to the shared log if provided
    if summary_log_file:
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "a") as sf:
            sf.write(f"Collapsed: {len(seq_counts):,} unique of {total:,} total reads\n")
            sf.write(f"FASTQ written: {col_gz}\n")
            sf.write(f"RC map saved: {rc_map_pickle}\n\n")

    return len(seq_counts), total, col_gz

def parse_read_id_for_rc(read_name, rc_map):
    return rc_map.get(read_name, 1)


# MAPPING STEPS
# INDEX
def build_bowtie_index(
    genome_fasta,
    index_base, 
    offrate=4,
    threads=1):
    """
    The indexing step will consider the offrate option for large/repetitive genomes
    """
    # cmd to define bowtie-build
    cmd = ["bowtie-build", "--offrate", str(offrate), "--threads", str(threads), genome_fasta, index_base]
    log_file = f"{index_base}_bowtie_build.log"
    print(f"\nBuilding Bowtie index with command:\n{' '.join(cmd)}")
    print(f"Logging Bowtie output to: {log_file}")
    # Run and save both st outputs
    with open(log_file, "w") as lf:
        process = subprocess.Popen(cmd, stdout=lf, stderr=lf)
        process.wait()

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)
    print(f"Index built successfully: {index_base}.*.ebwt\n")

# Alignment
def run_bowtie(
    fastq_file,
    index_base,
    output_bam=None,
    mismatches=1,
    save_non_mappers=False,
    non_mappers_file=None,
    threads=4,
    sort_mem="2G",
    index_output_bam=True,
    summary_log_file=None,
    ):
    """
    Map reads with Bowtie1 and produce a sorted BAM.
    Optionally save unmapped reads using Bowtie --un
    Summarize unique vs multi from the final BAM.
    """
    # Files will be save using basename
    base_name = os.path.basename(fastq_file)
    # Safe step for extensions
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base_name.endswith(ext):
            base_name = base_name[: -len(ext)]
            break

    if output_bam is None:
        output_bam = base_name + ".bam"

    # Clean bowtie index ID with basename (strip .fa/.fasta if passed)
    index_base = os.path.splitext(index_base)[0] if index_base.endswith((".fa", ".fasta")) else index_base

    # Define logs
    raw_log_file = f"{base_name}_raw_bowtie.log"
    if summary_log_file is None:
        summary_log_file = f"{base_name}_summary.log"  # default only if not provided

    # Build Bowtie command
    cmd = [
        "bowtie",
        "-S",
        "-q",
        "-v", str(mismatches),
        "-a",
        "--best", "--strata",
        "-p", str(threads),
    ]

    # Save unmaped
    tmp_un_plain = None
    unmapped_gz_path = None
    if save_non_mappers and non_mappers_file:
        # ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(non_mappers_file)) or ".", exist_ok=True)
        unmapped_gz_path = non_mappers_file if non_mappers_file.endswith(".gz") else (non_mappers_file + ".gz")
        tmp_un_plain = unmapped_gz_path[:-3]  # write plain FASTQ first
        cmd += ["--un", tmp_un_plain]

    # Add index and reads
    cmd += [index_base, fastq_file]

    # Note: I have had issues with samtools before when using ShortStack 3.8.5 becasue of the sort mem thus:
    sort_mem_value = sort_mem.upper()
    if sort_mem_value.endswith("G"):
        total_mem_bytes = int(float(sort_mem_value[:-1]) * 1024 * 1024 * 1024)
    elif sort_mem_value.endswith("M"):
        total_mem_bytes = int(float(sort_mem_value[:-1]) * 1024 * 1024)
    else:
        total_mem_bytes = int(sort_mem_value)
    # per-thread mem for samtools sort
    sort_mem_bytes = total_mem_bytes // max(1, threads)

    printable_cmd = " ".join(cmd)
    print("\nRunning Bowtie mapping:")
    print(printable_cmd + f" | samtools view -b - | samtools sort -m {sort_mem_bytes} -@ {threads} -o {output_bam}\n")

    with open(raw_log_file, "w") as lf:
        lf.write(f"Start time: {time.strftime('%c')}\n")
        lf.write("Running command: " + printable_cmd + "\n")

        # bowtie stdout = SAM, stderr = summary
        bowtie = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        view   = subprocess.Popen(["samtools", "view", "-b", "-@", str(threads), "-"],
                                  stdin=bowtie.stdout, stdout=subprocess.PIPE)
        sortp  = subprocess.Popen(["samtools", "sort", "-@", str(threads), "-m", str(sort_mem_bytes), "-o", output_bam, "-"],
                                  stdin=view.stdout)

        # Make sure intermediate pipes aren’t blocking
        bowtie.stdout.close()
        view.stdout.close()

        # Capture bowtie stderr into the log in real time
        for line in iter(bowtie.stderr.readline, b""):
            s = line.decode(errors="replace")
            lf.write(s)
        bowtie.stderr.close()

        # Wait for the pipeline
        rc_sort = sortp.wait()
        rc_view = view.wait()
        rc_bt   = bowtie.wait()

        lf.write(f"\nEnd time: {time.strftime('%c')}\n")
        lf.flush()

    # Check return codes
    if rc_bt != 0 or rc_view != 0 or rc_sort != 0:
        raise RuntimeError(f"bowtie/view/sort pipeline failed (bowtie={rc_bt}, view={rc_view}, sort={rc_sort}). See {raw_log_file}")

    print(f"Sorted BAM ready: {output_bam}")

    if tmp_un_plain is not None:
        print(f"Compressing unmapped reads with pigz → {unmapped_gz_path}")
        _pigz_compress(tmp_un_plain, unmapped_gz_path, threads=threads)

     # Optional bam index
    if index_output_bam:
        try:
            pysam.index(output_bam)
        except Exception:
            pass

    # Summarize alignment 
    print("\nGenerating alignment summary")
    bam_in = pysam.AlignmentFile(output_bam, "rb")

    # Let's start with collapsed reads
    alignments_per_read = defaultdict(int)
    for aln in tqdm(bam_in.fetch(until_eof=True), desc="Count alns per read ID"):
        if not aln.is_unmapped:
            alignments_per_read[aln.query_name] += 1
    bam_in.close()

    # Totals by read ID
    aligned_read_ids = len(alignments_per_read)                              
    unique_read_ids  = sum(1 for n in alignments_per_read.values() if n == 1) 
    multimapper_read_ids = sum(1 for n in alignments_per_read.values() if n > 1)

    # Read Bowtie stderr summary for total processed, this will be easier
    bt_total_processed = None
    bt_mapped_any = None
    try:
        with open(raw_log_file, "r") as lf:
            for line in lf:
                if line.startswith("# reads processed"):
                    bt_total_processed = int(line.strip().split(":")[-1].strip())
                elif line.startswith("# reads with at least one alignment"):
                    bt_mapped_any = int(line.strip().split(":")[-1].split()[0].strip())
    except Exception:
        pass

    # Add aln summary to the shared summary log (this will be helpful for nextflow)
    os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
    with open(summary_log_file, "a") as sf: 

        # Bowtie
        if bt_total_processed is not None and bt_mapped_any is not None:
            bt_unmapped = bt_total_processed - bt_mapped_any
            sf.write(
                "Bowtie (input-level): processed {0:,}, ≥1 alignment {1:,} ({2:.1f}%), "
                "unmapped {3:,} ({4:.1f}%).\n".format(
                    bt_total_processed,
                    bt_mapped_any,
                    (bt_mapped_any / max(1, bt_total_processed)) * 100.0,
                    bt_unmapped,
                    (bt_unmapped / max(1, bt_total_processed)) * 100.0
                )
            )

        # Final summary
        sf.write("BAM summary (by collapsed read ID / QNAME):\n")
        sf.write(
            "  Unique mappers: {0:,} / {1:,} ({2:.1f}%)\n".format(
                unique_read_ids, aligned_read_ids,
                (unique_read_ids / max(1, aligned_read_ids)) * 100.0
            )
        )
        sf.write(
            "  Multimappers:  {0:,} / {1:,} ({2:.1f}%)\n".format(
                multimapper_read_ids, aligned_read_ids,
                (multimapper_read_ids / max(1, aligned_read_ids)) * 100.0
            )
        )

        # If need save unmappers
        if save_non_mappers and non_mappers_file:
            out_un = non_mappers_file if non_mappers_file.endswith('.gz') else (non_mappers_file + '.gz')
            sf.write(f"Unmapped reads saved to: {out_un}\n")

        sf.write("\n")

    print(f"Mapping summary saved: {summary_log_file}\n")


# Add tags based on RC  files but also NH
def add_NH_RC_tags_from_collapse(
    input_bam,
    rc_map_pickle,
    output_bam
    ):
    """
    Add NH (alignments per collapsed read id) and RC using the saved collapse map.
    """
    print(f"\nAdding NH/RC (from collapse map) to {input_bam} ...")
    with open(rc_map_pickle, "rb") as pf:
        rc_map = pickle.load(pf)

    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    # Count #alignments per read id
    read_counts = Counter()
    for read in bam_in.fetch(until_eof=True):
        if not read.is_unmapped:
            read_counts[read.query_name] += 1
    bam_in.reset()

    for read in tqdm(bam_in.fetch(until_eof=True), desc="Writing reads with NH/RC", unit="reads"):
        if not read.is_unmapped:
            read.set_tag("NH", read_counts[read.query_name]) # use read_counts to add NH
            rc = parse_read_id_for_rc(read.query_name, rc_map) # read RC
            read.set_tag("RC", int(rc)) # add tags
        bam_out.write(read)

    bam_in.close()
    bam_out.close()
    print(f"NH/RC tags added and saved to {output_bam}\n")


# Build index with the bam list (must be provided)
# This will be done by nextflow
def build_weighted_unique_index_from_bams(
    bam_list_file,
    output_index="unique_index.pkl"
    ):
    """
    Unique mappers (NH==1) with RC weights.
    For each (chrom,strand) store sorted arrays + prefix sums (RC)
      - starts_pos, starts_pref
      - ends_pos,   ends_pref
    """
    # Read bam files from the list (this will be provided by nextflow)
    with open(bam_list_file) as f:
        bam_paths = [l.strip() for l in f if l.strip()]

    tmp = defaultdict(lambda: {"starts": [], "ends": []})  # key=(chrom,strand)

    for bf in bam_paths:
        # read bam
        bam = pysam.AlignmentFile(bf, "rb")
        # For each read ....
        for r in tqdm(bam.fetch(until_eof=True), desc=f"Indexing uniques in {os.path.basename(bf)}"):
            if r.is_unmapped:
                continue  # skip unmapped just in case
            try:
                nh = r.get_tag("NH") # how many possible positions?
            except KeyError:
                nh = 1    # Get reads with NH tag = 1 if not just drop them
            """
            To does not overestimate the local read context I will only admit
            unique mappers for the index. This is close to the -u mode of ShotStack
            Also reads a bit off the window can also count.
            """
            if nh != 1:
                continue
            # Now my index must have a weight based on the repetitivines of the reads
            # So use the RC tag. This is the most relevant step since the index will be used to score
            try:
                rc = int(r.get_tag("RC"))
            except KeyError:
                rc = 1
            # Save chr, strand, start, end
            chrom = r.reference_name
            strand = '-' if r.is_reverse else '+'
            s = max(0, r.reference_start)
            e = r.reference_end
            if e <= s:
                continue
            tmp[(chrom, strand)]["starts"].append((s, rc)) # add the rc
            tmp[(chrom, strand)]["ends"].append((e, rc))
        bam.close()

    # I need to build arrays
    index = {}
    for key, d in tmp.items():
        # Sort by position
        s_sorted = sorted(d["starts"])  # list of (pos, rc)
        e_sorted = sorted(d["ends"])
        starts_pos = [p for p, _ in s_sorted]
        ends_pos   = [p for p, _ in e_sorted]
        starts_pref = []
        ends_pref = []
        acc = 0
        for _, w in s_sorted:
            acc += w
            starts_pref.append(acc)
        acc = 0
        for _, w in e_sorted:
            acc += w
            ends_pref.append(acc)
        index[key] = {
            "starts_pos": starts_pos, "starts_pref": starts_pref,
            "ends_pos": ends_pos, "ends_pref": ends_pref
        }

    with open(output_index, "wb") as f:
        pickle.dump(index, f)

    print(f"\nWeighted unique index saved to {output_index}\n")


def weighted_unique_overlap(index, chrom, strand, start, end):
    """
    Sum of RC weights of NH==1 reads overlapping [start, end) on (chrom,strand)
    computed via prefix sums over starts/ends. Exact, O(log N)
    """
    # Checkpoints to skip non-overlapping regions or missing buckets
    if end <= start:
        return 0
    d = index.get((chrom, strand))
    if not d:
        return 0
    # Arrays to give weights based on RC
    # overlap = (∑ start < end) − (∑ end ≤ start)
    sp, sw = d["starts_pos"], d["starts_pref"]
    ep, ew = d["ends_pos"], d["ends_pref"]

    i = bisect.bisect_left(sp, end)     # starts < end
    sum_starts = sw[i-1] if i > 0 else 0

    j = bisect.bisect_right(ep, start)  # ends <= start
    sum_ends = ew[j-1] if j > 0 else 0
    # So apply the eq 
    return max(0, sum_starts - sum_ends)

# Helpers to expand the bams
def _rc_from_read(read):
    try:
        return int(read.get_tag("RC"))
    except Exception:
        return 1


# This helper function is essential becasue it helps to exapnd bams based on RC (n times)
def _propagate_n_copies(bam_out, aln, n):
    """
    Write 'n' copies of 'aln' each with a hash appended to the read ID to avoid posible duplicated values
    """
    base_qname = aln.query_name
    hdr = bam_out.header # Shared metadata
    aln_txt = aln.to_string()  # snapshot to clone safely
    for k in range(1, max(1, n) + 1):
        aln_cp = pysam.AlignedSegment.fromstring(aln_txt, hdr) # parse each line
        # Add the hash to the read qname and save
        h = hashlib.blake2b(f"{base_qname}:{k}".encode(), digest_size=3).hexdigest()
        aln_cp.query_name = f"{base_qname}_{h}"
        bam_out.write(aln_cp) # write

# Assigning multimappers
# RANDOM
def assign_multimappers_randomly(
    bam_file,
    output_bam,
    seed=1,
    sort_index=True,
    summary_log_file=None
    ):
    """
    Randomly assign a position to each multimapping collapsed read
    Propagate each read RC times (inspired on Cei's code but with out tally).
    """
    # helper to normalize tags on output records (res)
    def _sanitize_for_out(alignment):
        """
        clean the tags after propagating RC tags in the resolved BAM file
        """
        # Force NH and RC 1 per expanded instance
        try:
            alignment.set_tag("NH", 1, value_type='i')
        except Exception:
            pass

        # Same for RC
        try:
            if alignment.has_tag("RC"):
                alignment.set_tag("RC", 1, value_type='i')
        except Exception:
            pass

        # This is just in case and to standardize the out bams
        # take ride of multimapper and secondary tags 
        try:
            # pysam returns a list of (tag_name, tag_value, value_type)
            all_tags = alignment.get_tags(with_value_type=True)
            tags_to_remove = {"IH", "HI", "XA", "SA"}
            kept_tags = [
                (tag_name, tag_value, value_type)
                for (tag_name, tag_value, value_type) in all_tags
                if tag_name not in tags_to_remove
            ]
            alignment.set_tags(kept_tags)
        except Exception:
            pass

        return alignment


    # define a random seed
    # if we will assign random this must be reproducible
    start_time = time.time()
    random.seed(seed) 
    print("\nassigning multimappers randomly")

    # Read bam
    bam_in = pysam.AlignmentFile(bam_file, "rb")

    # build groups based on unique instances
    # if an instance has 1 aln it is a unique mapper
    # if an instance has >1 alignment it is multimapper
    instance = defaultdict(list)

    # skip unmapped secondary and supplementary aln (just in case)
    for aln in tqdm(bam_in.fetch(until_eof=True), desc="Reading alignments"):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        instance[aln.query_name].append(aln)

    # Prepare an out bam
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    # Counters for summary
    total_instances = 0        
    unique_instances = 0       
    multimapper_instances = 0  #

    # Assign multimappers (This is the key part of the function)
    # For each read ID resolve multimapping randomly and propagate reads based on RC
    for read_id, alignments in tqdm(instance.items(), desc="Randomly resolve multimappers"):
        total_instances += 1

        # Take the read count per ID (RC)
        rc = _rc_from_read(alignments[0])
        # If unique just propagate RC times
        if len(alignments) == 1:
            unique_instances += 1  # withe unique reads with a singe RC
            _propagate_n_copies(bam_out, _sanitize_for_out(alignments[0]), rc)
        else:
        # Else select a random isntance and then propagare RC times   
            multimapper_instances += 1
            # Take a random read for the MM set
            ran_alignment = random.choice(alignments)
            try:
                ran_alignment.set_tag("ZT", "random")  # Add a tag to mark randoms
            except Exception:
                pass
            # Call the helper to expand
            _propagate_n_copies(bam_out, _sanitize_for_out(ran_alignment), rc)

    # Close bams
    bam_out.close()
    bam_in.close()

    # Sort bams
    if sort_index:
        tmp = output_bam.replace(".bam", ".sorted.bam")
        subprocess.run(["samtools", "sort", "-o", tmp, output_bam], check=True)
        os.replace(tmp, output_bam)
        pysam.index(output_bam)

    # Print report based on the counters for terminal
    rn_time = time.time() - start_time
    print("\nSummary:")
    print(f"- Collapsed instances analyzed: {total_instances}")
    print(f"- Unique instances kept:         {unique_instances}")
    print(f"- Multimappers assigned:         {multimapper_instances}")
    print(f"Extended BAM: {output_bam}")
    print(f"Done in {rn_time:.2f} s.\n")
    
    # Now for log for nextflow
    if summary_log_file:
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "a") as sf:
            sf.write(f"Collapsed instances analyzed: {total_instances}\n")
            sf.write(f"Unique instances kept: {unique_instances}\n")
            sf.write(f"Multimappers assigned: {multimapper_instances}\n")
            sf.write(f"Extended BAM: {output_bam}\n")
            sf.write(f"Random seed: {seed}\n")
            sf.write(f"Run time: {rn_time:.2f}\n")



# Unique Weighted Mode (UWM) similar to ShotStack
def resolve_multimappers_with_index(
    bam_file,
    unique_index_path,
    output_bam,
    window_size=250,
    seed=1,
    sort_index=True,
    summary_log_file=None,
    consider_strand=False,
    ):
    """
    Aim of the function:
    Resolve multimappers using the unique window mode based on index, then expand by RC
    Inputs
    ------------------------------------------------------------------------------------
    bam_file          : BAM with collapsed read IDs
                        (NH and RC tags are expected, since they are added in the mapping step)
    unique_index_path : Index built from unique mappers (NH==1) across
                        selected BAMs using RC weights
    output_bam        : Bam file with resolved multimappers
    window_size       : Span (nt) of the local neighborhood around the alignment center
                        used to score each candidate
    seed              : Random seed to make tie-breaking reproducible
    sort_index        : Run samtools sort and index
    consider_strand   : If True, only count unique support on the same strand;
                        if False (default), ignore strand and use both strands.

    Step by step
    -------------------------------------------------------------------------------------
    For each read ID:
      * If it has exactly 1 possible origin (unique): keep it and write RC copies
      * If it has >1 alignments (multimapper):
        For each candidate alignment, compute a neighborhood score based in unique mappers
        score = sum of RC (from unique mappers) overlapping a window
                centered at the candidate alignment (same chrom & STRAND)
        Select the alignment with the highest score (best local unique context)
        If multiple candidates tie for best score, break ties randomly (this must use a seed)
        Tag chosen alignments with ZT:Z:uwm or random
        Write RC copies of the chosen alignment

    Results:
    -------------------------------------------------------------------------------------
    The resulting BAM has exactly sum(RC over mapped read IDs) records
    Multimappers are placed where uniquely mapped evidence is strongest nearby
    """
    # helper to normalize tags on output records (res)
    def _sanitize_for_out(alignment):
        """
        clean the tags after propagating RC tags in the resolved BAM file
        """
        # Force NH and RC 1 per expanded instance
        try:
            alignment.set_tag("NH", 1, value_type='i')
        except Exception:
            pass

        # Same for RC
        try:
            if alignment.has_tag("RC"):
                alignment.set_tag("RC", 1, value_type='i')
        except Exception:
            pass

        # This is just in case and to standardize the out bams
        # take ride of multimapper and secondary tags 
        try:
            # pysam returns a list of (tag_name, tag_value, value_type)
            all_tags = alignment.get_tags(with_value_type=True)
            tags_to_remove = {"IH", "HI", "XA", "SA"}
            kept_tags = [
                (tag_name, tag_value, value_type)
                for (tag_name, tag_value, value_type) in all_tags
                if tag_name not in tags_to_remove
            ]
            alignment.set_tags(kept_tags)
        except Exception:
            pass

        return alignment


    # Set timer and random seed
    start_time = time.time()
    random.seed(seed)
    # Open the index
    print(f"\nResolving {bam_file} using unique index: {unique_index_path}")
    with open(unique_index_path, "rb") as f:
        uniq_index = pickle.load(f)

    # Also IN and OUT BAMs
    bam_in = pysam.AlignmentFile(bam_file, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    # recycled from rand
    # Group all mapped alignments by read ID and skip non-mappers just in case
    reads = defaultdict(list)
    for r in tqdm(bam_in.fetch(until_eof=True), desc="Grouping reads"):
        if not r.is_unmapped:
            reads[r.query_name].append(r)

    # Divide the window size then center the read
    half_win = window_size // 2 # // make this integer ALWAYS
    # Set counters
    total = 0
    uniq =  0
    multi = 0
    ties = 0 # this new counter in case the read has ties in context

    # This is the most important for in the script
    # score multimappers based on local context (win size)
    for read_id, alns in tqdm(reads.items(), desc="Resolving multimappers"):
        total += 1
        # RC is the same for all alignments of this read ID; grab it from the first
        rc = _rc_from_read(alns[0])

        if len(alns) == 1:
            uniq += 1 # if the rc = 1 just write it into the bam
            _propagate_n_copies(bam_out, _sanitize_for_out(alns[0]), rc)
            continue

        # Now multi mappers
        multi += 1
        scores = []
        for a in alns:
            chrom = a.reference_name # chromosome/counting
            # strand is a key point because I will only consider um on the same strand
            # meaning same biological origin
            center = (a.reference_start + a.reference_end) // 2 # middle of the read
            # read start and end (considering the window size) to ask for overlap and assign score
            s = center - half_win
            e = center + half_win

            if consider_strand:
                strand = '-' if a.is_reverse else '+'
                score = weighted_unique_overlap(uniq_index, chrom, strand, s, e)
            else:
                # ignore strand: sum unique support on both strands
                score = (
                    weighted_unique_overlap(uniq_index, chrom, '+', s, e)
                    + weighted_unique_overlap(uniq_index, chrom, '-', s, e)
                )

            scores.append(score) # append the score after weighted_unique_overlap

        # Choose the candidate with the highest score (best unique support nearby)
        best = max(scores)
        # brake ties randomly
        tied_idxs = [i for i, sc in enumerate(scores) if sc == best] # collects all indices that reach that max
        # If length > 1 random.choice else the first
        if len(tied_idxs) > 1:
            ties += 1
            choice_idx = random.choice(tied_idxs)
        else:
            choice_idx = tied_idxs[0]

        chosen = alns[choice_idx]
        # assign tag for downstream processing
        aln = _sanitize_for_out(chosen)
        try:
            aln.set_tag("ZT", "random" if len(tied_idxs) > 1 else "uwm", value_type='Z')
        except Exception:
            pass
        _propagate_n_copies(bam_out, aln, rc)

    bam_out.close()
    bam_in.close()
    # Sort bams
    if sort_index:
        tmp = output_bam.replace(".bam", ".sorted.bam")
        subprocess.run(["samtools", "sort", "-o", tmp, output_bam], check=True)
        os.replace(tmp, output_bam)
        pysam.index(output_bam)
    # Print summary using counters
    rn_time = time.time() - start_time
    print(f"\nSummary:")
    print(f"- Collapsed sequences processed: {total}")
    print(f"- Unique sequences:             {uniq}")
    print(f"- Multimappers:                 {multi}")
    print(f"- Random ties (uwm):            {ties}")
    print(f"Extended BAM written: {output_bam}")
    print(f"Run time: {rn_time:.2f} s\n")

    # For nextflow I will need a log (base res from rand)
    if summary_log_file:
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "a") as sf:
            sf.write("UWM summary:\n")
            sf.write(f"  Collapsed sequences processed: {total:,}\n")
            sf.write(f"  Unique sequences:             {uniq:,}\n")
            sf.write(f"  Multimappers:                 {multi:,}\n")
            sf.write(f"  Random ties (uwm):            {ties:,}\n")
            sf.write(f"  Random seed:                  {seed:,}\n")
            sf.write(f"  Extended BAM: {output_bam}\n")
            sf.write(f"  Run time: {rn_time:.2f} s\n\n")

# MAIN
def main():
    start_global = time.time()

    parser = argparse.ArgumentParser(
        description="sRNA mapping & multimapper resolution (UWM and Random)"
    )
    parser.add_argument(
        "--mode", required=True,
        choices=["build", "map", "index", "uwm", "random"],
        help=("build: bowtie index | map: collapse FASTQ with bowtie -q BAM + NH/RC tags | "
              "indexing: build weighted unique index from user defined bam-list | "
              "uwm: resolve using unique index (extended BAM) | "
              "random: random assignment (extended BAM)")
    )

    # bt build / map
    parser.add_argument("--genome", help="Genome FASTA file for alignment (Bowtie index)")
    parser.add_argument("--offrate", type=int, default=4, help="bowtie-build offrate (default= 4)")
    parser.add_argument("--fastq", help="FASTQ(.gz) file for mode=map")
    parser.add_argument("--index", help="Bowtie1 index basename for mode=map")
    parser.add_argument("--out", help="Output BAM for mode=map (mapped collapsed reads)")
    parser.add_argument("--mismatches", type=int, default=1, help="Bowtie mismatches (default= 1)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for bowtie/samtools")
    parser.add_argument("--sort-mem", default="2G", help="samtools sort memory per thread (default = 2G)")
    parser.add_argument("--save-np", choices=["yes", "no"], default="no",
                        help="Save unmapped reads via Bowtie's --un/--un-gz (default=no)")
    parser.add_argument("--non-mappers",
                        help="Output fastq(.gz) for unmapped (used when --save-np yes)")

    # collapse temp files (mode=map)
    parser.add_argument("--collapsed-fastq", help="Collapsed fastq file (mode=map) Default: <fastq>.collapsed.fq[.gz]")
    parser.add_argument("--rc-map", help="RC (pickle) map file (mode=map) Default: <fastq>.rc.pkl")

    # index build and resolution
    parser.add_argument("--bam-list", help="Txt file with one BAM path per line")
    parser.add_argument("--unique-index", help="Path to unique index pickle (index output / uwm input)")
    parser.add_argument("--bam", help="Target BAM to resolve (uwm/random)")
    parser.add_argument("--resolved-out", help="Output BAM after resolution (extended BAM)")
    parser.add_argument("--window-size", type=int, default=250, help="Window size for unique-window scoring (default=250)")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for tie-breaks / random mode")

    # Support for stranded umw index. This mode will force to build the index ONLY with um on the same strand
    parser.add_argument(
        "--consider-strand",
        action="store_true",
        help="If provided, use strand-specific unique support when resolving (default: ignore strand)."
    )
    
     
    args = parser.parse_args()
    
    verify_tools()

    if args.mode == "build":
        if not args.genome:
            raise ValueError("--genome is required for mode=build.")
        index_base = os.path.splitext(os.path.basename(args.genome))[0]
        build_bowtie_index(args.genome, index_base, offrate=args.offrate, threads=args.threads)

    elif args.mode == "map":
        if not (args.fastq and args.index):
            raise ValueError("--fastq and --index are required for mode=map.")
        # collapsed outputs
        root = args.fastq[:-3] if args.fastq.endswith(".gz") else args.fastq
        root = os.path.splitext(root)[0]

        # define outputs
        n_collapsed = args.collapsed_fastq or (root + ".collapsed.fastq.gz")
        rc_map_path       = args.rc_map or (root + ".rc.pkl")
        mapped_bam        = args.out or (root + ".bam")

        # log files per sample
        summary_log_file  = mapped_bam.replace(".bam", ".map.summary.log")
        printable_cmd     = f"bowtie -q {args.index} {os.path.basename(n_collapsed)}"

        # write a small header once
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "w") as sf:
            sf.write(f"{time.strftime('%c')}\n")
            sf.write("Command:\n\t" + printable_cmd + " | samtools view -b - | samtools sort ...\n")
            sf.write("Completed.\n\n")

        # collapse fastq gz FASTQ + RC map
        _, _, collapsed_fq = collapse_fastq_to_fastq(
            args.fastq, n_collapsed, rc_map_path, pigz_threads=args.threads,
            summary_log_file=summary_log_file
        )

        # map collapsed FASTQ with Bowtie
        tmp_bam = mapped_bam.replace(".bam", ".tmp.bam")
        run_bowtie(
            fastq_file=collapsed_fq,
            index_base=args.index,
            output_bam=tmp_bam,
            mismatches=args.mismatches,
            save_non_mappers=(args.save_np == "yes"),
            non_mappers_file=args.non_mappers,
            threads=args.threads,
            sort_mem=args.sort_mem,
            index_output_bam=False,
            summary_log_file=summary_log_file,
        )

        # add NH RC using collapse map
        add_NH_RC_tags_from_collapse(tmp_bam, rc_map_path, mapped_bam)
        os.remove(tmp_bam)
        print(f"Mapped (collapsed) BAM with NH/RC: {mapped_bam}")

        # final append
        with open(summary_log_file, "a") as sf:
            sf.write(f"Output BAM: {mapped_bam}\n\n")

    elif args.mode == "index":
        if not args.bam_list:
            raise ValueError("--bam-list is required for mode=index.")
        out_index = args.unique_index or "unique_index.pkl"
        build_weighted_unique_index_from_bams(args.bam_list, out_index)

    elif args.mode == "uwm":
        if not (args.bam and args.unique_index and args.resolved_out):
            raise ValueError("--bam, --unique-index, and --resolved-out are required for mode=uwm.")
        base = os.path.splitext(os.path.basename(args.bam))[0]
        summary_log_file = base + ".map.summary.log"
        resolve_multimappers_with_index(
            bam_file=args.bam,
            unique_index_path=args.unique_index,
            output_bam=args.resolved_out,
            window_size=args.window_size,
            seed=args.seed,
            sort_index=True,
            summary_log_file=summary_log_file,
            consider_strand=args.consider_strand,   # True only if --consider-strand was provided
        )

    elif args.mode == "random":
        if not (args.bam and args.resolved_out):
            raise ValueError("--bam and --resolved-out are required for mode=random.")
        base = os.path.splitext(os.path.basename(args.bam))[0]
        summary_log_file = base + ".map.summary.log"
        assign_multimappers_randomly(
            bam_file=args.bam,
            output_bam=args.resolved_out,
            seed=args.seed,
            sort_index=True,
            summary_log_file=summary_log_file,
        )

    else:
        raise ValueError(f"Unknown mode: {args.mode}")

    rnt_global = time.time() - start_global
    print(f"\nTotal pipeline execution time: {rnt_global:.2f} seconds.\n")


if __name__ == "__main__":
    main()
