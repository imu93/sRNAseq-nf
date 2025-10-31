#!/usr/bin/env python3

# Code writen by: Isaac Martinez Ugalde PhD
# Aim:
# siRmap for sRNA alignment and multimapping placement
# Last update: Oct 29-10-2025
# This is the siRmap v1.0.3 module for mapping small RNA seq data
# The biigest changes are the use of numba to speed up the UWM assignment


# Also some changes to the way the unique index is stored to improve memory usage
# and multimapper assignment speed using the random functions from numba.
# For this siRNmap now creates both a pickle index and a numpy directory with individual
# arrays saved as .npy files for  better handeling of the memory 


# And other improvements to the bowtie mapping step is the addition of the max_multimaps
# parameter to cap the number of reported alignments per read. The idea is to mimic the
# beahvior of ShortStack 3.8.5 where usesrs were able to cap the number of reported alignments


# In addition, have updated the CLI to be more user friendly help messages and defaults 
# This version also adds suppor for saving indices using numpy binary format instead of pickle

# This last feature is under development and may require re-building existing indices
# but is suposed to improve loading times and memory usage acording to numpy docs    

#########################################################################################
# Personal giude for numpy dtypes
# Numeric dtype notes (choose and rationale)
# np.int32 : 32-bit signed integer (4 bytes). Range: -2_147_483_648 .. 2_147_483_647.
#     Good for genomic coordinates (start/end positions) 
#     Use when you have millions of records and positions
#########################################################################################
# np.int64 : 64-bit signed integer (8 bytes). Range: very large (~9e18).
#     Use for cumulative counters, prefix sums and any aggregated counts (RC) because
#       sums over many reads can exceed int32 limits. Also safer for total read counts
#       across many libraries or the whole dataset
#########################################################################################
# np.int8 / np.uint8 : 8-bit integers (1 byte). Use for compact small counters or flags
# strand: store as -1/+1 in int8 (or 0/1 in uint8) to save memory
#########################################################################################
# object dtype : used for strings (chrom names). This is flexible but
# memory-inefficient and slower for numeric operations 
#########################################################################################
# float64 / float32 : floating point types used for probabilities
#########################################################################################
# Guidelines for siRmap index storage:
#  starts/ends stored as np.int32 to save memory and because chromosome coordinates fit
#  per-read RC stored and prefix-sums computed as np.int64 (safe for summed counts)
#  strands stored as np.int8
#  chrom names are kept as object dtype but consider producing a chrom-code (np.int32)
#  and using that for sorting/lexsort to speed operations and reduce memory
# In case I need more info: https://numpy.org/doc/stable/user/basics.types.html

import shutil
import os
import time
import random
import hashlib
import gzip
import argparse
import pysam
import pickle
import numpy as np
from collections import defaultdict
from numba import njit, prange
import sys
from tqdm import tqdm
import subprocess
from collections import defaultdict, Counter
import fnmatch
sys.dont_write_bytecode = True


# Tool verification
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


# I will use pigz to zip files
def _pigz_compress(
    in_file,
    out_file,
    threads=0
    ):
    """
    Wrapper for pigz
    """
    # Check outs are gz just in case
    if not out_file.endswith(".gz"):
        raise ValueError("out_file must end with .gz")
    # open infile and out file  
    with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
        cmd = ["pigz", "-c"] # base pigz line
        # if to check the numebr of threads
        if threads > 0:
            cmd = ["pigz", "-p", str(threads), "-c"]
        # Call pigz
        proc = subprocess.Popen(cmd, stdin=fin, stdout=fout)
        # wait till pigz ends and capture rc as a control
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"pigz failed compressing {in_file} to {out_file}")
    os.remove(in_file)

# Support function to parse the RC map from collapse, no more dependency on siRmap for this task
# The C collapse script must be used to produce the map
def _load_rc_map(path):
    """
    Load RC map produced with collapse module.
    """
    # Some ifs for control
    if not path:
        return None

    # Allow multiple trailing .gz extensions in filename (tolerant to pipeline naming mistakes)
    is_gz = False
    test_name = path
    while test_name.endswith(".gz"):
        test_name = test_name[:-3]
        is_gz = True

    # Require underlying file to be TSV
    if not test_name.endswith(".tsv"):
        raise ValueError("This is not a valid MAP input. Expected .tsv or .tsv(.gz).")

    opener = gzip.open if is_gz else open

    rc = {}
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) != 2:
                    # skip malformed lines silently (could be logged)
                    continue
                rid, count = cols[0], cols[1]
                try:
                    rc[rid] = int(count)
                except ValueError:
                    # skip non-integer counts
                    continue
    except FileNotFoundError:
        raise
    return rc

def parse_read_id_for_rc(read_name, rc_map):
    if rc_map is None:
        return 1
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
# I have modified this function to allow capping the number of multimappers reported per read
# And also to remove the sorting step from the bowtie pipeline, this will be done later if needed in the asignment step
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
    max_multimaps=None,
):
    """
    Map reads with Bowtie1 and produce a sorted BAM.
    Optionally save unmapped reads using Bowtie
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

    # Split threads across stages to avoid issues with local nextflow
    T = max(1, int(threads))
    view_threads  = 1
    sort_threads  = max(1, T // 3)
    bt_threads    = max(1, T - view_threads - sort_threads)

    # Build Bowtie command (keep your original -v mismatches behavior)
    cmd = ["bowtie", "-S", "-q", "-v", str(mismatches)]

    # Reporting mode: all vs capped
    if max_multimaps is None:
        cmd += ["-a"]                       # report all alignments (existing default)
    else:
        try:
            n = int(max_multimaps)
            if n <= 0:
                cmd += ["-a"]
            else:
                cmd += ["-k", str(n)]      # report up to n alignments per read
        except Exception:
            cmd += ["-a"]

    cmd += ["--best", "--strata", "-p", str(bt_threads)]
    # Save unmaped
    tmp_un_plain = None
    unmapped_gz_path = None
    if save_non_mappers and non_mappers_file:
        # ensure directory exists
        # Validate directory, if not build
        os.makedirs(os.path.dirname(os.path.abspath(non_mappers_file)) or ".", exist_ok=True)
        # Now let's force the non mappers as gz to save space
        unmapped_gz_path = non_mappers_file if non_mappers_file.endswith(".gz") else (non_mappers_file + ".gz")
        tmp_un_plain = unmapped_gz_path[:-3]  # write plain FASTQ first
        cmd += ["--un", tmp_un_plain]

    # Add index and reads
    cmd += [index_base, fastq_file]

    # Note: I have had issues with samtools before when using ShortStack 3.8.5 because of the sort mem thus:
    sort_mem_value = sort_mem.upper()
    if sort_mem_value.endswith("G"):
        total_mem_bytes = int(float(sort_mem_value[:-1]) * 1024 * 1024 * 1024)
    elif sort_mem_value.endswith("M"):
        total_mem_bytes = int(float(sort_mem_value[:-1]) * 1024 * 1024)
    else:
        total_mem_bytes = int(sort_mem_value)
    # per-thread mem for samtools sort
    
    # Compute per-thread mem using sort_threads (not total threads)
    sort_mem_bytes = total_mem_bytes // max(1, sort_threads)

    printable_cmd = " ".join(cmd)
    print("\nRunning Bowtie mapping:")
    # Change: remove sort, just pipe to view
    print(printable_cmd + f" | samtools view -b -@ {view_threads} - -o {output_bam}\n")

    with open(raw_log_file, "w") as lf:
        lf.write(f"Start time: {time.strftime('%c')}\n")
        lf.write("Running command: " + printable_cmd + "\n")

        # Change: remove sort from pipeline
        bowtie = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        view   = subprocess.Popen(["samtools", "view", "-b", "-@", str(view_threads), "-", "-o", output_bam],
                                stdin=bowtie.stdout)

        # Make sure pipe isn't blocking
        bowtie.stdout.close()

        # Capture bowtie stderr into the log in real time
        for line in iter(bowtie.stderr.readline, b""):
            s = line.decode(errors="replace")
            lf.write(s)
        bowtie.stderr.close()

        # Wait for the pipeline
        rc_view = view.wait()
        rc_bt   = bowtie.wait()

        lf.write(f"\nEnd time: {time.strftime('%c')}\n")
        lf.flush()

    # Check return codes
    if rc_bt != 0 or rc_view != 0:
        raise RuntimeError(f"bowtie/view pipeline failed (bowtie={rc_bt}, view={rc_view}). See {raw_log_file}")

    print(f"BAM ready: {output_bam}")

    if tmp_un_plain is not None:
        print(f"Compressing unmapped reads with pigz → {unmapped_gz_path}")
        _pigz_compress(tmp_un_plain, unmapped_gz_path, threads=threads)

     # bam index ?
    if index_output_bam:
        try:
            pysam.index(output_bam)
        except Exception:
            pass

    # Summarize alignment 
    print("\nGenerating alignment summary")
    bam_in = pysam.AlignmentFile(output_bam, "rb")

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



# Add tags based on RC files but also NH
def add_NH_RC_tags_from_collapse(
    input_bam,
    rc_map_path,
    output_bam
    ):
    """
    Add NH (alignments per collapsed read id) and RC using the saved collapse map.
    If rc_map_path is None, RC will be set to 1 for all mapped reads.
    """
    print(f"\nAdding NH/RC (from collapse map) to {input_bam} ...")
    rc_map = _load_rc_map(rc_map_path) if rc_map_path else None

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
            rc = parse_read_id_for_rc(read.query_name, rc_map) # read RC (or 1 if no map)
            read.set_tag("RC", int(rc)) # add tags
        bam_out.write(read)

    bam_in.close()
    bam_out.close()
    print(f"NH/RC tags added and saved to {output_bam}\n")



# I will build an support function to read my list of bams
def _read_bam_paths(bam_list_or_bam):
    # If I have a single bam for the unique index not usual because
    # I have replicates per library (debug)
    if bam_list_or_bam.endswith(".bam"):
        return [bam_list_or_bam]
    # Now save paths a list
    paths = []
    # Open the txt
    with open(bam_list_or_bam, "rt") as file:
        for line in file:
            line = line.strip()            # remove jumps of line \n
            if not line or line.startswith("#"):
                continue
            paths.append(line)
    return paths


# This function will save the index in numpy format
# This is under development so that I can improve loading times and memory usage
def _save_index_numpy_dir(index: dict, out_prefix: str):
    """
    Save index as a directory of .npy files, one file per array
    Directory name will be prefix.npidx unless out_prefix already ends with .npidx
    """
    # choose directory name
    out_dir = out_prefix if out_prefix.endswith(".npidx") else os.path.splitext(out_prefix)[0] + ".npidx"
    # if a file with that name exists, append suffix to avoid error
    if os.path.exists(out_dir) and not os.path.isdir(out_dir):
        out_dir = out_dir + ".d"
    os.makedirs(out_dir, exist_ok=True)
    
    # Save each block as separate .npy files
    for (chrom, strand), block in index.items():
        # sanitize chrom to avoid path separators in filenames
        safe_chrom = str(chrom).replace(os.sep, "_").replace(":", "_")
        keyname = f"{safe_chrom}__{strand}"
        # each block has four arrays: starts_pos, starts_pref, ends_pos, ends_pref
        np.save(os.path.join(out_dir, f"{keyname}__starts_pos.npy"),  block["starts_pos"])
        np.save(os.path.join(out_dir, f"{keyname}__starts_pref.npy"), block["starts_pref"])
        np.save(os.path.join(out_dir, f"{keyname}__ends_pos.npy"),    block["ends_pos"])
        np.save(os.path.join(out_dir, f"{keyname}__ends_pref.npy"),   block["ends_pref"])
    return out_dir


# Now I need to tweak my function to read multiple bams 
def build_weighted_unique_index_from_bams(bam_list_or_bam, output_index="unique_index.pkl"):
    # Lets build lists to save the meta data
    chr_list    = []
    start_list  = []
    end_list    = []
    strand_list = []
    rc_list     = []
    
    # Read my list of bams
    bam_paths = _read_bam_paths(bam_list_or_bam)
    if not bam_paths:
        raise ValueError("ERROR: No BAMs found, please pass a valid bam list")
    # Now per library I will read each alignment
    for lib in bam_paths:
        with pysam.AlignmentFile(lib, "rb") as bam:
            # iterate reads with an inner progress bar (no total available)
            for read in tqdm(bam.fetch(until_eof=True), desc=f"Reading {os.path.basename(lib)}", unit="reads", leave=False):
                # Skip unmapped reads (just in case)
                if read.is_unmapped:
                    continue
                nh = read.get_tag("NH")
                # if not unique skip
                if nh != 1:
                    continue
            
                # And the read count per collapsed read
                rc_list.append(int(read.get_tag("RC")))
                    
                # Auxiliar flags for the index
                chr_list.append(read.reference_name)
                start_list.append(read.reference_start)
                end_list.append(read.reference_end)
                strand_list.append(-1 if read.is_reverse else 1)

 
    # Now I'll build indexes in numpy
    chroms = np.array(chr_list,    dtype=object)
    starts = np.array(start_list,  dtype=np.int32)
    ends   = np.array(end_list,    dtype=np.int32)
    strands= np.array(strand_list, dtype=np.int8)
    rc     = np.array(rc_list,     dtype=np.int64)
    # Set a mask just to save valid reads 
    mask = (ends > starts)
    # I'll apply the same mask for all the vectors
    chroms  = chroms[mask]
    starts  = starts[mask]
    ends    = ends[mask]
    strands = strands[mask]
    rc      = rc[mask]
    # Now the idea is to have an ordered array so for this I need to order
    # By chromosome start and strand
    # To asses this I will use a "factor-like"
    # In python, _, helps to ignore the first output,
    # which in this case are the "levels" of a factor in R
    _, chrom_code = np.unique(chroms, return_inverse=True)
    order = np.lexsort((strands, starts, chrom_code))
    # Let's use this order to re-index 
    chroms_sorted  = chroms[order]
    starts_sorted  = starts[order]
    ends_sorted    = ends[order]
    strands_sorted = strands[order]
    rc_sorted      = rc[order]
    
    # Let's define the chr that can play
    index = {} 
    chroms_unique = np.unique(chroms_sorted)
    for chrom in chroms_unique:
        for strand_i in (-1, +1):   # it could be + or -
            # Crate a new mask
            mask_ch_str = (chroms_sorted == chrom) & (strands_sorted == strand_i)
            
            # Take only "CONTIGUOUS" chucks of this Chr and strand
            start_block  = starts_sorted[mask_ch_str]  # Block start
            rcs_block    = rc_sorted[mask_ch_str]      # RC for starts in the block
            end_block    = ends_sorted[mask_ch_str]    # Ends of the block        


            # Now is time to do a trick with the windows to precompute the RC
            # This potentially can help to save time when counting overlaps around a 
            # window defined by each multimapper. 
            
            # Here I'll build to lists of cumulative sum of RC one for start and one of the end
            # The idea is to say how many read counts I have before the end - the reads I have 
            # before start. So in other worlds lets imagine I have a read with to parts R (end of the read)
            # and L (left of the read. How many RC I have within a defined window based on L and R
            
            # If the read is in this direction ---- L ---------> R ----
            # <= R - < L
            
            # So lets use the index of cumsums this will be L
            starts_pos  = start_block.astype(np.int32, copy=False)
            starts_pref = np.cumsum(rcs_block, dtype=np.int64)

            if end_block.size:
                # Now the same for the end. But first sort
                ord_end     = np.argsort(end_block, kind="stable")
                # Reorder the end and now as init32
                ends_pos  = end_block[ord_end].astype(np.int32, copy=False)
                # and estimate the second list of cumsums so R
                ends_pref = np.cumsum(rcs_block[ord_end], dtype=np.int64)
                
            else:
                # if there are no reads in the block produce empty arrays
                ends_pos  = np.empty(0, dtype=np.int32)
                ends_pref = np.empty(0, dtype=np.int64)
                
            # Now add strand and close the index
            strand_ch = '-' if strand_i < 0 else '+'
            index[(chrom, strand_ch)] = {
                "starts_pos":  starts_pos,
                "starts_pref": starts_pref,
                "ends_pos":    ends_pos,
                "ends_pref":   ends_pref,
            }
    # Save the index using pickle, like an RDS in R
    with open(output_index, "wb") as f:
        pickle.dump(index, f)
    print(f"Weighted unique index saved to {output_index}")
    
    
    # Under development: save numpy-backed index
    # Also save a numpy-backed directory .npidx for fast mmap-able loads
    try:
        npidx_dir = _save_index_numpy_dir(index, output_index)
        print(f"Also saved numpy index dir to {npidx_dir}")
    except Exception as e:
        print(f"Warning: failed to write .npidx index: {e}")
    
#######################################################################################################
#   Aiming to speed up assigmnet of uwm I will use numba, which is supossed to rich C speed levels
#   The logic is this:
#   For a window [L, R] and a (chrom, strand) block of UNIQUE reads,
#   weighted_overlaps = (RC start < R) - (RC end ≤ L)
#
#   The "left" search finds how many starts are strictly BEFORE R.
#   The "right" search finds how many ends are AT OR BEFORE L.
#######################################################################################################
# This first function will help to get all the elements of the array
# before R limit
@njit(cache=False, fastmath=True)
def leftmost_start_at_or_after(positions_sorted, right_edge_R):
    """
    For half-open windows -----L----|---->R----
                                 included   excluded
    Get all index positions before R
    """
    # The index always [left_bound, right_bound),
    # where right_bound is exclusive
    left_bound = 0
    # R = the length of the array of sorted positions
    right_bound = positions_sorted.size
    
    while left_bound < right_bound:
        # Choose the midpoint of the current candidate range.
        # -----L----|---->R----
        midpoint = (left_bound + right_bound) // 2
        #   value at midpoint is LEFT of the bar (value < R)
        #   ----L----|----R----
        #      … < … |----R----
        if positions_sorted[midpoint] < right_edge_R:
            left_bound = midpoint + 1
        else:
            # Else is R
            right_bound = midpoint
    # left_bound now points at the first position >= query_pos.
    return left_bound

# All elements AFTER L, so we can later compute counts inside [L, R)
@njit(cache=False, fastmath=True)
def leftmost_end_after(positions_sorted, left_edge_L):
    left_bound = 0
    right_bound = positions_sorted.size
    while left_bound < right_bound:
        midpoint = (left_bound + right_bound) // 2
        if positions_sorted[midpoint] <= left_edge_L:
            left_bound = midpoint + 1
        else:
            right_bound = midpoint
    return left_bound


# Single-window overlap using prefix sums
@njit(cache=False, fastmath=True)  # Added cache and fastmath
def weighted_overlap_in_window(
    left_edge_L,
    right_edge_R,
    start_coords,      # sorted starts
    start_prefix_rc,   # prefix sum of RC aligned to start_coords
    end_coords,        # sorted ends
    end_prefix_rc      # prefix sum of RC aligned to end_coords
):
    """
    Half-open window   -----L----|---->R----
                               included   excluded

    Compute the RC-weighted overlap of UNIQUE reads with [L, R):
      overlap = (RC of starts where start < R) - (RC of ends where end <= L)
    """
    if right_edge_R <= left_edge_L:
        return np.int64(0)

    # Use faster binary search for both operations
    idx_first_start_ge_R = np.searchsorted(start_coords, right_edge_R, side='left')
    idx_first_end_gt_L = np.searchsorted(end_coords, left_edge_L, side='right')

    # Optimized RC calculations with bounds checking
    rc_starts_before_R = start_prefix_rc[idx_first_start_ge_R - 1] if idx_first_start_ge_R > 0 else np.int64(0)
    rc_ends_at_or_before_L = end_prefix_rc[idx_first_end_gt_L - 1] if idx_first_end_gt_L > 0 else np.int64(0)

    return rc_starts_before_R - rc_ends_at_or_before_L


# Many windows at once for one (chrom,strand) block
@njit(cache=False, fastmath=True)
def weighted_overlaps_batch(
    left_edges_L,       # array of L
    right_edges_R,      # array of R 
    start_coords,       # sorted starts
    start_prefix_rc,    # RC for starts
    end_coords,         # sorted ends
    end_prefix_rc       # RC for ends
):
    # How many windows we have 
    n_windows = left_edges_L.size
    
    # I need an "empthy" array collect the index
    overlaps_rc = np.zeros(n_windows, dtype=np.int64)
    for window_idx in range(n_windows):
        # and now call the single index function to collect the RCs one by one
        overlaps_rc[window_idx] = weighted_overlap_in_window(
            left_edges_L[window_idx],
            right_edges_R[window_idx],
            start_coords,
            start_prefix_rc,
            end_coords,
            end_prefix_rc,
        )

    return overlaps_rc

########################################################################################################
# I will use this function to build the windows arround the multiple loci per multimapper
def score_windows_block(index, chrom, strand, centers, width):
    """
    Build [L, R) windows around "centers" and score them against the UNIQUE index block
    """
    # Let's use the index of chrom and strand to retrieve the ordered array
    block = index.get((chrom, strand))
    # prevent worng of empty indices
    if block is None or block["starts_pos"].size == 0:
        return np.zeros_like(centers, dtype=np.int64)
        
    # turn centers into window edges [L, R)
    half = width // 2
    Ls = (centers - half).astype(np.int64, copy=False)
    Rs = (centers + half).astype(np.int64, copy=False)
    # ensure contiguous arrays for numba
    start_pos_bp    = np.ascontiguousarray(block["starts_pos"])
    start_prefix_rc = np.ascontiguousarray(block["starts_pref"])
    end_pos_bp      = np.ascontiguousarray(block["ends_pos"])
    end_prefix_rc   = np.ascontiguousarray(block["ends_pref"])
    # and return the array
    return weighted_overlaps_batch(Ls, Rs, start_pos_bp, start_prefix_rc, end_pos_bp, end_prefix_rc)

# Finally the main function to choose the best locus based on UNIQUE context
# Here I will basiclly reuse the logic from siRmap but adapted to use the numba functions
def choose_candidate_with_index(index, read_id, candidates, window_size, consider_strand, seed=1):
    """
    Decide the most likely locus for a multimapping read using UNIQUE support
    """
    # For random placement of reads use a fixed seed
    rnd = random.Random(seed + _stable_int_from_str(read_id))
    # Number of candidate loci
    n = len(candidates)
    if n == 1:
        return (0, False, 1.0)
    # Precompute per-candidate centers, fist i need some empty arrays
    # This will be like empty matrices in R
    centers = np.empty(n, dtype=np.int64)
    chroms  = []
    strands = []
    # Now for the possible loci
    for i, (chrom, start, end, is_rev) in enumerate(candidates):
        # save the center an, chr and strand
        centers[i] = (start + end) // 2
        chroms.append(chrom)
        strands.append('-' if is_rev else '+')
    # No empth array for scores
    scores = np.zeros(n, dtype=np.int64)
    # I would like to keep the possibility of consider strand so if consider strand true
    # I will group by (chrom, strand) otherwise just by chrom like in ShortStack
    if consider_strand:
        # Group by (chrom, strand) and score each group in one batch call
        groups = defaultdict(list)
        for i in range(n):
            groups[(chroms[i], strands[i])].append(i)
        # Call score index block and pass the centers
        # Fill the score array
        for key, idxs in groups.items():
            idxs_arr         = np.array(idxs, dtype=np.int64) # call the index array
            sub_centers      = centers[idxs_arr] # I have my centers now
            sub_scores       = score_windows_block(index, key[0], key[1], sub_centers, window_size) # collect scores
            scores[idxs_arr] = sub_scores # and pass to the array of scores
    else:
        # Strand-agnostic like in ShortStack
        by_chrom = defaultdict(list)
        for i in range(n):
            # Group by chrom only and append index to thet by_chrom list
            by_chrom[chroms[i]].append(i)
        # Now per chrom call both strands and sum
        for chrom, idxs in by_chrom.items():
            idxs_arr    = np.array(idxs, dtype=np.int64)
            sub_centers = centers[idxs_arr]
            # So get the scores per strand first
            strand_plus  = score_windows_block(index, chrom, '+', sub_centers, window_size)
            strand_minus = score_windows_block(index, chrom, '-', sub_centers, window_size)
            # Sum the estimated scores and sum both strands
            scores[idxs_arr] = strand_plus + strand_minus

    # Convert scores to probabilities
    # if negative scores (should not happen) set to zero 
    scores = np.maximum(scores, 0)
    total = float(scores.sum())
    # Normalize to get probabilities that sum 1
    if total > 0.0:
        probs = (scores / total).astype(np.float64, copy=False)
    else:
        probs = np.full(n, 1.0 / n, dtype=np.float64)

    # Pick best; break ties deterministically with seed
    best = probs.max()
    tied = [i for i, p in enumerate(probs) if p == best]
    had_tie = (len(tied) > 1)
    best_idx = rnd.choice(tied) if had_tie else tied[0]
    return (best_idx, had_tie, float(probs[best_idx]))


# Now I will reuse sanitize function from my original siRmap the basic idea of the function
# is to expand the number or reads per RC this function is essential for downstream analysis 
# like featureCounts, because I must keep the right number of reads after collapsing and 
# align the reads
def _sanitize_for_out_inplace(alignment):
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
    
# I will use hashes to properly name read (avoid redundancy)
def _stable_int_from_str(s: str) -> int:
    return int.from_bytes(hashlib.md5(s.encode('utf-8')).digest()[:8], 'little', signed=False)

# I will add a support function for RC counting
def _get_rc_safe(read: pysam.AlignedSegment) -> int:
    """Collapsed count if present, else 1"""
    try:
        return int(read.get_tag("RC"))
    except Exception:
        return 1

def _clear_secondary_and_supplementary(aln: pysam.AlignedSegment):
    """Make sure we output only primary-like alignments"""
    # 0x100 = secondary, 0x800 = supplementary
    aln.flag &= ~0x100
    aln.flag &= ~0x800
    # Make sure it's marked as mapped
    aln.flag &= ~0x4

def _expand_write(read: pysam.AlignedSegment, out_bam: pysam.AlignmentFile, rc: int,  zt_value: str | None = None):
    """
    Write 'rc' copies of aln to out_bam after sanitizing
    Each copy has NH=1 and RC=1 and no multimapper/supplementary tags
    """
    # sanitize in-place (okay to mutate a temp object)
    _clear_secondary_and_supplementary(read)
    # add ZT tag whcih can be used to track unique/multimapper status
    # also I will use it to  track probabilities associated with UWM
    # when radom, the probability is 1.0
    if zt_value is not None:
        try:
            read.set_tag("ZT", zt_value, value_type='Z')
        except Exception:
            pass
          
    _sanitize_for_out_inplace(read)
    for _ in range(int(rc)):
        out_bam.write(read)


# Group alignments by query name to process multimappers together
def _group_alignments_by_qname(in_bam_path: str):
    # defaultdict of lists
    groups = defaultdict(list)
    # Read the BAM and group alignments
    with pysam.AlignmentFile(in_bam_path, "rb") as bam:
        header = bam.header  # take header here
        for read in bam.fetch(until_eof=True):
            # Skip unmapped and supplementary alignments
            if read.is_unmapped or read.is_supplementary:
                continue
            # Group by query name
            groups[read.query_name].append(read)
    return groups, header


# This function will load the numpy index saved with _save_index_numpy_dir
# I will use memory mapping to keep memory low
def _load_index_numpy_dir(index_dir: str, use_mmap: bool = True):
    """
    Load an index saved with _save_index_numpy_dir
    """
    if not os.path.isdir(index_dir):
        raise FileNotFoundError(index_dir)
    index = {}
    # Check the endings to identify blocks of arrays
    for fname in os.listdir(index_dir):
        # lets find only the starts_pos files to identify blocks
        if not fname.endswith("__starts_pos.npy"):
            continue
        keyname = fname[:-len("__starts_pos.npy")]
        # reconstruct filenames
        f_starts_pos  = os.path.join(index_dir, f"{keyname}__starts_pos.npy")
        f_starts_pref = os.path.join(index_dir, f"{keyname}__starts_pref.npy")
        f_ends_pos    = os.path.join(index_dir, f"{keyname}__ends_pos.npy")
        f_ends_pref   = os.path.join(index_dir, f"{keyname}__ends_pref.npy")
        mmap_mode     = "r" if use_mmap else None
        starts_pos    = np.load(f_starts_pos, mmap_mode=mmap_mode)
        starts_pref   = np.load(f_starts_pref, mmap_mode=mmap_mode)
        ends_pos      = np.load(f_ends_pos, mmap_mode=mmap_mode)
        ends_pref     = np.load(f_ends_pref, mmap_mode=mmap_mode)
        # recover chrom and strand (split on last '__' to allow chrom names with '__')
        parts = keyname.rsplit("__", 1)
        if len(parts) == 2:
            chrom, strand = parts
        else:
            chrom, strand = keyname, "+"
        chrom = chrom.replace("_", "/") if "_" in chrom else chrom  # best-effort reverse of sanitize
        index[(chrom, strand)] = {
            "starts_pos":  starts_pos,
            "starts_pref": starts_pref,
            "ends_pos":    ends_pos,
            "ends_pref":   ends_pref,
        }
    return index

# Since I remove the samtools sort from the mapping function I need my memory helper
# to parse human readable memory strings like 2G, 500M, 100K
def _parse_mem_to_bytes(s: str) -> int:
    s = s.strip().upper()
    if s.endswith("G"):
        return int(float(s[:-1]) * 1024 * 1024 * 1024)
    if s.endswith("M"):
        return int(float(s[:-1]) * 1024 * 1024)
    if s.endswith("K"):
        return int(float(s[:-1]) * 1024)
    return int(s)
# Now I need a new function to sort and index bams using samtools
def sort_and_index_bam(
    bam_path: str,
    threads: int = 2,
    sort_mem: str = "2G",
    inplace: bool = True,
):
    """
    Modular samrtools sortfunction  
    """
    tmp_sorted = bam_path if not inplace else bam_path.replace(".bam", ".sorted.bam")
    # be conservative with threads
    sort_threads = max(1, min(int(threads), 4))
    total_bytes = _parse_mem_to_bytes(sort_mem)
    # keep per-thread memory modest to avoid OOMs
    mem_per_thread_mb = max(100, (total_bytes // max(1, sort_threads * 4)) // (1024 * 1024))

    subprocess.run(
        [
            "samtools", "sort",
            "-@", str(sort_threads),
            "-m", f"{mem_per_thread_mb}M",
            "-o", tmp_sorted,
            bam_path,
        ],
        check=True
    )
    if inplace:
        os.replace(tmp_sorted, bam_path)
    # Now use try except to avoid crashes, not critical but ill keep this structure
    try:
        pysam.index(bam_path if inplace else tmp_sorted)
    except Exception:
        pass


# Unique Window Mode UWM
def resolve_multimappers_uwm(
    index_path: str,
    query_bam_path: str,
    output_bam_path: str,
    window_size: int = 100,
    consider_strand: bool = False,
    seed: int = 1,
    progress: bool = True,
    summary_log_file: str | None = None,
    do_sort: bool = True,
    sort_threads: int = 2,
    sort_mem: str = "2G",
):
    start_time = time.time()
    """
    Resolve multimapping reads using Unique Window Mode (UWM) based on a weighted unique index
    """

    # Load index (support both old .pkl and new .npidx directory)
    if os.path.isdir(index_path) or index_path.endswith(".npidx"):
        index = _load_index_numpy_dir(index_path, use_mmap=True)
    else:
        # fallback to pickle for backward compatibility
        with open(index_path, "rb") as f:
            index = pickle.load(f)

    # Call my function to group alignments by read 
    groups, header = _group_alignments_by_qname(query_bam_path)
    n_groups = len(groups)

    # Counters for progress and report
    total_instances = 0
    unique_instances = 0
    multimapper_instances = 0
    tie_breaks = 0
    prob_sum = 0.0  
    total_reads_written = 0
    unique_reads_written = 0
    mm_reads_written = 0

    # I'll use pysam to open the query bam
    # Need to be sure that the out_dir exist (this comes from siRmap)
    os.makedirs(os.path.dirname(output_bam_path) or ".", exist_ok=True)

    # Open out bam
    with pysam.AlignmentFile(output_bam_path, "wb", header=header) as out_bam:
        # I want a nice progress bar like in the random resolver
        iterator = groups.items()
        if progress:
            iterator = tqdm(
                iterator,
                total=n_groups,
                desc="Resolving (UWM)",
                unit="readID",
                smoothing=0.05,
            )

        # Go locus by locus for the possible sites
        for read_id, reads in iterator:
            total_instances += 1

            # Filter to candidates
            candidates = []
            rc_from_any = 1
            for loci in reads:
                # record RC (if multiple, they *should* be same for the group)
                rc_from_any = _get_rc_safe(loci) if rc_from_any == 1 else rc_from_any
                chrom = loci.reference_name
                r_start = loci.reference_start
                r_end = loci.reference_end
                is_rev = bool(loci.is_reverse)
                candidates.append((chrom, r_start, r_end, is_rev))

            if not candidates:
                # nothing to do for this read
                continue

            # Decide if this read is unique:
            # Seek for only one candidate record OR NH tag == 1 on at least one candidate
            is_unique_by_count = (len(candidates) == 1)
            is_unique_by_tag = False
            for aln in reads:
                try:
                    if aln.get_tag("NH") == 1:
                        is_unique_by_tag = True
                        break
                except Exception:
                    pass

            if is_unique_by_count or is_unique_by_tag:
                # choose the primary if present, else first
                primaries = [a for a in reads if not a.is_secondary]
                chosen = primaries[0] if primaries else reads[0]
                rc_from_any = _get_rc_safe(chosen)
                # mark as unique (keep my ZT style)
                _expand_write(chosen, out_bam, rc_from_any, zt_value="unique:1")

                unique_instances += 1
                unique_reads_written += int(rc_from_any)
                total_reads_written += int(rc_from_any)

                # keep tqdm snappy but not spammy
                if progress and (total_instances % 1000 == 0 or total_instances == n_groups):
                    avgp = (prob_sum / multimapper_instances) if multimapper_instances else 0.0
                    iterator.set_postfix(
                        uniques=unique_instances,
                        multis=multimapper_instances,
                        ties=tie_breaks,
                        avgp=f"{avgp:.3f}",
                    )
                continue

            # Resolve using the UW index
            # choose_candidate_with_index returns (best_idx, had_tie, prob_best)
            best_idx, had_tie, prob = choose_candidate_with_index(
                index=index,
                read_id=read_id,
                candidates=candidates,
                window_size=window_size,
                consider_strand=consider_strand,
                seed=seed,
            )
            # add counters for report
            multimapper_instances += 1
            mm_reads_written += int(rc_from_any)
            total_reads_written += int(rc_from_any)
            if had_tie:
                tie_breaks += 1
            prob_sum += float(prob)

            # pick the corresponding pysam alignment object
            chosen_aln = reads[best_idx]
            rc_from_any = _get_rc_safe(chosen_aln)

            # >>> ZT tag for multimappers (uwm vs random)
            zt = f"{'random' if had_tie else 'uwm'}:{prob:.6f}"
            _expand_write(chosen_aln, out_bam, rc_from_any, zt_value=zt)

            # update the progress bar periodically
            if progress and (total_instances % 1000 == 0 or total_instances == n_groups):
                avgp = (prob_sum / multimapper_instances) if multimapper_instances else 0.0
                iterator.set_postfix(
                    uniques=unique_instances,
                    multis=multimapper_instances,
                    ties=tie_breaks,
                    avgp=f"{avgp:.3f}",
                )
    if do_sort:
        print("Sorting and indexing resolved BAM...")
        sort_and_index_bam(output_bam_path, threads=sort_threads, sort_mem=sort_mem, inplace=True)

    # Compute summary stats once we are done
    run_time = time.time() - start_time
    pct_uni = (unique_instances / max(1, total_instances)) * 100.0
    pct_mm = (multimapper_instances / max(1, total_instances)) * 100.0
    avgp = (prob_sum / max(1, multimapper_instances))

    # Console summary (keep my style)
    print("\nSummary (UWM):")
    print(f" Collapsed instances analyzed: {total_instances:,}")
    print(f" Unique instances kept:        {unique_instances:,} ({pct_uni:.1f}%)")
    print(f" Multimappers resolved:        {multimapper_instances:,} ({pct_mm:.1f}%)")
    print(f" Ties broken:                  {tie_breaks:,}")
    print(f" Mean UWM confidence:          {avgp:.3f}")
    print(f" Output BAM: {output_bam_path}")
    print(f" Window={window_size}, Strand-aware={consider_strand}, Seed={seed}")
    print(f" Reads written (expanded):      {total_reads_written:,}")
    print(f"   from uniques:              {unique_reads_written:,}")
    print(f"   from multimappers:         {mm_reads_written:,}")
    print(f"Done in {run_time:.2f} s.\n")

    # Now save the same info to a log for nextflow
    if summary_log_file:
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "a") as sf:
            sf.write("=== UWM Resolution ===\n")
            sf.write(f"Input BAM: {query_bam_path}\n")
            sf.write(f"Index: {index_path}\n")
            sf.write(f"Window: {window_size}\n")
            sf.write(f"Strand-aware: {consider_strand}\n")
            sf.write(f"Seed: {seed}\n")
            sf.write(f"Collapsed instances analyzed: {total_instances}\n")
            sf.write(f"Unique instances kept: {unique_instances} ({pct_uni:.1f}%)\n")
            sf.write(f"Multimappers resolved: {multimapper_instances} ({pct_mm:.1f}%)\n")
            sf.write(f"Ties broken: {tie_breaks}\n")
            sf.write(f"Mean UWM confidence: {avgp:.3f}\n")
            sf.write(f"Output BAM: {output_bam_path}\n")
            sf.write(f"Run time: {run_time:.2f} s\n")
            sf.write(f"Reads written (expanded): {total_reads_written:,}\n")
            sf.write(f"Expanded unique: {unique_reads_written:,}\n")
            sf.write(f"Expanded multimappers: {mm_reads_written:,}\n\n")

    return output_bam_path


# Support function for numpy/numba random assignment
@njit(cache=False, fastmath=True)
def _random_assign_numba(n_aligns):
    """Fast random pick using numba"""
    return np.random.randint(0, n_aligns)


# Radom placement of multimappers
def assign_multimappers_randomly(
    bam_file,
    output_bam,
    seed=1,
    summary_log_file=None,
    do_sort: bool = True,
    sort_threads: int = 2,
    sort_mem: str = "2G",
):
    """
    Random assignment of one alignment per QNAME with Numba
    Writes a cleaned BAM with NH=1, RC expanded, and ZT tag
    This functuon incorportaes all the changes I have made 
    using numpy and numba instead of pure python
    """
    start_time = time.time()
    np.random.seed(seed)
    print("\nassigning multimappers randomly")

    groups, header = _group_alignments_by_qname(bam_file)

    total_instances = 0
    unique_instances = 0
    multimapper_instances = 0
    total_reads_written = 0
    unique_reads_expanded = 0          
    multimapper_reads_expanded = 0 
    

    with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for rid, alns in tqdm(groups.items(), desc="Randomly resolve multimappers"):
            total_instances += 1

            if len(alns) == 1:
                chosen = alns[0]
                unique_instances += 1
                rc = _get_rc_safe(chosen)
                unique_reads_expanded += int(rc)
            else:
                # FIX: numba helper takes one arg
                chosen_idx = _random_assign_numba(len(alns))
                chosen = alns[chosen_idx]
                multimapper_instances += 1
                rc = _get_rc_safe(chosen)
                multimapper_reads_expanded += int(rc)

            rc = _get_rc_safe(chosen)
            # Sanitize + expand like UWM, with ZT tag
            _expand_write(chosen, out_bam, rc, zt_value="random:1")
            total_reads_written += int(rc)

    # sort & index after expansion (to mirror UWM)
    if do_sort:
        print("Sorting and indexing resolved BAM...")
        sort_and_index_bam(output_bam, threads=sort_threads, sort_mem=sort_mem, inplace=True)

    run_time = time.time() - start_time
    print("\nSummary:")
    print(f"- Collapsed instances analyzed: {total_instances}")
    print(f"- Unique instances kept:         {unique_instances}")
    print(f"- Multimappers assigned:         {multimapper_instances}")
    print(f"- Reads written (expanded):      {total_reads_written}")
    print(f"Extended BAM: {output_bam}")
    print(f"Done in {run_time:.2f} s.\n")

    if summary_log_file:
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "a") as sf:
            sf.write(f"Collapsed instances analyzed: {total_instances}\n")
            sf.write(f"Unique instances kept: {unique_instances}\n")
            sf.write(f"Multimappers assigned: {multimapper_instances}\n")
            sf.write(f"Reads written (expanded): {total_reads_written}\n")
            sf.write(f"Unique reads (expanded): {unique_reads_expanded}\n")
            sf.write(f"Multimapper reads (expanded): {multimapper_reads_expanded}\n")
            sf.write(f"Extended BAM: {output_bam}\n")
            sf.write(f"Random seed: {seed}\n")
            sf.write(f"Run time: {run_time:.2f}\n")

    
def main():
    # Subadd parser for each module, this looks more organized in comparaison to my prevoius CLI
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    
    # Bowtie-build arguments
    p_bt_build = sub.add_parser("bowtie-build", help = "Build bowtie index")
    p_bt_build.add_argument("--genome", help="Genome FASTA file for alignment (Bowtie index)")
    p_bt_build.add_argument("--index", help="Bowtie1 index basename")
    p_bt_build.add_argument("--offrate", type=int, default=4, help="bowtie-build offrate (default= 4)")
    p_bt_build.add_argument("--threads", type=int, default=4, help="Threads for bowtie-build")
    
    # Bowtie-aln arguments
    p_bt_aln = sub.add_parser("bowtie-aln", help = "Build alignemnet")
    p_bt_aln.add_argument("--fastq", help="FASTQ(.gz) file for bowtie-aln")
    p_bt_aln.add_argument("--index", help="Bowtie1 index basename for bowtie-aln")
    p_bt_aln.add_argument("--out", help="Output BAM)")
    p_bt_aln.add_argument("--mismatches", type=int, default=1, help="Bowtie mismatches (default= 1)")
    p_bt_aln.add_argument("--threads", type=int, default=4, help="Threads for bowtie/samtools")
    p_bt_aln.add_argument("--sort-mem", default="2G", help="samtools sort memory per thread (default = 2G)")
    p_bt_aln.add_argument("--rc-map",
                        help="RC map produced collapse (TSV or PKL)")
    p_bt_aln.add_argument("--save-np", choices=["yes", "no"], default="no",
                        help="Save unmapped reads via Bowtie's --un/--un-gz (default=no)")
    p_bt_aln.add_argument("--non-mappers",
                        help="Output fastq(.gz) for unmapped (used when --save-np yes)")
    p_bt_aln.add_argument(
        "--max-multimaps",
        type=int,
        default=None,
        help="Cap reported alignments per read: use N to report up to N alignments (-k N). "
             "If omitted the script uses Bowtie -a (report all multimappers).",
    )

    p_idx = sub.add_parser("build-index-uwm", help="Build UNIQUE weighted index from BAMs")
    p_idx.add_argument("--bams", required=True, help=".bam path OR a text file of .bam paths")
    p_idx.add_argument(
        "--out",
        required=True,
        help=(
            "Output index path. By default a pickle (.pkl) is written; "
            "the script also writes a <prefix>.npidx directory with uncompressed "
            "NumPy (.npy) arrays (memory-map friendly). You may pass either "
            "a .pkl filename or a directory name ending with .npidx to use/load "
            "the NumPy-backed index directly."
        ),
    )
    
    # UWM resolve arguments
    p_res_uwm = sub.add_parser("resolve-uwm", help="Resolve multimappers using the UWM mode and write expanded BAMs by RC (read counts of collapsed reads)")
    p_res_uwm.add_argument(
        "--index",
        required=True,
        help="Path to UNIQUE index: a .pkl file or a .npidx directory with .npy arrays can be used. npidx is preferred for low‑memory loads"
    )
    p_res_uwm.add_argument("--in-bam", required=True, help="Query BAM with uniques+multis")
    p_res_uwm.add_argument("--out-bam", required=True, help="Output BAM (resolved)")
    p_res_uwm.add_argument("--window", type=int, default=250, help="Window width for UWM scoring (default=250)")
    p_res_uwm.add_argument("--strand", action="store_true", help="Consider strand when scoring. If not set, strand-agnostic scoring is used (default=False)")
    p_res_uwm.add_argument("--seed", type=int, default=123, help="Random seed for tie-breaking (default=123)")
    p_res_uwm.add_argument("--threads", type=int, default=4, help="Threads for sorting (default=4)")
    p_res_uwm.add_argument("--sort-mem", default="2G", help="samtools sort memory per thread (default = 2G)")
    p_res_uwm.add_argument("--summary-log", help="Summary log file")
    p_res_uwm.add_argument("--do-sort", type=bool, default=True, help="Disable sorting/indexing after resolution")
    
    
    # Random resolve arguments
    p_res_rand = sub.add_parser("resolve-random", help="Resolve multimappers and write cleaned BAM using random placement")
    p_res_rand.add_argument("--in-bam", required=True, help="Query BAM with uniques+multis")
    p_res_rand.add_argument("--out-bam", required=True, help="Output BAM (resolved random)")
    p_res_rand.add_argument("--seed", type=int, default=1)
    p_res_rand.add_argument("--summary-log", help="Summary log file")
    p_res_rand.add_argument("--threads", type=int, default=4, help="Threads for sorting (default=4)")
    p_res_rand.add_argument("--sort-mem", default="2G", help="samtools sort memory per thread (default = 2G)")
    p_res_rand.add_argument("--do-sort", type=bool, default=True, help="Disable sorting/indexing after resolution")
    
    
    # parse args    
    args = ap.parse_args()
    

    # Commands to run
    if args.cmd == "bowtie-build":
        if not args.genome:
            raise ValueError("--genome is required for bowtie-build")
        index_base = os.path.splitext(os.path.basename(args.genome))[0]
        build_bowtie_index(args.genome, index_base, offrate=args.offrate, threads=args.threads)
    
    # And elif for the rest of the commands
    elif args.cmd == "bowtie-aln":
        if not (args.fastq and args.index):
            raise ValueError("--fastq and --index are required for bowtie-aln")
        # collapsed outputs
        root = args.fastq[:-3] if args.fastq.endswith(".gz") else args.fastq
        root = os.path.splitext(root)[0]

        # define outputs
        mapped_bam        = args.out or (root + ".bam")
        rc_map_path       = args.rc_map  # may be None; then RC=1

        # log files per sample
        summary_log_file  = mapped_bam.replace(".bam", ".map.summary.log")
        printable_cmd     = f"bowtie -q {args.index} {os.path.basename(args.fastq)}"

        # write a small header once
        os.makedirs(os.path.dirname(os.path.abspath(summary_log_file)) or ".", exist_ok=True)
        with open(summary_log_file, "w") as sf:
            sf.write(f"{time.strftime('%c')}\n")
            sf.write("Command:\n\t" + printable_cmd + " | samtools view -b - | samtools sort ...\n")
            sf.write("Completed.\n\n")

        # map (already-collapsed) FASTQ with Bowtie
        tmp_bam = mapped_bam.replace(".bam", ".tmp.bam")
        run_bowtie(
            fastq_file=args.fastq,
            index_base=args.index,
            output_bam=tmp_bam,
            mismatches=args.mismatches,
            save_non_mappers=(args.save_np == "yes"),
            non_mappers_file=args.non_mappers,
            threads=args.threads,
            sort_mem=args.sort_mem,
            index_output_bam=False,
            summary_log_file=summary_log_file,
            max_multimaps=args.max_multimaps,
        )

        # add NH RC using collapse map if provided; else RC=1
        add_NH_RC_tags_from_collapse(tmp_bam, rc_map_path, mapped_bam)
        os.remove(tmp_bam)
        print(f"Added NH/RC tags to: {mapped_bam}")

        # final append
        with open(summary_log_file, "a") as sf:
            sf.write(f"Output BAM: {mapped_bam}\n\n")
    

    elif args.cmd == "build-index-uwm":
        build_weighted_unique_index_from_bams(args.bams, args.out)
        
    elif args.cmd == "resolve-uwm":
        summary_log = args.summary_log or args.out_bam.replace(".bam", ".uwm.summary.log")
        resolve_multimappers_uwm(
            index_path=args.index,
            query_bam_path=args.in_bam,
            output_bam_path=args.out_bam,
            window_size=args.window,
            consider_strand=args.strand,
            seed=args.seed,
            progress=True,
            summary_log_file=summary_log,
            do_sort=True,
            sort_threads=args.threads if hasattr(args, "threads") else 2,
            sort_mem=args.sort_mem,
        )

    elif args.cmd == "resolve-random":
        summary_log_file = args.summary_log or (os.path.splitext(os.path.basename(args.in_bam))[0] + ".map.summary.log")
        assign_multimappers_randomly(
        bam_file=args.in_bam,
        output_bam=args.out_bam,
        seed=args.seed,
        summary_log_file=summary_log_file,
        do_sort=True,
        sort_threads=args.threads,
        sort_mem=args.sort_mem,
        )



if __name__ == "__main__":
    main()

