#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def CYAN  = "\u001B[36m"
def RESET = "\u001B[0m"

def BANNER = $/
       ___  _  _____                         ___
  ___ / _ \/ |/ / _ | ___ ___ ___ _    ___  / _/
 (_-</ , _/    / __ |(_-</ -_) _ `/___/ _ \/ _/ 
/___/_/|_/_/|_/_/ |_/___/\__/\_, /   /_//_/_/   
                              /_/               
/$

log.info("\n${CYAN}${BANNER}${RESET}")


// Include the scripts for feature + fisrt nt and test


params.reads       = params.reads       ?: "${launchDir}/test/*.fastq.gz"
params.threads     = params.threads     ?: 2
params.preproc     = params.preproc     ?: 'legacy'   // 'fastp' or 'legacy'
params.adapter     = params.adapter     ?: "AGATCGGAAGAG"
params.fastp_use_adapter = (params.fastp_use_adapter != null ? params.fastp_use_adapter : false)
params.minlen      = params.minlen      ?: 18
params.maxlen      = params.maxlen      ?: 27
params.map_gate    = params.map_gate    ?: 'none' // 'trure' to save resources, none to skip the gate and increse speed and resorce usge
params.genome      = params.genome      ?: "${launchDir}/genome/caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa"
params.annotation  = params.annotation  ?: "${launchDir}/annotation/caenorhabditis_elegans.PRJNA13758.WBP19.overlapping_annotation.gff3"
params.max_multimaps = params.max_multimaps ?: null
if (params.max_multimaps != null) {
    try {
        int value = Integer.parseInt(params.max_multimaps.toString())
        if (value < 1) {
            params.max_multimaps = null // Reset to null if not a positive integer
        } else {
            params.max_multimaps = value // Keep the valid integer
        }
    } catch (NumberFormatException e) {
        params.max_multimaps = null // Reset to null if not a valid integer
    }
}
params.make_bedgraph = (params.make_bedgraph != null ? params.make_bedgraph : true) ?: null
if (params.make_bedgraph instanceof String) {
    params.make_bedgraph = params.make_bedgraph.toLowerCase() in ['true','1','yes','y']
}
params.expand_unmapped = params.expand_unmapped ?: false //boolean
params.n_mismatch    = params.n_mismatch    ?: 1 // integer value (bowtie allows up to 3 mismatching bases)
params.offrate_sm  = params.offrate_sm  ?: 4
params.thr_sm      = params.thr_sm      ?: 12
params.smem_sm     = params.smem_sm     ?: "12G"
params.wins_sm     = params.wins_sm     ?: 250
params.consider_strand = params.consider_strand ?: false
params.assign_mode    = params.assign_mode    ?: 'uwm'   // 'uwm' or 'random'
params.minoverlap     = params.minoverlap     ?: 0.7
params.contrast       = params.contrast       ?: "${launchDir}/contrast.txt"
params.lfc            = params.lfc            ?: 1
params.fdr            = params.fdr            ?: 0.05
params.threshold_inc   = params.threshold_inc   ?: false   // read FC/FDR from contrast file?
params.hk_norm	      = params.hk_norm        ?: false // boolean
params.norm_feature   = params.norm_feature   ?: "rRNA_S"
params.stringent_tmm  = params.stringent_tmm  ?: false   // boolean
params.use_rds        = params.use_rds        ?: true   // boolean
params.first_nt       = params.first_nt       ?: ""
params.apply_first_nt_downstream = params.apply_first_nt_downstream ?: false
params.srcDir         = "${workflow.projectDir}/bin"
reads_ch              = Channel.fromPath(params.reads)
genome_ch             = Channel.value( file(params.genome) )
annotation_ch         = Channel.fromPath(params.annotation)
contrast_ch           = Channel.fromPath(params.contrast)
siRmap_script_ch      = Channel.value( file("${params.srcDir}/siRmap.py") )
collapse_script_ch    = Channel.value( file("${params.srcDir}/collapse") )
summary_script_uwm_ch   = Channel.fromPath("${params.srcDir}/02.summary_uwm.R")
summary_script_rand_ch  = Channel.fromPath("${params.srcDir}/02.summary_rand.R")
bam2Rds_script_ch	= Channel.fromPath("${params.srcDir}/01.bam2Rds.R")
fnmtx_script_ch         = Channel.fromPath("${params.srcDir}/02.get_fn_mtx.R")
featureCounts_script_ch = Channel.fromPath("${params.srcDir}/02.featureCounts.R")
dea_script_ch           = Channel.fromPath("${params.srcDir}/04.edgeR.R")


def results_dir = "Results_${new Date().format('yyyyMMdd_HHmmss')}"
new File(results_dir).mkdirs()


process fastp {
  label 'preproc'

  input:
  path read


  output:
  path "*.html", emit: qc_html
  path "*.json", emit: qc_json
  path "*.ps.fastq.gz", emit: fastq

  tag "${read.simpleName}"

  publishDir "${results_dir}/01.fastp_qc", mode: 'copy'

  script:
  """
  base='${read.simpleName}'
  base=\${base%.fastq.gz}; base=\${base%.fq.gz}; base=\${base%.fastq}; base=\${base%.fq}

  ADAPTER_OPT=""
  if [ "${params.fastp_use_adapter}" = "true" ] && [ -n "${params.adapter}" ]; then
    ADAPTER_OPT="--adapter_sequence ${params.adapter}"
  fi

  fastp \
  -i ${read} \
  -o "\${base}.ps.fastq.gz" \
  --length_required ${params.minlen} \
  --length_limit ${params.maxlen} \
  -e 25 \
  -q 20 \
  -u 10 \
  -n 1 \
  -w ${params.threads} \
  -h "\${base}.html" \
  -j "\${base}.json" \
  \$ADAPTER_OPT

  """
}

process multiqc_fastp {
  label 'preproc'
  tag   "MultiQC (fastp)"

  input:
    path fastp_reports  // html/json from fastp

  output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

  publishDir "${results_dir}/01.fastp_qc", mode: 'move'

  script:
  """
  multiqc . --outdir . --module fastp
  """
}


process fastqc { label 'preproc'
    input:
    path read  
   
    output:
    path "*.html", emit: qc_html    
    path "*.zip" , emit: qc_zip
    
    publishDir "${results_dir}/01.raw_qc", mode: 'copy'

    script:
    """
    fastqc --threads ${params.threads} ${read}
    """
}

process multiqc { label 'preproc'
    tag "MultiQC summary"

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

    publishDir "${results_dir}/01.raw_qc", mode: 'move'

    script:
    """
    multiqc . --outdir .
    """
}

process cutadapt { label 'preproc'
    input:
    path read

    output:
    path "*trimmed.fastq.gz", emit: fastq

    publishDir "${results_dir}/02.cut_adapt", mode: 'copy'

    tag "${read.simpleName}"

    script:
    """
    cutadapt -j ${params.threads} \
             -a ${params.adapter} \
             -m ${params.minlen} \
             --max-n 0.05 \
             --discard-untrimmed \
             -o ${read.simpleName}.trimmed.fastq.gz ${read}
    """
}

process fastqc_trimm { label 'preproc'
    input:
    path read  
   
    output:
    path "*.html", emit: qc_html    
    path "*.zip" , emit: qc_zip
    
    publishDir "${results_dir}/03.trimmed_qc", mode: 'copy'

    script:
    """
    fastqc --threads ${params.threads} ${read}
    """
}

process multiqc_tr { label 'preproc'
    tag "MultiQC summary"

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

    publishDir "${results_dir}/03.trimmed_qc", mode: 'move'

    script:
    """
    multiqc . --outdir .
    """
}

process collapse {
  label 'preproc'
  tag { fastq.simpleName }

  input:
  path fastq
  path collapse_script

  output:
  path "${fastq.simpleName}.collapsed.fastq.gz", emit: collapsed_fq
  path "${fastq.simpleName}.map.tsv.gz",            emit: map_tsv

  publishDir "${results_dir}/05.map", mode: 'copy'

  script:
  """
  ${collapse_script} \
    -i ${fastq} \
    -f ${fastq.simpleName}.collapsed.fastq.gz \
    -M ${fastq.simpleName}.map.tsv
  pigz ${fastq.simpleName}.map.tsv -p ${params.threads}
  """
}

process pullseq {
  label 'preproc'

  input:
    tuple path(read), path(rc_map)

  output:
    tuple path("*ps.fastq.gz"), path(rc_map), emit: fastq_map

  tag "${read.simpleName}"
  publishDir "${results_dir}/04.pullseq", mode: 'copy'

  script:
  """
  base='${read.simpleName}'
  base=\${base%.fastq}

  pullseq -i ${read} -m ${params.minlen} -a ${params.maxlen} > "\${base}.ps.fastq"
  pigz "\${base}.ps.fastq"
  """
}

process build_index {
  tag "Bowtie index"
  input:
    path siRmap_script
    path genome
  output:
    path "*.ebwt", emit: ebwt
  script:
  """
  IDX_BASE=\$(basename "${genome}")
  IDX_BASE=\${IDX_BASE%.fa}
  IDX_BASE=\${IDX_BASE%.fasta}
  IDX_BASE=\${IDX_BASE%.fna}
  IDX_BASE=\${IDX_BASE%.fa.gz}
  IDX_BASE=\${IDX_BASE%.fasta.gz}
  IDX_BASE=\${IDX_BASE%.fna.gz}

  python ${siRmap_script} \
          bowtie-build \
          --genome ${genome} \
          --index \${IDX_BASE} \
          --offrate ${params.offrate_sm} \
          --threads ${params.thr_sm}
  """
}

process map_collapse {
  tag { fastq.baseName }
  cpus params.thr_sm
  memory params.smem_sm

  input:
    tuple path(fastq), path(rc_map_tsv)
    path genome
    path siRmap_script
    path bowtie_idx         

  output:
    path "*.collapsed.bam",               emit: collapsed_bam
    path "*.unmapped.collapsed.fastq.gz", emit: unmapped_fq
    path "*.log",                         emit: logs

  publishDir "${results_dir}/05.map", mode: 'copy'

  script:
  """
  IDX_BASE=\$(basename "${genome}")
  IDX_BASE=\${IDX_BASE%.fa}
  IDX_BASE=\${IDX_BASE%.fasta}
  IDX_BASE=\${IDX_BASE%.fna}
  IDX_BASE=\${IDX_BASE%.fa.gz}
  IDX_BASE=\${IDX_BASE%.fasta.gz}
  IDX_BASE=\${IDX_BASE%.fna.gz}

  BASENAME=\$(basename "${fastq}")
  BASENAME=\${BASENAME%.gz}
  BASENAME=\${BASENAME%.fastq}
  BASENAME=\${BASENAME%.fq}
  BASENAME=\${BASENAME%.collapsed}  
  BASENAME=\${BASENAME%.ps}  

  python ${siRmap_script} bowtie-aln \
    --fastq ${fastq} \
    --index "\$IDX_BASE" \
    --mismatches ${params.n_mismatch} \
    --out "\${BASENAME}.collapsed.bam" \
    --save-np yes  \
    --non-mappers "\${BASENAME}.unmapped.collapsed.fastq.gz" \
    --threads ${params.thr_sm} \
    --sort-mem ${params.smem_sm} \
    --rc-map ${rc_map_tsv} \
    ${params.max_multimaps ? "--max-multimaps ${params.max_multimaps}" : ""}
  """
}

process expand_unmapped {

  input:
  tuple path(unmapped_fq), path(rc_map_tsv)
  path siRmap_script

  output:
  path "*.unmapped.expanded.fastq.gz"
  publishDir "${results_dir}/05.map", mode: 'copy'

  script:
  """
  BASENAME=\$(basename "${unmapped_fq}")
  BASENAME=\${BASENAME%.gz}
  BASENAME=\${BASENAME%.fastq}
  BASENAME=\${BASENAME%.fq}
  BASENAME=\${BASENAME%.collapsed}  
  BASENAME=\${BASENAME%.unmapped}  
  
  python ${siRmap_script} expand-fq \
    --fastq ${unmapped_fq} \
    --rc-map ${rc_map_tsv} \
    --out "\${BASENAME}.unmapped.expanded.fastq.gz"
  """

}

process resolve_random {
  tag { bam.baseName }
  cpus params.thr_sm
  memory params.smem_sm

  input:
    path bam
    path siRmap_script
  when:
    params.assign_mode == 'random'

  output:
    path "*.expanded.bam", emit: expanded_bam
    path "*.log",          emit: logs

  publishDir "${results_dir}/06.resolved_random", mode: 'copy'

  script:
  """
  SAMPLE=\$(basename "${bam}" .bam)
  python ${siRmap_script} resolve-random \
    --in-bam ${bam} \
    --out-bam "\${SAMPLE}.expanded.bam" \
    --seed 123 \
    --threads ${params.threads} \
    --sort-mem ${params.smem_sm} \
    --do-sort true
  """
}

process build_unique_index {
  tag "unique-index"

  input:
    path siRmap_script
    path bam_paths

  output:
    path "unique_index.pkl", emit: uniq_idx


  publishDir "${results_dir}/06.uwm_index", mode: 'copy'
  script:
  """
  printf "%s\n" ${bam_paths.join(' ')} | tr ' ' '\\n' > bam_list.txt
  python ${siRmap_script} build-index-uwm --bams bam_list.txt --out unique_index.pkl
  """
}

process resolve_uwm { 
  tag { bam.baseName }
  cpus params.thr_sm
  memory params.smem_sm

  input:
    path bam
    path uniq_idx
    path siRmap_script
  when:
    params.assign_mode == 'uwm'

  output:
    path "*.expanded.bam", emit: expanded_bam
    path "*.log",          emit: logs
  publishDir "${results_dir}/06.resolved_uwm", mode: 'copy'

  script:
  """
  SAMPLE=\$(basename "${bam}" .bam)
  python ${siRmap_script} resolve-uwm \
    --in-bam ${bam} \
    --index ${uniq_idx} \
    --window ${params.wins_sm} \
    --out-bam "\${SAMPLE}.expanded.bam" \
    --threads ${params.thr_sm} \
    --sort-mem ${params.smem_sm} \
    --seed 123 ${ params.consider_strand ? '--strand' : '' } 
  """
}

process summarize_rand {
    tag "summary plots rand"
    input:
    path log_files
    path summary_script_rand
    output:
    path "*.png", emit: png
    publishDir "${results_dir}/06.summary", mode: 'copy'
    script:
    """
    Rscript ${summary_script_rand}
    """
}

process summarize_uwm {
    tag "summary plots uwm"
    input:
    path log_files
    path summary_script_uwm
    output:
    path "*.png", emit: png
    publishDir "${results_dir}/06.summary", mode: 'copy'
    script:
    """
    Rscript ${summary_script_uwm}
    """
}


process bam2Rds {
    tag "bam2Rds"
    input:
    path bam2Rds_script
    path bam_files

    output:
    path "*.Rds", emit: matrices

    publishDir "${results_dir}/07.rds", mode: 'copy'

    script:
    """
    Rscript ${bam2Rds_script} ${bam_files.join(' ')}
    """
}

process firstnt_counts {
  tag { bam.simpleName }

  input:
    path bam

  output:
    path "*.firstnt.tsv", emit: counts

  publishDir "${results_dir}/09.fn_counts", mode: 'copy'

  script:
  """
  base=\$(basename "${bam}" .bam)
  out="\${base}.firstnt.tsv"
  samtools view -h -F 0x904 -@ ${params.threads} "${bam}" \
  | awk -v MIN=${params.minlen} -v MAX=${params.maxlen} '
      BEGIN { }
      /^@/ { next }
      {
        seq = toupper(\$10)
        l = length(seq)
        if (l < MIN || l > MAX) next

        flag = \$2 + 0
        b = substr(seq,1,1)

        if ( int(flag/16) % 2 == 1 ) {
          last = substr(seq,l,1)
          b = (last=="A"?"T": last=="C"?"G": last=="G"?"C": last=="T"?"A": last)
        }

        if (b=="A"||b=="C"||b=="G"||b=="T") {
          cnt[l ":" b]++
          tot[l]++
        }
      }
      END {
        print "sample\\tlength\\tA\\tC\\tG\\tT\\ttotal"
        for (l=MIN; l<=MAX; l++) {
          a = cnt[l ":A"] + 0
          c = cnt[l ":C"] + 0
          g = cnt[l ":G"] + 0
          t = cnt[l ":T"] + 0
          print "__SAMPLE__\\t" l "\\t" a "\\t" c "\\t" g "\\t" t "\\t" (tot[l]+0)
        }
      }
    ' > "\$out.tmp"

  sed "s/^__SAMPLE__/\${base}/" "\$out.tmp" > "\$out"
  rm -f "\$out.tmp"
  """
}

process fn_mtx {
    tag "Get first-nucleotide matrices"

    input:
    path fnmtx_script
    path rds_files

    output:
    path "first_nt_rlength.Rds", emit: matrices
    publishDir "${results_dir}/08.fn_mtx", mode: 'copy'

    script:
    """
    Rscript ${fnmtx_script}
    """
}

process plot_firstnt {
  tag "Plot first nucleotide distributions"

  input:
    path tables   // <- a list; Nextflow will stage all of them into the work dir

  output:
    path "length_dit_fn_percentage.png", emit: png_percentage

  publishDir "${results_dir}/09.fn_plots", mode: 'copy'

  script:
  """
  Rscript ${params.srcDir}/03.plot_fn.R ${params.minlen} ${params.maxlen} "*.firstnt.tsv"
  """
}


process filter_firstnt {
  when:
    params.first_nt != null && params.first_nt.toString().trim().length() > 0

  tag { "${bam_file.simpleName}" }

  input:
    path bam_file

  output:
    path "*.nt*.trim.mapped.bam", emit: filtered_bams

  publishDir "${results_dir}/07.bams_firstnt", mode: 'copy'

  script:
  """
  letters_raw='${params.first_nt}'

  letters=\$(echo "\$letters_raw" | tr 'u' 't' | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')

  if [ -z "\$letters" ] || ! echo "\$letters" | grep -Eq '^[ACGT]+\$'; then
    echo "ERROR: first_nt must be A|C|G|T (e.g. 'G' or 'AT'). Got: '\$letters_raw'" >&2
    exit 1
  fi

  base=\$(basename "${bam_file}" .bam)
  outFile="\${base}.nt\${letters}.trim.mapped.bam"

  LC_ALL=C samtools view -h -F 0x904 -@ ${params.threads} "${bam_file}" \
  | awk -v set="\$letters" '
      BEGIN {
        for (i=1; i<=length(set); i++) allow[toupper(substr(set,i,1))]=1
      }
      function comp(b){ return (b=="A"?"T": b=="C"?"G": b=="G"?"C": b=="T"?"A": b) }
      /^@/ { print; next }
      {
        seq = toupper(\$10)
        b = substr(seq,1,1)
        if ( and(\$2,16) ) b = comp(substr(seq,length(seq),1))
        if (allow[b]) print
      }
    ' \
  | samtools view -b -@ ${params.threads} -o "\$outFile" -
  """
}

process featureCounts {
    tag "featureCounts"
    input:
    path bam_files
    path annotation
    path featureCounts_script
    
    output:
    path "*.txt", emit: table_of_counts
    path "*.featureCounts", emit: featc_read_tables   
    path "*.RDS", emit: featurecount_rds
    publishDir "${results_dir}/11.featureCounts", mode: 'copy'
    
    script:
    """
    Rscript ${featureCounts_script} ${annotation} ${params.minoverlap} ${params.thr_sm} ${bam_files.join(" ")}
    """
}

process bam2bedgraph {
    tag "${bam_file.simpleName}"

    input:
    path bam_file

    output:
    path "*.bedGraph.gz", emit: bedgraphs

    publishDir "${results_dir}/12.bedGraphs", mode: 'copy'

    script:
    """
    for file in *.bam; do
        echo "Processing \$file"
        outFile=\$(echo \$file | sed 's/\\..*//')
        
        total=\$(samtools view -c -F 0x904 \$file)
        normScale=\$(bc <<< "scale=4;1000000/\$total")

        bedtools genomecov -ibam \$file -bg -scale \$normScale -strand '+' > \$outFile.us
        bedtools genomecov -ibam \$file -bg -scale \$normScale -strand '-' | sed -r 's/(\\S+)\$/-\\1/' >> \$outFile.us

        bedtools sort -i \$outFile.us | pigz -c > \$outFile.bedGraph.gz
        rm \$outFile.us
    done
    """
}

process edgeR_dea {
    tag "edgeR DEA"

    input:
    path dea_script
    path counts_txt                               
    path contrast_file

    output:
    path "tmm",     emit: tmm_dir
    path "norm",    emit: norm_dir
    path "dge",     emit: dge_dir
    path "figures", emit: figs_dir
    publishDir "${results_dir}/13.edgeR", mode: 'copy'

    script:
    def STRICT = params.stringent_tmm ? 'TRUE' : 'FALSE'
    def HKNORM = params.hk_norm ? 'TRUE' : 'FALSE'
    def TRESH  = params.threshold_inc ? 'TRUE' : 'FALSE' 
    """
    Rscript ${dea_script} \
      --counts_file ${counts_txt} \
      --contrast_file ${contrast_file} \
      --threshold_inc ${TRESH} \
      --logfold_change ${params.lfc} \
      --fdr ${params.fdr} \
      --hk_norm ${HKNORM} \
      --norm_feature '${params.norm_feature}' \
      --stringent_tmm ${STRICT}
    """
}

workflow {  

  def map_inputs_ch
  if ( params.preproc == 'fastp' ) {
   
    fp_out = fastp(reads_ch)
    multiqc_fastp( fp_out.qc_json.collect() )

    def collapsed = collapse(fp_out.fastq, collapse_script_ch)

    def fq_join = collapsed.collapsed_fq.map { fq ->
      def id = fq.simpleName.replace('.collapsed.fastq','')
      tuple(id, fq)
    }

    def map_join = collapsed.map_tsv.map { m ->
      def id = m.simpleName.replace('.map.tsv','')
      tuple(id, m)
    }

    map_inputs_ch = fq_join.join(map_join).map { sid, fq, m -> tuple(fq, m) }
    // dup map_join
    rc_map_ch = map_join

  } else {
   
    fastqc_raw = fastqc(reads_ch)
    multiqc( fastqc_raw.qc_zip.collect() )

    trimmed   = cutadapt(reads_ch)
    fastqc_tr = fastqc_trimm(trimmed)
    multiqc_tr( fastqc_tr.qc_zip.collect() )

    def collapsed  = collapse(trimmed, collapse_script_ch)

    def fq_join = collapsed.collapsed_fq.map { fq ->
      def id = fq.simpleName.replace('.collapsed.fastq','')
      tuple(id, fq)
    }

    def map_join = collapsed.map_tsv.map { m ->
      def id = m.simpleName.replace('.map.tsv','')
      tuple(id, m)
    }

    def collapsed_pairs = fq_join.join(map_join).map { sid, fq, m -> tuple(fq, m) }

    def pulled = pullseq(collapsed_pairs)
    map_inputs_ch = pulled.fastq_map
    // dup map_join 
    rc_map_ch = map_join
  }

  
  def fastq_map_for_map
  if( params.map_gate == 'all' ) {
    def all_done = map_inputs_ch.collect().map { true }
    fastq_map_for_map = map_inputs_ch
                      .combine(all_done)
                      .map { pair, _ -> pair }
  } else {
    fastq_map_for_map = map_inputs_ch
  }

  idx = build_index(siRmap_script_ch, genome_ch)

  mapped = map_collapse(
            fastq_map_for_map,
            genome_ch,
            siRmap_script_ch,
            idx.ebwt.collect()
          )
  // new if to expand unpaed reads by ID in case user need this
  if( params.expand_unmapped ) {
    def unmapped_join = mapped.unmapped_fq.map { f ->
      def id = f.simpleName.replace('.unmapped.collapsed.fastq','')
      tuple(id, f)
    }

    def unmapped_pairs = unmapped_join.join(rc_map_ch).map { sid, unmapped_fq, rc_map_tsv ->
      tuple(unmapped_fq, rc_map_tsv)
    }

    expand_unmapped(unmapped_pairs, siRmap_script_ch)
  }

  resolved_bams_ch = Channel.empty()
  uniq_idx_ch      = Channel.empty()
  if ( params.assign_mode == 'uwm' ) {
    uniq_idx_ch = build_unique_index(
                    siRmap_script_ch,
                    mapped.collapsed_bam.collect()
                  ).uniq_idx

    def uwm_out = resolve_uwm(
                    mapped.collapsed_bam,
                    uniq_idx_ch,
                    siRmap_script_ch
                  )

    resolved_bams_ch = uwm_out.expanded_bam

    summarize_uwm(
      uwm_out.logs.collect(),
      summary_script_uwm_ch
    )

  } else if ( params.assign_mode == 'random' ) {
    def rand_out = resolve_random(
                     mapped.collapsed_bam,
                     siRmap_script_ch
                   )

    resolved_bams_ch = rand_out.expanded_bam

    summarize_rand(
      rand_out.logs.collect(),
      summary_script_rand_ch
    )

  } else {
    resolved_bams_ch = mapped.collapsed_bam
  }
  
  if( params.use_rds ) {
    rds_out   = bam2Rds(bam2Rds_script_ch, resolved_bams_ch.collect())
  }

  def counts_out = firstnt_counts(resolved_bams_ch)
    def all_counts = counts_out.counts.collect()
    plot_firstnt(all_counts)
  

  if( !params.apply_first_nt_downstream ) {
    log.info "Stopping after first-nucleotide plots (apply_first_nt_downstream=false)."
    return
  }

  def want_filter = params.first_nt != null && params.first_nt.toString().trim().length() > 0
  def bams_for_quant = resolved_bams_ch

  if( want_filter ) {
    def filtered = filter_firstnt(resolved_bams_ch)
    bams_for_quant = filtered.filtered_bams
  }
  fc = featureCounts(bams_for_quant.collect(), annotation_ch, featureCounts_script_ch)
  counts_only = fc.table_of_counts.filter { it.name == 'featureCounts_counts.txt' }
  if (params.make_bedgraph) {
  bam2bedgraph(bams_for_quant)
  log.info "bedGraph generation: enabled"
  } else {
  log.info "bedGraph generation: disabled (make_bedgraph=NULL)"
  }
  edgeR_dea(dea_script_ch, counts_only, contrast_ch)
}

workflow.onComplete {
    println "Run name : ${workflow.runName}"
    println "Duration : ${workflow.duration}"
    println "Succeeded: ${workflow.success}"
}
