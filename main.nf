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


params.reads       = params.reads       ?: "${launchDir}/test/*.fastq.gz"
params.threads     = params.threads     ?: 2
params.adapter     = params.adapter     ?: "AGATCGGAAGAG"
params.minlen      = params.minlen      ?: 18
params.maxlen      = params.maxlen      ?: 27
params.genome      = params.genome      ?: "${launchDir}/genome/caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa"
params.annotation  = params.annotation  ?: "${launchDir}/annotation/caenorhabditis_elegans.PRJNA13758.WBP19.overlapping_annotation.gff3"
params.offrate_sm  = params.offrate_sm  ?: 4
params.thr_sm      = params.thr_sm      ?: 12
params.smem_sm     = params.smem_sm     ?: "12G"
params.wins_sm     = params.wins_sm     ?: 200
params.consider_strand = params.consider_strand ?: false
params.assign_mode = params.assign_mode ?: 'uwm'   // 'uwm' or 'random'
params.minoverlap       = params.minoverlap     ?: 0.7
params.contrast       = params.contrast       ?: "${launchDir}/contrast.txt"
params.lfc            = params.lfc            ?: 1
params.fdr            = params.fdr            ?: 0.05
params.treshold_inc   = params.treshold_inc   ?: false   // read FC/FDR from contrast file?
params.hk_norm	      = params.hk_norm        ?: false // boolean
params.norm_feature   = params.norm_feature   ?: "rRNA_S"
params.stringent_tmm  = params.stringent_tmm  ?: false   // boolean

params.first_nt = params.first_nt ?: ""
params.apply_first_nt_downstream = params.apply_first_nt_downstream ?: false
params.srcDir = "${workflow.projectDir}/bin"
reads_ch = Channel.fromPath(params.reads)
genome_ch        = Channel.value( file(params.genome) )
annotation_ch = Channel.fromPath(params.annotation)
contrast_ch = Channel.fromPath(params.contrast)
siRmap_script_ch = Channel.value( file("${params.srcDir}/siRmap.py") )
summary_script_uwm_ch   = Channel.fromPath("${params.srcDir}/02.summary_uwm.R")
summary_script_rand_ch  = Channel.fromPath("${params.srcDir}/02.summary_rand.R")
bam2Rds_script_ch	= Channel.fromPath("${params.srcDir}/01.bam2Rds.R")
fnmtx_script_ch         = Channel.fromPath("${params.srcDir}/02.get_fn_mtx.R")
featureCounts_script_ch = Channel.fromPath("${params.srcDir}/02.featureCounts.R")
dea_script_ch = Channel.fromPath("${params.srcDir}/04.edgeR.R")
process fastqc {
    input:
    path read  
   
    output:
    path "*.html", emit: qc_html    
    path "*.zip" , emit: qc_zip
    
    publishDir "01.raw_qc", mode: 'copy'

    script:
    """
    fastqc --threads ${params.threads} ${read}
    """
}

process multiqc {
    tag "MultiQC summary"

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

    publishDir "01.raw_qc", mode: 'move'

    script:
    """
    multiqc . --outdir .
    """
}

process cutadapt {
    input:
    path read

    output:
    path "*trimmed.fastq.gz", emit: fastq

    publishDir "02.cut_adapt", mode: 'copy'

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

process fastqc_trimm {
    input:
    path read  
   
    output:
    path "*.html", emit: qc_html    
    path "*.zip" , emit: qc_zip
    
    publishDir "03.trimmed_qc", mode: 'copy'

    script:
    """
    fastqc --threads ${params.threads} ${read}
    """
}

process multiqc_tr {
    tag "MultiQC summary"

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

    publishDir "03.trimmed_qc", mode: 'move'

    script:
    """
    multiqc . --outdir .
    """
}

process pullseq {
    input:
    path read

    output:
    path "*ps.fastq.gz", emit: fastq
  
    tag "${read.simpleName}"

    publishDir "04.pullseq", mode: 'copy'
    script:
    """
    pullseq -i ${read} -m ${params.minlen} -a ${params.maxlen} > ${read.simpleName}.ps.fastq
    pigz ${read.simpleName}.ps.fastq
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
  python ${siRmap_script} \
          --mode build \
          --genome ${genome} \
          --offrate ${params.offrate_sm} \
          --threads ${params.thr_sm}
  """
}


process map_collapse {
  tag { fastq.baseName }
  cpus params.thr_sm
  memory params.smem_sm

  input:
    path fastq            
    path genome           
    path siRmap_script    
    path bowtie_idx      

  output:
    path "*.collapsed.bam",                 emit: collapsed_bam
    path "*.collapsed.fastq.gz",            emit: collapsed_fq
    path "*.rc.pkl",                        emit: rc_pkl
    path "*.unmapped.collapsed.fastq.gz",   emit: unmapped_fq
    path "*.log",                           emit: logs

  publishDir "05.map", mode: 'copy'

  script:
  """
  IDX_BASE=\$(basename "${genome}")
  IDX_BASE=\${IDX_BASE%.fa}
  IDX_BASE=\${IDX_BASE%.fasta}
  IDX_BASE=\${IDX_BASE%.fna}
  IDX_BASE=\${IDX_BASE%.fa.gz}
  IDX_BASE=\${IDX_BASE%.fasta.gz}
  IDX_BASE=\${IDX_BASE%.fna.gz}

  python ${siRmap_script} --mode map \
    --fastq ${fastq} \
    --mismatches 1 \
    --index "\$IDX_BASE" \
    --out ${fastq.baseName}.collapsed.bam \
    --save-np yes \
    --non-mappers ${fastq.baseName.replace('.fastq','')}.unmapped.collapsed.fastq.gz \
    --threads ${params.thr_sm} \
    --sort-mem ${params.smem_sm}
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

  publishDir "06.resolved_random", mode: 'copy'

  script:
  """
  SAMPLE=\$(basename "${bam}" .bam)
  python ${siRmap_script} --mode random \
    --bam ${bam} \
    --resolved-out "\${SAMPLE}.expanded.bam" \
    --seed 1
  """
}

process build_unique_index {
  tag "unique-index"

  input:
    path siRmap_script
    path bam_paths

  output:
    path "unique_index.pkl", emit: uniq_idx


  publishDir "06.uwm_index", mode: 'copy'
  script:
  """
  printf "%s\n" ${bam_paths.join(' ')} | tr ' ' '\\n' > bam_list.txt

  python ${siRmap_script} --mode index --bam-list bam_list.txt --unique-index unique_index.pkl
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
  publishDir "06.resolved_uwm", mode: 'copy'

  script:
  """
  SAMPLE=\$(basename "${bam}" .bam)
  python ${siRmap_script} --mode uwm \
    --bam ${bam} \
    --unique-index ${uniq_idx} \
    --window-size ${params.wins_sm} \
    --resolved-out "\${SAMPLE}.expanded.bam" \
    --seed 123 ${ params.consider_strand ? '--consider-strand' : '' }
  """
}

process summarize_rand {
    tag "summary plots rand"
    input:
    path log_files
    path summary_script_rand
    output:
    path "*.png", emit: png
    publishDir "06.summary", mode: 'copy'
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
    publishDir "06.summary", mode: 'copy'
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

    publishDir "07.rds", mode: 'copy'

    script:
    """
    Rscript ${bam2Rds_script} ${bam_files.join(' ')}
    """
}

process fn_mtx {
    tag "Get first-nucleotide matrices"

    input:
    path fnmtx_script
    path rds_files

    output:
    path "first_nt_rlength.Rds", emit: matrices
    publishDir "08.fn_mtx", mode: 'copy'

    script:
    """
    Rscript ${fnmtx_script}
    """
}

process plot_firstnt {
    tag "Plot first nucleotide distributions"
    input:
    path txt_files

    output:
    path "length_dit_fn_percentage.png", emit: png_percentage

    publishDir "09.fn_plots", mode: 'copy'

    script:
    """
    Rscript ${params.srcDir}/03.plot_fn.R ${params.minlen} ${params.maxlen}
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
  publishDir "10.bams_firstnt", mode: 'copy'

  script:
  """
  letters_raw='${params.first_nt}'
  letters=\$(echo "\$letters_raw" | tr 'u' 't' | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')

  if [ -z "\$letters" ] || ! echo "\$letters" | grep -Eq '^[ACGT]+\$'; then
    echo "ERROR: --first_nt must be A/C/G/T (e.g. 'G' or 'AT'). Got: '\$letters_raw'" >&2
    exit 1
  fi

  base=\$(basename "${bam_file}" .bam)
  outFile="\${base}.nt\${letters}.trim.mapped.bam"

  samtools view -h -F 4 -@ ${params.threads} "${bam_file}" \
  | awk -v set="\$letters" '
      BEGIN { for (i=1;i<=length(set);i++) allow[toupper(substr(set,i,1))]=1 }
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
    publishDir "11.featureCounts", mode: 'copy'
    
    script:
    """
    Rscript ${featureCounts_script} ${annotation} ${params.minoverlap} ${params.threads} ${bam_files.join(" ")}
    """
}

process bam2bedgraph {
    tag "${bam_file.simpleName}"

    input:
    path bam_file

    output:
    path "*.bedGraph.gz", emit: bedgraphs

    publishDir "12.bedGraphs", mode: 'copy'

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
    publishDir "13.edgeR", mode: 'copy'

    script:
    def STRICT = params.stringent_tmm ? 'TRUE' : 'FALSE'
    def HKNORM = params.hk_norm ? 'TRUE' : 'FALSE'
    def TRESH  = params.treshold_inc ? 'TRUE' : 'FALSE' 
    """
    Rscript ${dea_script} \
      --counts_file ${counts_txt} \
      --contrast_file ${contrast_file} \
      --treshold_inc ${TRESH} \
      --logfold_change ${params.lfc} \
      --fdr ${params.fdr} \
      --hk_norm ${HKNORM} \
      --norm_feature '${params.norm_feature}' \
      --stringent_tmm ${STRICT}
    """
}

workflow {
  fastqc_raw = fastqc(reads_ch)
  multiqc( fastqc_raw.qc_zip.collect() )
  trimmed   = cutadapt(reads_ch)
  fastqc_tr = fastqc_trimm(trimmed)
  multiqc_tr( fastqc_tr.qc_zip.collect() )
  pulled = pullseq(trimmed)
  idx = build_index(siRmap_script_ch, genome_ch)
  mapped = map_collapse(
            pulled.fastq,
            genome_ch,
            siRmap_script_ch,
            idx.ebwt.collect()
          )

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

  rds_out   = bam2Rds(bam2Rds_script_ch, resolved_bams_ch.collect())
  fnmtx_out = fn_mtx(fnmtx_script_ch, rds_out.matrices)
  plot_firstnt(fnmtx_out.matrices)
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
  bam2bedgraph(bams_for_quant)
  edgeR_dea(dea_script_ch, counts_only, contrast_ch)
}

workflow.onComplete {
    println "Run name : ${workflow.runName}"
    println "Duration : ${workflow.duration}"
    println "Succeeded: ${workflow.success}"
}
