#!/usr/bin/env nextflow

params.reads       = params.reads       ?: "${launchDir}/test/*.fastq.gz"
params.threads     = params.threads     ?: 2
params.adapter     = params.adapter     ?: "AGATCGGAAGAG"
params.minlen      = params.minlen      ?: 18
params.maxlen      = params.maxlen      ?: 27
params.genome      = params.genome      ?: "${launchDir}/genome/caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa"
params.annotation  = params.annotation  ?: "${launchDir}/annotation/caenorhabditis_elegans.PRJNA13758.WBP19.overlapping_annotation.gff3"
params.thr_ss      = params.thr_ss      ?: 12
params.sm_ss       = params.sm_ss       ?: "12G"
params.srcDir = "${workflow.projectDir}/src"

reads_ch = Channel.fromPath(params.reads)
genome_ch = Channel.fromPath(params.genome)
annotation_ch = Channel.fromPath(params.annotation)
summary_script_ch       = Channel.fromPath("${params.srcDir}/01.ss3_summary.R")
bam2Rds_script_ch	= Channel.fromPath("${params.srcDir}/01.bam2Rds.R")
fnmtx_script_ch         = Channel.fromPath("${params.srcDir}/02.get_fn_mtx.R")
featureCounts_script_ch = Channel.fromPath("${params.srcDir}/02.featureCounts.R")
    
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

process shortstack {
    tag "ShortStack analysis"

    input:
    path genome
    path reads
    
    output:
    path "ShortStack_*", emit: ss_out
    path "Log.txt", emit: log
    path "merged.bam", emit: merged_bam

    publishDir "05.ShortStack", mode: 'copy'

    script:
    """
    cp ${genome} genome.fa

    FQFILES=\$(echo *.ps.fastq.gz)

    ShortStack \\
      --genomefile genome.fa \\
      --nohp \\
      --mincov 5 \\
      --pad 1 \\
      --mismatches 1 \\
      --mmap u \\
      --bowtie_m all \\
      --ranmax none \\
      --bowtie_cores ${params.thr_ss} \\
      --sort_mem ${params.sm_ss} \\
      --readfile \$FQFILES
    cp ShortStack_*/Log.txt .
    cp ShortStack_*/merged_alignments.bam merged.bam
    """
}

process summarize_shortstack {
    tag "ShortStack summary plots"

    input:
    path ss_dir
    path log_file
    path summary_script

    output:
    path "ShortStack_summary.pdf", emit: pdf

    publishDir "06.summary", mode: 'copy'

    script:
    """
    Rscript ${summary_script} Log.txt
    """
}

process split_bam {
  tag "Splitting BAM"

  input:
  path bam_file

  output:
  path "*.trim.mapped.bam", emit: split_bams

  publishDir "07.bams", mode: 'copy'

  script:
  """
  samtools split -f %!.%. -@ ${params.threads} -v ${bam_file}

  for file in *.ps.bam; do
      outFile=\$(echo \$file | sed 's|.ps.bam|.trim.mapped.bam|')
      samtools view -F 4 -b -@ ${params.threads} \$file > \$outFile
      rm \$file
  done
  """
}

process bam2Rds {
    tag "${bam_file.simpleName}"
    input:
    path bam2Rds_script
    path bam_file

    output:
    path "*.Rds", emit: matrices

    publishDir "08.rds", mode: 'copy'

    script:
    """
    Rscript ${bam2Rds_script}
    """
}

<<<<<<< HEAD
process fn_mtx {
    tag "Get first-nucleotide matrices"

    input:
    path fnmtx_script
    path rds_files

    output:
    path "first_nt_rlength.Rds", emit: matrices

    publishDir "09.fn_mtx", mode: 'copy'

    script:
    """
    Rscript ${fnmtx_script}
    """
=======
  script:
  """
  sample=\$(basename "${bam_file}" .bam | sed 's/.trim.mapped\$//')

  samtools view "${bam_file}" | \\
  awk '
  function revcomp(seq,    rev, i, base) {
    rev = ""
    for (i = length(seq); i > 0; i--) {
      base = substr(seq, i, 1)
      if (base == "A") base = "T"
      else if (base == "C") base = "G"
      else if (base == "G") base = "C"
      else if (base == "T") base = "A"
      rev = rev base
    }
    return rev
  }
  {
    flag = \$2
    seq = \$10
    if (flag == 16) {
      seq = revcomp(seq)
    }
    len = length(seq)
    first_nt = substr(seq, 1, 1)
    print len, first_nt
  }' > "\${sample}.length_firstnt.txt"
  """
>>>>>>> e0d1cecd1577b4c4976b7d6e7b98d7f2aa0a3986
}

process plot_firstnt {
    tag "Plot first nucleotide distributions"
    input:
    path txt_files

    output:
    path "length_dit_fn_percentage.pdf", emit: pdf_percentage

    publishDir "11.fn_plots", mode: 'copy'

    script:
    """
    Rscript ${params.srcDir}/03.plot_fn.R ${params.minlen} ${params.maxlen}
    """
}


process featureCounts {
    tag "${bam_file.simpleName}"
    input:
    path bam_file
    path annotation
    path featureCounts_script
    
    output:
    path "*.txt", emit: table_of_counts
    path "*.featureCounts", emit: featc_read_tables   
    path "*.RDS", emit: featurecount_rds
    publishDir "09.featureCounts", mode: 'copy'
    
    script:
    """
    Rscript ${featureCounts_script} ${annotation} ${bam_file.join(" ")}
    """
}

process bam2bedgraph {
    tag "${bam_file.simpleName}"

    input:
    path bam_file

    output:
    path "*.bedGraph.gz", emit: bedgraphs

    publishDir "10.bedGraphs", mode: 'copy'

    script:
    """
    for file in *.bam; do
        echo "Processing \$file"
        outFile=\$(echo \$file | sed 's/\\..*//')

        normScale=\$(bc <<< "scale=4;1000000/\$(samtools view -f 0 -c \$file)")

        bedtools genomecov -ibam \$file -bg -scale \$normScale -strand '+' > \$outFile.us

        bedtools genomecov -ibam \$file -bg -scale \$normScale -strand '-' | sed -r 's/([0-9]+)\$/-\\1/' >> \$outFile.us

        bedtools sort -i \$outFile.us | pigz -c > \$outFile.bedGraph.gz
        rm \$outFile.us
    done
    """
}

workflow {
    fastqc_out = fastqc(reads_ch)
    multiqc(fastqc_out.qc_zip.collect())
    trimmed_reads = cutadapt(reads_ch)
    fastq_tr = fastqc_trimm(trimmed_reads)
    multiqc_tr(fastq_tr.qc_zip.collect())
    pulled_reads = pullseq(trimmed_reads)
    (ss_out, log_file, merged_bam) = shortstack(genome_ch, pulled_reads.collect())
    summarize_shortstack(ss_out, log_file, summary_script_ch)
    split_bam_result = split_bam(merged_bam)
    rds = bam2Rds(bam2Rds_script_ch, split_bam_result.split_bams.flatten().collect())
    matrix = fn_mtx(fnmtx_script_ch, rds.matrices.collect())
    plots = plot_firstnt(matrix.matrices)
    featureCounts(split_bam_result.split_bams.collect(),annotation_ch,featureCounts_script_ch)
    bam2bedgraph(split_bam_result.split_bams.flatten())
}
