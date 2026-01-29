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


// Note: Include the scripts for feature + fisrt nt and test


// NEW params for annotations chkp
params.cfeat	=	params.cfeat	?:	"rRNA,tRNA,gene,DNA,LINE,LTR"
params.mfeat	=	params.mfeat	?:	1
params.reads       = params.reads       ?: "${launchDir}/test/*.fastq.gz"
params.threads     = params.threads     ?: 2
params.preproc     = params.preproc     ?: 'legacy'   // 'fastp' or 'legacy'
params.adapter     = params.adapter     ?: "AGATCGGAAGAG"
params.fastp_use_adapter = (params.fastp_use_adapter != null ? params.fastp_use_adapter : false)
params.disable_length_filter = (params.disable_length_filter != null   ? params.disable_length_filter.toString().toLowerCase() in ['true','1','yes','y'] : false)
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


params.raw_mode = params.raw_mode ?: false
// Only for benchmaks
params.uwm_mmap_max = params.uwm_mmap_max ?: null
if (params.uwm_mmap_max != null) {
    try {
        int value = params.uwm_mmap_max as int
        if (value < 1) {
            params.uwm_mmap_max = null
        } else {
            params.uwm_mmap_max = value
        }
    } catch (Exception e) {
        params.uwm_mmap_max = null
    }
}

params.uwm_suppress_equal_prob = (params.uwm_suppress_equal_prob != null
    ? params.uwm_suppress_equal_prob.toString().toLowerCase() in ['true','1','yes','y']
    : false)


params.make_bedgraph = (params.make_bedgraph != null ? params.make_bedgraph : true) ?: null
if (params.make_bedgraph instanceof String) {
    params.make_bedgraph = params.make_bedgraph.toLowerCase() in ['true','1','yes','y']
}

params.expand_unmapped = params.expand_unmapped ?: false //boolean
params.n_mismatch    = params.n_mismatch    ?: 1 // integer value (bowtie allows up to 3 mismatching bases)
params.index_thr     = params.index_thr ?: 1 // the max number will be the number of libraries
params.offrate_sm  = params.offrate_sm  ?: 4
params.thr_sm      = params.thr_sm      ?: 12
params.save_non_mappers = (
    params.save_non_mappers != null
        ? params.save_non_mappers.toString().toLowerCase() in ['true','1','yes','y']
        : true   // default: mismo comportamiento que antes (sí los guarda)
)
params.smem_sm     = params.smem_sm     ?: "12G"
params.wins_sm     = params.wins_sm     ?: 250
params.consider_strand = params.consider_strand ?: false
params.assign_mode    = params.assign_mode    ?: 'uwm'   // 'uwm' or 'random'
params.minoverlap     = params.minoverlap     ?: 0.7
params.contrast       = params.contrast       ?: "${launchDir}/contrast.txt"

// and the same for the contrast
if( !params.contrast ) {
  error "\n[INPUT ERROR] params.contrast is not set\n"
}

def contrast_file = file(params.contrast)

if( !contrast_file.exists() ) {
  error "\n[INPUT ERROR] Contrast file not found: ${contrast_file}\n"
}

if( contrast_file.size() == 0 ) {
  error "\n[INPUT ERROR] Contrast file is empty: ${contrast_file}\n"
}
params.fn_reps        = (params.fn_reps != null ? params.fn_reps.toString().toLowerCase() in ['true','1','yes','y'] : false)
params.fn_all_lengths = (params.fn_all_lengths != null ? params.fn_all_lengths.toString().toLowerCase() in ['true','1','yes','y'] : false)
params.lfc            = params.lfc            ?: 1
params.fdr            = params.fdr            ?: 0.05
params.threshold_inc   = params.threshold_inc   ?: false   // read FC/FDR from contrast file?
params.hk_norm	      = params.hk_norm        ?: false // boolean
params.norm_feature   = params.norm_feature   ?: "rRNA_S"
params.disable_tmm_for_fbn = params.disable_tmm_for_fbn  ?: false   // if true: use all features instead of TMM or pseudo-TMM
params.stringent_tmm  = params.stringent_tmm  ?: false   // boolean
params.use_rds        = params.use_rds        ?: true   // boolean
params.first_nt       = params.first_nt       ?: ""
params.make_complex_plots = (params.make_complex_plots != null ? params.make_complex_plots.toString().toLowerCase() in ['true','1','yes','y'] : false)
params.require_assigned = (params.require_assigned != null ? params.require_assigned.toString().toLowerCase() in ['true','1','yes','y'] : false)
params.apply_first_nt_downstream = params.apply_first_nt_downstream ?: false
params.srcDir         = "${workflow.projectDir}/bin"
reads_ch              = Channel.fromPath(params.reads)
genome_ch             = Channel.value( file(params.genome) )
annotation_ch         = Channel.fromPath(params.annotation)
contrast_ch           = Channel.fromPath(params.contrast)
siRmap_script_ch      = Channel.value( file("${params.srcDir}/siRmap.py") )
collapse_script_ch    = Channel.value( file("${params.srcDir}/collapse") )
ann_chkp_ch	 	= Channel.fromPath("${params.srcDir}/00.annotation_checkp.R")
summary_script_uwm_ch   = Channel.fromPath("${params.srcDir}/02.summary_uwm.R")
summary_script_rand_ch  = Channel.fromPath("${params.srcDir}/02.summary_rand.R")
bam2Rds_script_ch	= Channel.fromPath("${params.srcDir}/01.bam2Rds.R")
fnmtx_script_ch         = Channel.fromPath("${params.srcDir}/02.get_fn_mtx.R")
featureCounts_script_ch = Channel.fromPath("${params.srcDir}/02.featureCounts.R")
fcComplex_script_ch     = Channel.fromPath("${params.srcDir}/03.fcComplex.R")
cls_mtx_script_ch       = Channel.fromPath("${params.srcDir}/03.cls_mtx.py")
fn_cls_script_ch        = Channel.fromPath("${params.srcDir}/03.complex_plot.R")
dea_script_ch           = Channel.fromPath("${params.srcDir}/04.edgeR.R")



def results_dir = "Results_${new Date().format('yyyyMMdd_HHmmss')}"
new File(results_dir).mkdirs()

process validate_annotation {
  tag "validate_annotation"
  label 'preproc'

  input:
    path ann_chkp
    path annotation

  output:
    path "annotation_summary.txt", emit: report
    path "annotation.validation.ok",  emit: ok
  
  
  publishDir "${results_dir}/00.ann_rep", mode: 'copy'
    
    
  script:
  """
  Rscript ${ann_chkp} -i ${annotation} \
  		      -c "${params.cfeat}" \
  		      -m ${params.mfeat} \
  		      -o annotation_summary.txt

  touch annotation.validation.ok
  """
}

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

  LEN_OPT=""
if [ "${params.disable_length_filter}" != "true" ]; then
  LEN_OPT="--length_required ${params.minlen} --length_limit ${params.maxlen}"
fi

  fastp \
    -i ${read} \
    -o "\${base}.ps.fastq.gz" \
    \$LEN_OPT \
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

process cutadapt {
  label 'preproc'

  input:
    path read

  output:
    path "*trimmed.fastq.gz", emit: fastq

  publishDir "${results_dir}/02.cut_adapt", mode: 'copy'
  tag "${read.simpleName}"

  script:
  """
  MINLEN_OPT=""
  if [ "${params.disable_length_filter}" = "true" ]; then
    # minimo tecnico para evitar reads vacíos (len=0) que crashean bowtie
    MINLEN_OPT="-m 10"
  else
    # filtro normal (biologico)
    MINLEN_OPT="-m ${params.minlen}"
  fi

  cutadapt -j ${params.threads} \
           -a ${params.adapter} \
           \$MINLEN_OPT \
           --max-n 0.05 \
           --discard-untrimmed \
           -o ${read.simpleName}.trimmed.fastq.gz \
           ${read}
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

process pullseq_raw {
  label 'preproc'

  input:
    path read

  output:
    path "*ps.fastq.gz", emit: fastq

  tag "${read.simpleName}"
  publishDir "${results_dir}/04.pullseq_raw", mode: 'copy'

  script:
  """
  base='${read.simpleName}'
  base=\${base%.fastq}
  base=\${base%.fastq.gz}
  base=\${base%.fq}
  base=\${base%.fq.gz}

 
  if [ "${params.disable_length_filter}" = "true" ]; then
    read MIN MAX < <(
      zcat ${read} | awk 'NR%4==2{l=length(\$0); if(min==""||l<min)min=l; if(max==""||l>max)max=l} END{print min, max}'
    )
    if [ -z "\$MIN" ] || [ -z "\$MAX" ]; then
      echo "ERROR: could not infer read lengths" >&2
      exit 1
    fi
    # pullseq requiere MAX > MIN
    if [ "\$MAX" -le "\$MIN" ]; then
      MAX=\$((MIN+1))
    fi
  else
    MIN=${params.minlen}
    MAX=${params.maxlen}
  fi

  # --- SIEMPRE darle a pullseq un FASTQ plano
  zcat ${read} > "\${base}.in.fastq"
  pullseq -i "\${base}.in.fastq" -m \$MIN -a \$MAX > "\${base}.ps.fastq"
  rm -f "\${base}.in.fastq"
  pigz "\${base}.ps.fastq"
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
  base=\${base%.fastq.gz}
  base=\${base%.fq}
  base=\${base%.fq.gz}

 
  if [ "${params.disable_length_filter}" = "true" ]; then
    read MIN MAX < <(
      zcat ${read} | awk 'NR%4==2{l=length(\$0); if(min==""||l<min)min=l; if(max==""||l>max)max=l} END{print min, max}'
    )
    if [ -z "\$MIN" ] || [ -z "\$MAX" ]; then
      echo "ERROR: could not infer read lengths" >&2
      exit 1
    fi
    # pullseq requiere MAX > MIN
    if [ "\$MAX" -le "\$MIN" ]; then
      MAX=\$((MIN+1))
    fi
  else
    MIN=${params.minlen}
    MAX=${params.maxlen}
  fi

  # --- SIEMPRE darle a pullseq un FASTQ plano
  zcat ${read} > "\${base}.in.fastq"
  pullseq -i "\${base}.in.fastq" -m \$MIN -a \$MAX > "\${base}.ps.fastq"
  rm -f "\${base}.in.fastq"
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

process map_raw {
  tag { fastq.baseName }
  cpus params.thr_sm
  memory params.smem_sm

  input:
    path fastq
    path genome
    path siRmap_script
    path bowtie_idx         

  output:
    path "*.bam",               emit: raw_bam
    path "*.unmapped.fastq.gz", emit: unmapped_fq, optional: true
    path "*.log",               emit: logs

  publishDir "${results_dir}/05.map_raw", mode: 'copy'

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
  BASENAME=\${BASENAME%.ps}

  python ${siRmap_script} bowtie-aln \
    --fastq ${fastq} \
    --index "\$IDX_BASE" \
    --mismatches ${params.n_mismatch} \
    --out "\${BASENAME}.bam" \
    --save-np ${params.save_non_mappers ? 'yes' : 'no'}  \
    --non-mappers "\${BASENAME}.unmapped.fastq.gz" \
    --threads ${params.thr_sm} \
    --sort-mem ${params.smem_sm} \
    ${params.max_multimaps ? "--max-multimaps ${params.max_multimaps}" : ""}
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
    path "*.unmapped.collapsed.fastq.gz", emit: unmapped_fq, optional: true
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
    --save-np ${params.save_non_mappers ? 'yes' : 'no'}  \
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
    --do-sort true \
    ${ params.raw_mode ? '--raw' : '' }
  """
}

process build_unique_index {
  tag "unique-index"
  cpus params.index_thr
  memory params.smem_sm


  input:
    path siRmap_script
    path bam_paths
    

  output:
    path "unique_index.pkl", emit: uniq_idx


  publishDir "${results_dir}/06.uwm_index", mode: 'copy'
  script:
  """
  printf "%s\n" ${bam_paths.join(' ')} | tr ' ' '\\n' > bam_list.txt
  python ${siRmap_script} build-index-uwm \
                          --bams bam_list.txt \
                          --threads ${params.index_thr} \
                          --out unique_index.pkl 
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
    --seed 123 ${ params.consider_strand ? '--strand' : '' }  \
    ${ params.raw_mode ? '--raw' : '' } \
    ${ (params.raw_mode && params.uwm_mmap_max != null) ? "--mmap-max ${params.uwm_mmap_max}" : "" } \
    ${ (params.raw_mode && params.uwm_suppress_equal_prob) ? "--suppress-equal-prob" : "" }
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
  samtools view -h -F 0x904 -@ ${params.threads} "${bam}" \\
| awk -v DISABLE="${params.disable_length_filter}" -v MIN=${params.minlen} -v MAX=${params.maxlen} '
    function comp(b){ return (b=="A"?"T": b=="C"?"G": b=="G"?"C": b=="T"?"A": b) }
    /^@/ { next }
    {
      seq=toupper(\$10); l=length(seq)
      if (DISABLE!="true" && (l<MIN || l>MAX)) next

      b=substr(seq,1,1)
      if (and(\$2,16)) b=comp(substr(seq,l,1))

      if (b=="A"||b=="C"||b=="G"||b=="T") { cnt[l ":" b]++; tot[l]++ }
    }
    END{
      print "sample\\tlength\\tA\\tC\\tG\\tT\\ttotal"
      # gather lengths
      for (x in tot) lens[x]=1
      n=0
      for (x in lens) { L[++n]=x }
      # sort numeric
      for (i=1;i<=n;i++) for (j=i+1;j<=n;j++) if (L[i]+0 > L[j]+0) { tmp=L[i]; L[i]=L[j]; L[j]=tmp }
      for (i=1;i<=n;i++){
        l=L[i]
        a=cnt[l ":A"]+0; c=cnt[l ":C"]+0; g=cnt[l ":G"]+0; t=cnt[l ":T"]+0
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
  label 'preproc'

  input:
    path tables

  output:
    path "length_dit_fn_percentage.png",       emit: avg_png
    path "length_dit_fn_percentage.reps.png",  emit: reps_png, optional: true

  publishDir "${results_dir}/09.fn_plots", mode: 'copy'

  script:
  def repsFlag = params.fn_reps ? '--reps TRUE' : ''
  def allFlag  = params.fn_all_lengths ? '--all_lengths TRUE' : ''
  def minFlag  = "--min ${params.minlen}"
  def maxFlag  = "--max ${params.maxlen}"

  """
  Rscript ${params.srcDir}/03.plot_fn.R \
    ${minFlag} ${maxFlag} \
    ${repsFlag} \
    ${allFlag}
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


process featureCounts_tagBam {
  tag "featureCounts (tagged BAM)"
  when:
    params.make_complex_plots

  input:
    path bam_files
    path annotation
    path fcComplex_script

  output:
    path "*.featureCounts.bam", emit: tagged_bams
    path "*.log", optional: true

  publishDir "${results_dir}/11.featureCounts_tagged", mode: 'copy'

  script:
  """
  Rscript ${fcComplex_script} ${annotation} ${params.minoverlap} ${params.thr_sm} ${bam_files.join(" ")}
  """
}

process cls_mtx {
  tag "Class matrix (cls_mtx)"
  when:
    params.make_complex_plots

  input:
    path bams
    path cls_mtx_script

  output:
    path "*.cls_mtx.tsv", emit: cls_tables

  publishDir "${results_dir}/11.featureCounts_tagged", mode: 'copy'

  script:
    def lenArgs     = params.disable_length_filter ? "--all-lengths" : "--minlen ${params.minlen} --maxlen ${params.maxlen}"
    def assignedArg = params.require_assigned ? "--require-assigned" : ""
    """
    python ${cls_mtx_script} \
      --bam ${bams.join(' ')} \
      --suffix ".cls_mtx.tsv" \
      ${lenArgs} \
      ${assignedArg}
    """
}

process fn_cls{
  tag "Complex plot (fn_cls)"

  when:
    params.make_complex_plots

  input:
    tuple val(id), path(fn_tables), path(cls_tables)

  output:
    path "fn_cls_cpm.png", emit: avg_png

  publishDir "${results_dir}/09.fn_plots", mode: 'copy'

  script:
  """
  Rscript ${params.srcDir}/03.complex_plot.R ${fn_tables.join(" ")} ${cls_tables.join(" ")}
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
    def fbnFlag = params.disable_tmm_for_fbn ? 'TRUE' : 'FALSE'
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
      --stringent_tmm ${STRICT} \
      --disable_tmm_for_fbn ${fbnFlag}
      
    """
}

workflow {
  
  validate_annotation(ann_chkp_ch, annotation_ch)
   
  def map_inputs_ch
  def rc_map_ch
  if ( params.preproc == 'fastp' ) {
   
    fp_out = fastp(reads_ch)
    multiqc_fastp( fp_out.qc_json.collect() )

    if( params.raw_mode ) {
      map_inputs_ch = fp_out.fastq
      rc_map_ch     = Channel.empty()
    } else {
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
    }

  } else {
     
    fastqc_raw = fastqc(reads_ch)
    multiqc( fastqc_raw.qc_zip.collect() )

    trimmed   = cutadapt(reads_ch)
    fastqc_tr = fastqc_trimm(trimmed)
    multiqc_tr( fastqc_tr.qc_zip.collect() )

    if( params.raw_mode ) {
      def pulled_raw = pullseq_raw(trimmed)
      map_inputs_ch = pulled_raw.fastq
      rc_map_ch     = Channel.empty()
    } else {
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

  def mapped
  if( params.raw_mode ) {
    mapped = map_raw(
              fastq_map_for_map,
              genome_ch,
              siRmap_script_ch,
              idx.ebwt.collect()
            )
  } else {
    mapped = map_collapse(
              fastq_map_for_map,
              genome_ch,
              siRmap_script_ch,
              idx.ebwt.collect()
            )
  }

  def mapped_bam_ch = params.raw_mode ? mapped.raw_bam : mapped.collapsed_bam

  // new if to expand unpaed reads by ID in case user need this
  if( params.expand_unmapped && !params.raw_mode && params.save_non_mappers ) {
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
                    mapped_bam_ch.collect()
                  ).uniq_idx

    def uwm_out = resolve_uwm(
                    mapped_bam_ch,
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
                     mapped_bam_ch,
                     siRmap_script_ch
                   )

    resolved_bams_ch = rand_out.expanded_bam

    summarize_rand(
      rand_out.logs.collect(),
      summary_script_rand_ch
    )

  } else {
    resolved_bams_ch = mapped_bam_ch
  }
  
  if( params.use_rds ) {
    rds_out   = bam2Rds(bam2Rds_script_ch, resolved_bams_ch.collect())
  }

  def counts_out = firstnt_counts(resolved_bams_ch)
    def all_counts = counts_out.counts.collect()
    plot_firstnt(all_counts)
    def all_fn = all_counts.map { files -> tuple('all', files) }
  

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

  
  if (params.make_complex_plots) {
    fc_bam = featureCounts_tagBam(bams_for_quant.collect(), annotation_ch, fcComplex_script_ch)

    def cls_out = cls_mtx(fc_bam.tagged_bams.collect(), cls_mtx_script_ch)

    def all_cls = cls_out.cls_tables.collect()
                  .map { files -> tuple('all', files) }

    fn_cls( all_fn.join(all_cls) )
  }


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
