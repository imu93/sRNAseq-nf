# This script aims to prepare annotation files for sRNAseq-nf
pacman::p_load(rtracklayer, Rsamtools, dplyr, stringr, optparse, GenomicRanges)

option_list = list(
  make_option(c("-G", "--genome_file"), type="character", default=NULL,
              help="Genome FASTA file", metavar="character"),
  make_option(c("-P", "--pseudogene_inc"), type="logical", default=FALSE,
              help="Whether pseudogenes are provided [default: %default]"),
  make_option(c("-g","--gene_annotation"), type="character", default=NULL,
              help="GFF/GFF3 with protein-coding gene annotation", metavar="character"),
  
  make_option(c("-y", "--pseudogene_gff"), type="character", default=NULL,
              help="GFF/GFF3 with pseudogene annotation (optional)", metavar="character"),
  
  make_option(c("--rRNA"), type="character", default=NULL,
              help="GFF/GFF3 with rRNA annotation (optional)", metavar="character"),
  
  make_option(c("--tRNA"), type="character", default=NULL,
              help="GFF/GFF3 with tRNA annotation (optional)", metavar="character"),
  
  make_option(c("--miRNA"), type="character", default=NULL,
              help="GFF/GFF3 with miRNA annotation (optional)", metavar="character"),
  
  make_option(c("--piRNA"), type="character", default=NULL,
              help="GFF/GFF3 with piRNA annotation (optional)", metavar="character"),
  
  make_option(c("--lncRNA"), type="character", default=NULL,
              help="GFF/GFF3 with lncRNA annotation (optional)", metavar="character"),
  
  make_option(c("--lincRNA"), type="character", default=NULL,
              help="GFF/GFF3 with lincRNA annotation (optional)", metavar="character"),
  
  make_option(c("--snRNA"), type="character", default=NULL,
              help="GFF/GFF3 with snRNA annotation (optional)", metavar="character"),
  
  make_option(c("--snoRNA"), type="character", default=NULL,
              help="GFF/GFF3 with snoRNA annotation (optional)", metavar="character"),
  
  make_option(c("--circRNA"), type="character", default=NULL,
              help="GFF/GFF3 with circRNA annotation (optional)", metavar="character"),
  
  make_option(c("--other_ncRNA"), type="character", default=NULL,
              help="GFF/GFF3 with other_ncRNA annotation (optional)", metavar="character"),
  
  make_option(c("--sRNA"), type="character", default=NULL,
              help="GFF/GFF3 with sRNA annotation (optional)", metavar="character"),
  
  make_option(c("-R", "--repeats"), type="character", default=NULL,
              help="GFF/GFF3 with repeat annotation (optional)", metavar="character"),
  
  make_option(c("-u", "--uns_out_gff3"), type="character", default="uns_annotation.gff3",
              help="Output merged GFF3 [default: %default]", metavar="character"),
  
  make_option(c("-o", "--out_gff3"), type="character", default="annotation.gff3",
              help="Output merged GFF3 [default: %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

die = function(...) stop(paste0(...), call.=FALSE)
check_file = function(p, label) {
  if (is.null(p) || !nzchar(p)) die("Missing required input: ", label)
  if (!file.exists(p)) die("File not found for ", label, ": ", p)
  if (file.info(p)$size == 0) die("File is empty for ", label, ": ", p)
  normalizePath(p)
}

read_gff = function(p, label) {
  if (is.null(p) || !nzchar(p)) return(NULL)
  if (!file.exists(p)) die("File not found for ", label, ": ", p)
  if (file.info(p)$size == 0) die("File is empty for ", label, ": ", p)
  message("[load] ", label, ": ", p)
  import(p)
}

faFile = check_file(opt$genome_file, "--genome_file")

if (isTRUE(opt$pseudogene_inc)) {
  if (is.null(opt$pseudogene_gff)) {
    die("pseudogene_inc=TRUE but no pseudogene inputs were provided (--pseudogene_gff)")
  }
}

is_valid_gr = function(gr) {
  !is.null(gr) && inherits(gr, "GenomicRanges") && length(gr) > 0
}

gff_gene   = read_gff(opt$gene_annotation, "gene_annotation")
gff_psg    = read_gff(opt$pseudogene_gff, "pseudogene_gff")
gff_rrna   = read_gff(opt$rRNA, "rRNA")
gff_trna   = read_gff(opt$tRNA, "tRNA")
gff_mirna  = read_gff(opt$miRNA, "miRNA")
gff_pirna  = read_gff(opt$piRNA, "piRNA")
gff_lncrna = read_gff(opt$lncRNA, "lncRNA")
gff_lincr  = read_gff(opt$lincRNA, "lincRNA")
gff_snrna  = read_gff(opt$snRNA, "snRNA")
gff_snorna = read_gff(opt$snoRNA, "snoRNA")
gff_circ   = read_gff(opt$circRNA, "circRNA")
gff_sRNA   = read_gff(opt$sRNA, "sRNA")
gff_other  = read_gff(opt$other_ncRNA, "other_ncRNA")
gff_rep    = read_gff(opt$repeats, "repeats")


gffList = list(
  gene        = gff_gene,
  pseudogene  = gff_psg,
  rRNA        = gff_rrna,
  tRNA        = gff_trna,
  miRNA       = gff_mirna,
  piRNA       = gff_pirna,
  lncRNA      = gff_lncrna,
  lincRNA     = gff_lincr,
  snRNA       = gff_snrna,
  snoRNA      = gff_snorna,
  circRNA     = gff_circ,
  sRNA        = gff_sRNA,
  other_ncRNA = gff_other,
  repeats      = gff_rep
)

gffList = gffList[vapply(gffList, is_valid_gr, logical(1))]

if (length(gffList) == 0) die("No annotation inputs contained valid instances.")


sInfo = seqinfo(scanFaIndex(faFile))
gen_lv = seqlevels(sInfo)
allFeatures = list()
for (file in names(gffList)) {
  print(file)
  x = gffList[[file]]
  x = as.data.frame(x)
  x = x[!grepl("^(MT|chrM)$", x$seqnames),]
  x = makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  
  gff_lv = seqlevels(x)
  common = intersect(gff_lv, gen_lv)
  only_gff = setdiff(gff_lv, gen_lv)
  only_gen = setdiff(gen_lv, gff_lv)
  
  if (length(common) == 0) {
    message("[FAIL] ", file, " has 0 shared seqlevels with genome")
    message("only in GFF: ", paste(gff_lv, collapse=", "))
    die("No shared seqlevels for: ", file)
  }
  
  if (length(only_gff) > 0) {
    message("[MISMATCH] ", file, " only in GFF: ", paste(only_gff, collapse=", "))
  }
  
  x = keepSeqlevels(x, common, pruning.mode="coarse")
  seqinfo(x) = sInfo[seqlevels(x)]
  
  allFeatures[[file]] = x
}


names(allFeatures) = names(gffList)

if ("gene" %in% names(allFeatures)) {
  allFeatures$gene = allFeatures$gene[grepl("cds|CDS|UTR|utr", allFeatures$gene$type),]
}

if ("other_ncRNA" %in% names(allFeatures)) {
  allFeatures$other_ncRNA = allFeatures$other_ncRNA[allFeatures$other_ncRNA$type == "gene",]
}

if ("pseudogene" %in% names(allFeatures)) {
  pseudo_transcripts = allFeatures$pseudogene[allFeatures$pseudogene$type == "pseudogenic_transcript",]
  allFeatures$pseudogene = allFeatures$pseudogene[allFeatures$pseudogene$type == "exon",]
  allFeatures$pseudogene$ID = ifelse(
    is.na(allFeatures$pseudogene$ID),
    sub(".*:", "gene=", unlist(allFeatures$pseudogene$Parent)),
    allFeatures$pseudogene$ID
  )
}

ensure_id = function(gr, prefix) {
  id = mcols(gr)$ID
  if (is.null(id) || length(id) == 0) {
    mcols(gr)$ID = paste0(prefix, ":", as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  } else {
    bad = is.na(id) | !nzchar(id)
    if (any(bad)) id[bad] = paste0(prefix, ":", as.character(seqnames(gr))[bad], ":", start(gr)[bad], "-", end(gr)[bad])
    mcols(gr)$ID = id
  }
  gr
}

allFeatures = mapply(function(n, gr) ensure_id(gr, n), names(allFeatures), allFeatures, SIMPLIFY=FALSE)

str_list = list()
for (feature in names(allFeatures)) {
  tmp.feat = allFeatures[[feature]]
  tmp.feat$type = ifelse(feature != "repeats", feature, as.character(tmp.feat$type))
  
  if (feature == "repeats") {
    if (!is.null(mcols(tmp.feat)$class)) {
      tmp.feat$type = as.character(tmp.feat$class)
    } else {
      tmp.feat$type = as.character(tmp.feat$type)
    }
    tmp.feat$type = str_replace(tmp.feat$type, "_S|_As", "")
    strand(tmp.feat) = ifelse(strand(tmp.feat) == "*", "+", strand(tmp.feat))
  }
  
  tmp.feat$ID  = str_replace(tmp.feat$ID, "gene=", "Gene:")
  tmp.feat$ID  = str_replace(tmp.feat$ID, "Gene:", paste0(tmp.feat$type, ":"))
  tmp.feat$ID  = str_replace(tmp.feat$ID, ":", "_S:")
  
  tmp.AS = tmp.feat
  strand(tmp.AS) = ifelse(strand(tmp.AS) == "+", "-", "+")
  tmp.AS$ID = str_replace(tmp.AS$ID, "_S:", "_As:")
  
  str_list[[feature]] = c(tmp.feat, tmp.AS)
}
names(str_list) = NULL

str_annotation = do.call(c, str_list)
str_annotation = str_annotation[order(seqnames(str_annotation), start(str_annotation)),]



str_annotation$IGV = paste0(as.character(seqnames(str_annotation)), ":", start(str_annotation), "-", end(str_annotation))
str_annotation$color = "grey"

col2 =  c("#0B3BFA","#007CFF", "#73C2FF", "#00C3DB",
          "#00AB14",   "#90DB3F", "#D2F102", "#FFDA2C",
          "#F02765", "#F75616","#8846A3","#601573", "#807B77", "#DBDBDB")


str_annotation$color = ifelse(grepl("^gene", str_annotation$type), col2[1], str_annotation$color)
str_annotation$color = ifelse(grepl("^pseudo", str_annotation$type), col2[3], str_annotation$color)
str_annotation$color = ifelse(grepl("DNA|RC", str_annotation$type), col2[4], str_annotation$color)
str_annotation$color = ifelse(grepl("LTR", str_annotation$type), col2[5], str_annotation$color)
str_annotation$color = ifelse(grepl("LINE|PLE|Retrop", str_annotation$type), col2[6], str_annotation$color)
str_annotation$color = ifelse(grepl("Unknown", str_annotation$type), col2[7], str_annotation$color)
str_annotation$color = ifelse(grepl("miRNA", str_annotation$type), col2[8], str_annotation$color)
str_annotation$color = ifelse(grepl("rRNA", str_annotation$type), col2[9], str_annotation$color)
str_annotation$color = ifelse(grepl("tRNA", str_annotation$type), col2[10], str_annotation$color)
str_annotation$color = ifelse(grepl("lincRNA", str_annotation$type), col2[11], str_annotation$color)
str_annotation$color = ifelse(grepl("lncRNA|piRNA|snRNA|snoRNA|yRNA|other_ncRNA|circRNA", str_annotation$type), col2[12], str_annotation$color)
str_annotation$color = ifelse(grepl("Low|Simple|Sate|SINE", str_annotation$type), col2[13], str_annotation$color)

uns_annotation = str_annotation[grepl("_S", str_annotation$ID),]
uns_annotation$Parent = uns_annotation$ID

if (!is.null(mcols(uns_annotation)$class)) uns_annotation$Class = uns_annotation$class
if (is.null(mcols(uns_annotation)$biotype)) uns_annotation$biotype = NA_character_

uns_annotation$biotype = ifelse(grepl("gene_S", uns_annotation$ID), "gene", uns_annotation$biotype)
uns_annotation$biotype = ifelse(grepl("DNA|RC", uns_annotation$ID), "Transposon", uns_annotation$biotype)
uns_annotation$biotype = ifelse(grepl("LTR|PLE|LINE|SINE|Retro", uns_annotation$ID), "Retrotransposon", uns_annotation$biotype)
uns_annotation$biotype = ifelse(grepl("Unk|Sat|Simple|Low", uns_annotation$ID), "Other_repeat", uns_annotation$biotype)
uns_annotation$biotype = ifelse(grepl("circRNA", uns_annotation$type), "circRNA", uns_annotation$biotype)
uns_annotation$biotype = ifelse(grepl("lncRNA", uns_annotation$type), "lncRNA", uns_annotation$biotype)

export(uns_annotation, opt$uns_out_gff3, "gff3")
export(str_annotation, opt$out_gff3, "gff3")
