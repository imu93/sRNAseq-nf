# This script aims to produce Rds files from BAMs
pacman::p_load(rtracklayer, Rsamtools, parallel)

what_to_extract = c("rname", "pos", "qwidth", "strand", "flag", "seq", "qname")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  bams = args
} else {
  # More general than .*mapped.bam$ to catch *.expanded.bam, etc.
  bams = list.files(pattern = "\\.bam$")
}

if (length(bams) == 0) {
  stop("No BAM files found (none passed as arguments and none matched in current directory).")
}

bam2Rds = function(bamf, ...) {
  cat("Processing file:", bamf, "\n")
  example = as.data.frame(
    scanBam(
     bamf,
      param = ScanBamParam(
      what = what_to_extract,
      reverseComplement = TRUE,
      flag = scanBamFlag(isUnmappedQuery = FALSE))
    )
  )

  toGR = GRanges(
    seqnames = example$rname,
    IRanges(start = example$pos, width = example$qwidth),
    strand = example$strand
  )
  toGR$len  = width(toGR)
  toGR$Id   = as.character(example$qname)
  toGR$flag = as.numeric(example$flag)
  toGR$seq  = as.character(example$seq)

  saveRDS(toGR, file = sub("\\.bam$", ".Rds", bamf))
}

mclapply(bams, bam2Rds, mc.cores = 2)
