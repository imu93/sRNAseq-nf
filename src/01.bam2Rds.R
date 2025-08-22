# This script aims to produce Rds files from bams
# Load packages
pacman::p_load(rtracklayer, Rsamtools, parallel)
# Define columns to extract from bam files
what_to_extract = c("rname", "pos", "qwidth", "strand", "flag", "seq", "qname")
# Define input files
bams = list.files(pattern = ".*mapped.bam$")

# Check if any BAM files are present
if (length(bams) == 0) {
  stop("No BAM files found in the current directory")
}

# Create a function to produce files
bam2Rds = function(bamf, ...) {
  cat("Processing file:", bamf, "\n")
  print(bamf)
  # Extract columns as df with Rsamtools
  example = as.data.frame(scanBam(bamf, 
                                  param = ScanBamParam(what = what_to_extract,
                                                       reverseComplement = TRUE,
                                                       flag = scanBamFlag(isUnmappedQuery = FALSE))))
  # Now as GRanges
  toGR = GRanges(seqnames = example$rname, 
                 IRanges(start = example$pos, width = example$qwidth),
                 strand = example$strand)
  # Add metadata columns
  toGR$len = width(toGR)
  toGR$Id = as.character(example$qname)
  toGR$flag = as.numeric(example$flag)
  toGR$seq = as.character(example$seq)
  # Save as RDS
  saveRDS(toGR, file = gsub(".bam", ".Rds", bamf))
}
# run
mclapply(bams, bam2Rds, mc.cores = 2)
