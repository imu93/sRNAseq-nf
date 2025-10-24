#!/usr/bin/env Rscript
pacman::p_load(Rsubread, rtracklayer)

args = commandArgs(trailingOnly = TRUE)
gff_file = args[1]
min_ov = as.numeric(args[2])
threads = as.numeric(args[3])
bam_files = args[4:length(args)]


# Import GFF3
genome = import(gff_file)

# Build SAF data frame
df = data.frame(
  GeneID = genome$ID,
  Chr = as.character(seqnames(genome)),
  Start = start(genome),
  End = end(genome),
  Strand = as.character(strand(genome))
)

# Run featureCounts
tb = featureCounts(
  files = bam_files,
  annot.ext = df,
  allowMultiOverlap = TRUE,
  strandSpecific = 1,
  fracOverlap = min_ov,
  nthreads = threads,
  useMetaFeatures = TRUE,
  largestOverlap = TRUE,
  reportReads = "CORE"
)

# Save output
saveRDS(tb, "featureCounts_results.RDS")
write.table(tb$counts, "featureCounts_counts.txt", sep = "\t", quote = FALSE)
