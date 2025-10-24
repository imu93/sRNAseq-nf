# This script aims to produce fn matrices
# Load packages
pacman::p_load(rtracklayer,parallel)
args = commandArgs(trailingOnly = T)

# Read all Rds files produced with 01.bam2Rds.R
files = list.files(pattern = ".*.Rds$")

# Now I need a function to obtain fn
get_fn_mtx = function(file, ...){
  cat("Processing file:", file, "\n")
  x = readRDS(file)
  lib  = subset(x, seqnames(x) != "MtDNA") # 
  fnuc = substr(lib$seq,1,1)
  lens = split(fnuc, width(lib))
  lis = lapply(lens,table)
  lis = lapply(lis, function(x){x[c("C","T","A","G")]})
  counts = as.data.frame(do.call(rbind, lis))
  return(counts)
}

# Run over all the replicates
fn_mtx = mclapply(files, get_fn_mtx, mc.cores = 1)
names(fn_mtx) = files
# Rave the list
saveRDS(fn_mtx, "first_nt_rlength.Rds")
