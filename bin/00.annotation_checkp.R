pacman::p_load(rtracklayer, dplyr, tidyr,stringr,optparse)

option_list = list(
  make_option(c("-i", "--gff_file"), type="character", default=NULL,
              help= "Path to GFF file (supported formats GFF and GFF3)", metavar="character"),
  make_option(c("-c", "--required_categories"), type="character", default="rRNA,tRNA,gene,DNA,LINE,LTR",
               help = "Categories required in the GFF/GFF3 annotation file", metavar = "character"),
  make_option(c("-m", "--min_instances"), type = "numeric", default = 1,
              help = "Min number of genomic instances per feature", metavar = "numeric"),
  make_option(c("-o", "--out_txt_file"), type="character", default="annotation_summary.txt",
              help= "Path to out txt table with summary", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



if (is.null(opt$gff_file)){
  print_help(opt_parser)
  stop("A GFF file must be provided", call.=FALSE)
}

minFeat = opt$min_instances
cats = opt$required_categories
gffFile = opt$gff_file
outFile = opt$out_txt_file
################################################################################
#gffFile = "~/storage/Data/CELE_ANN/sirmap/sRNAseq-nf/example/annotation/caenorhabditis_elegans.PRJNA13758.WBP19.overlapping_annotation.gff3.gz"
gff = import(gffFile)
n_instances = length(gff)
lst_gff = split(gff, gff$type)

strand = sapply(lst_gff, function(x){grepl("_S|_As", x$ID) %>% unique()})
instances = sapply(lst_gff, length)
wdt_ar = sapply(lst_gff, function(x){width(x) %>% sum}) / 1e6 
wdt_rd = sapply(lst_gff, function(x){reduce(x) %>% width() %>% sum}) / 1e6
df = tibble("Class"=names(strand), "N_instances"=instances, "Mb"=wdt_ar, "NOV_Mb"=wdt_rd, "Bi_stranded"=strand)

c_classes = opt$required_categories %>% strsplit(., ",") %>% unlist()
#c_classes = "rRNA,tRNA,gene,DNA,LINE,LTR" %>% strsplit(.,",") %>% unlist()
na_clas = setdiff(c_classes, df$Class)

if (length(na_clas) >= 1) {
  print(paste("Error:", na_clas, "is missing in the annotation file"))
}

comp = df[match(c_classes, df$Class),]
#minFeat =1

if (!all(comp$N_instances >= minFeat)) {
  stop(paste("Not all features have at least:", minFeat, "instance"))
       
}

if (is.null(mcols(gff)$IGV)) {
  stop("GFF is missing required attribute column: 'IGV'", call.=FALSE)
}

igv = as.character(mcols(gff)$IGV)
if (all(is.na(igv)) || all(!nzchar(igv))) {
  stop("GFF column 'IGV' exists but is empty for all records", call.=FALSE)
}

igv = as.character(mcols(gff)$IGV)
str = sub(":.*", "", gff$ID) %>% sub(".*_", "", .)
id = paste(igv, str, sep= "_")

# classify each row
tag <- case_when(
  str_detect(id, "_S$")  ~ "S",
  str_detect(id, "_As$") ~ "As",
  TRUE ~ NA_character_
)

pair_df = tibble(IGV = igv, tag = tag) %>%
  filter(!is.na(IGV), nzchar(IGV), !is.na(tag)) %>%
  distinct() %>%                          # don't double count duplicates
  count(IGV, tag, name="n") %>%
  pivot_wider(names_from = tag, values_from = n, values_fill = 0) %>%
  mutate(ok = (S > 0 & As > 0))

# failing IGVs
bad = pair_df %>% filter(!ok)

if (nrow(bad) > 0) {
  # write a report / stop
  # write.table(bad, opt$out_txt_file, sep="\t", row.names=FALSE, quote=FALSE)
  stop(paste("Missing sense/antisense pairs for ", nrow(bad), " IGV groups."), call.=FALSE)
}

write.table(df, outFile, sep = "\t", quote = F, row.names = F)







