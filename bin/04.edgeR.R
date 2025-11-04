pacman::p_load(edgeR, ggplot2, ggpubr, dplyr, purrr, gplots, optparse, 
               rtracklayer, plyranges, reshape, kableExtra, 
               ggbreak, stringr, pals, stringi)

# Variables
option_list = list(
  make_option(c("-i", "--counts_file"), type="character", default=NULL,
              help= "Path to featureCounts counts table", metavar="character"),
  make_option(c("-c", "--contrast_file"), type="character", default=1, 
              help= "Single row txt file specifying the contrast", metavar="character"),
  make_option(c("-t", "--threshold_inc"), type="logical", default = FALSE,
              help = "For complex contrast files FC and FDR can be incleded as extra columns",
              metavar = "logical"),
  make_option(c("-l", "--logfold_change"), type="numeric", default=1, 
              help="log2 fold change for DEA", metavar="numeric"),
  make_option(c("-f", "--fdr"), type="numeric", default=0.05, 
              help= "FDR cutoff value for DEA", metavar="numeric"),
  make_option(c("-x", "--hk_norm"), type="logical", default=FALSE, 
              help= "Housekeeping like normalization", metavar="logical"),
  make_option(c("-n", "--norm_feature"), type="character", default="rRNA_S", 
              help= "Genomic feature to estimate normFactors", metavar="character"),
  make_option(c("-s", "--stringent_tmm"), type="logical", default=FALSE, 
              help= "Stringency level for feature specific TMM", metavar="logical"),
  make_option(c("-d", "--disable_tmm_for_fbn"), type = "logical", default = FALSE,
              help = "Set true in case you whant to diable TMM for FB normFactor estimation")
); 

# Parse args
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Perhaps in a for 
if (is.null(opt$logfold_change)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for logfold_change", call.=FALSE)
}

if (is.null(opt$fdr)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for fdr", call.=FALSE)
}

tresh_c = as.logical(opt$threshold_inc)
if (tresh_c == TRUE) {
  print("-l and -f will be omited and the contrast file will be used")
}


# In Files
FC=opt$logfold_change
FDR=opt$fdr
normFeature = opt$norm_feature
# If user supplied alternation with |, make sure it's grouped & anchored once
if (grepl("\\|", normFeature)) {
  # remove leading ^ in each part then group & anchor once
  parts = unlist(strsplit(normFeature, "\\|"))
  parts = sub("^\\^", "", parts)
  normFeature = paste0("^(", paste(parts, collapse="|"), ")")
} else {
  # if single feature is passed then ensure it is anchored at start
  if (!grepl("^\\^", normFeature)) {
    normFeature = paste0("^", normFeature)
  }
}

strict = as.logical(opt$stringent_tmm)
all_cat_norm =  as.logical(opt$disable_tmm_for_fbn)
cot_file = opt$contrast_file
hk_norm = as.logical(opt$hk_norm)

  
##################################################################################
# Inpath and featureCounts table
if (is.null(opt$counts_file)) stop("--counts_file is required", call.=FALSE)
countsFile = opt$counts_file

# Define out paths
outPath_1 = "tmm/"
outPath_2 = "norm/"
outPath_3 = "dge/"
outPath_6 = "figures/"

# Build directories
lst_dirs = list(outPath_1, outPath_2, outPath_3, outPath_6)
sapply(lst_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)


# Read table of contrast
conds = read.table(cot_file, header = T)
c_lst = split(conds, conds$contrast)

verify_cotrast = function(x) {
  tmp.x = x
  # First I need to check that each condition has at least 2 replicates
  if (!all((split(tmp.x, tmp.x$group) %>% sapply(nrow)) >= 2)) {
    stop("Need ≥2 replicates per group.", call. = FALSE)
  } else {
    print("N replicates: Pass")
  }
  # Now the order
  if (length(split(tmp.x, tmp.x$order)) != 2) {
    stop("Wrong order format, the user must provide at least two conditions", call. = FALSE)
  } else {
    print("Order format: Pass")
  }
  
  req = c("contrast", "type", "group", "c_name", "order")
  miss = setdiff(req, colnames(tmp.x))
  if (tresh_c == F && length(miss) > 0) {
    stop(sprintf("When -t/--treshold_inc is FALSE, contrast file must have columns: 5. Missing: %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  req = c("contrast", "type", "group", "c_name", "order", "FC", "FDR")
  miss = setdiff(req, colnames(tmp.x))
  if (tresh_c == T && length(miss) > 0) {
    stop(sprintf("When -t/--treshold_inc is TRUE, contrast file must have columns: FC, FDR. Missing: %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  if (!is.numeric(tmp.x$order) || !all(tmp.x$order %in% c(1,2))) {
    stop("'order' must be numeric and only 1 or 2.", call. = FALSE)
  }
  
  if (tresh_c == T) {
    if (!is.numeric(tmp.x$FC)) {
      stop("FC colummn must be numeric", call. = F)
    }
    if (!is.numeric(tmp.x$FDR)) {
      stop("FDR colummn must be numeric", call. = F)
    }
    # FDR range should match 
    if (any(tmp.x$FDR < 0 | tmp.x$FDR > 1))
      stop("FDR must be between 0 and 1.", call. = FALSE)
    # Each contrast ID must  be constant
    if (length(unique(tmp.x$FC)) != 1)
      stop(sprintf("FC must be constant within each contrast. Offenders: %s",
                   unique(tmp.x$contrast)), call. = FALSE)
    if (length(unique(tmp.x$FDR)) != 1)
      stop(sprintf("FDR must be constant within each contrast. Offenders: %s",
                   unique(tmp.x$contrast)), call. = FALSE)
  }
}


# Run 
for (i in names(c_lst)) {
  print(paste0("Verifing:", i))
  verify_cotrast(c_lst[[i]])
}


contrast_vec = unique(conds$contrast)
conts = setNames(as.list(contrast_vec), contrast_vec)

if (tresh_c == TRUE) {
  l_lfc = split(conds, conds$contrast) %>% sapply(function(x){unique(x$FC)})
  l_fdr = split(conds, conds$contrast) %>% sapply(function(x){unique(x$FDR)})
  names(l_lfc) = names(conts)
  names(l_fdr) = names(conts)
} else{
  l_lfc = rep(FC, length(conts))
  l_fdr = rep(FDR, length(conts))
  names(l_lfc) = names(conts)
  names(l_fdr) = names(conts)
}

##################################################################################
# Lists to store the results
dt_lst_nonorm = list()
dt_lst_normrrna = list()
mds_list = list()
norm_norm_simp_list = list()
mir_ecol_list = list()
tmm_ma = list()
gene_norm_list = list()
ma_norm_lst = list()



# This function will help to have more control on the contrast based on the order
# column and the c_name of the contrast table
# Also it will add support for contrast with "-" in the IDs 
order_contrast_from_cname = function(design, tmp.tb) {
  trim_last_token = function(x) sub("_[^_]+$", "", x)
  grpA = unique(trim_last_token(tmp.tb$c_name[tmp.tb$order == 1]))
  grpB = unique(trim_last_token(tmp.tb$c_name[tmp.tb$order == 2]))
  if (grepl("-", grpA)) {
    grpA =  str_replace(grpA,"-", ".")
  }
  if (grepl("-", grpB)) {
    grpB =  str_replace(grpB,"-", ".")
  }
  comp = setNames(numeric(ncol(design)), colnames(design))
  # Weight each side equally if multiple groups appear on a side
  if (length(grpA)) comp[intersect(grpA, names(comp))] <-  1 / max(1, length(grpA))
  if (length(grpB)) comp[intersect(grpB, names(comp))] <- -1 / max(1, length(grpB))
  if (all(comp == 0)) {
    stop("No matching design columns for groups derived from c_name. ",
         "Check that design colnames == sub('_[^_]+$', '', c_name).",
         call. = FALSE)
  }
  comp
}

for ( cont in names(conts)) {
  print(cont)
  contrast_id = cont
  # Get the condition IDs for plots
  tmp.tb = conds[conds$contrast == cont,]
  cond_A = tmp.tb[tmp.tb$order == 1,]$group %>% unique()
  cond_B = tmp.tb[tmp.tb$order == 2,]$group %>% unique()
  # I will use the type column to get the correct libraries
  libs = conds[conds$contrast == cont,]$type %>% unique()
  if (length(libs) > 1) {
    libs = paste0(libs, collapse = "|")
  }else{
    libs = libs
  }
  # Read table and create groups 
  counts = read.delim(countsFile, check.names = F)
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.collapsed.*", "", colnames(counts))
  # As for the order_contrast_from_cname, some library names can have a "-"
  # in the ID, I will add support by replacing the "-" by "."
  if (grepl("-", colnames(counts)) %>% sum() > 0) {
    colnames(counts) = str_replace( colnames(counts) ,"-", ".")
  }
  groups = colnames(counts) %>%  sub("_[^_]*$", "", .) %>% factor()
  print(groups)
  ##################################################################################
  lib.size = colSums(counts)
  # How many reads do we have?
  round((lib.size)/1e6, 2)
  # Create DGE and filter by expression
  dge = DGEList(counts=counts, group=groups, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 10)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  f.exp.counts = colSums(dge$counts)
  dge = calcNormFactors(dge, method = "TMM")
  ##################################################################################
  colors = rich.colors(2) %>% rep(2) %>% rev() %>% adjustcolor(alpha.f = .2)
  shapes = ifelse(grepl("input", levels(dge$samples$group)), 15, 16)
  names(colors) = levels(groups)
  names(shapes) = levels(groups)
  
  lab = function(x){
    x %>%
      str_replace_all("_", " ") %>%
      str_squish() %>%
      str_to_title()
  }
  
  groups = dge$samples$group            
  lvl    = levels(groups)
  labs_auto = setNames(lab(lvl), lvl)
  pal = setNames(adjustcolor(rich.colors(length(lvl)), alpha.f = 0.7), lvl)
  
  mds = plotMDS(dge, col = pal[groups], pch = 16, cex = 4)
  
  df_mds = data.frame(
    dim1 = mds$x,
    dim2 = mds$y,
    group = groups,
    exp   = cont
  )
  
  mds_gp = ggplot(df_mds, aes(dim1, dim2, color = group)) +
    geom_point(size = 12) +
    theme_light() +
    scale_color_manual(values = pal, breaks = lvl, labels = labs_auto, name = "") +
    xlab(sprintf("Leading logFC dim 1 %d%%", round(100 * mds$var.explained[1]))) +
    ylab(sprintf("Leading logFC dim 2 %d%%", round(100 * mds$var.explained[2]))) +
    facet_wrap(~exp) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.position="bottom",
          legend.text=element_text(size=16),
          strip.text=element_text(size=18, colour="black"),
          plot.margin=unit(c(.5,.5,.5,.5),"inches")) +
    guides(color = guide_legend(override.aes = list(size = 14)))
    
  mds_list[[cont]] = mds_gp
  ##################################################################################
  design = model.matrix(~0+dge$samples$group)
  colnames(design) =  levels(dge$samples$group)
  
  ##################################################################################
  dge = estimateDisp(dge, design=design, robust=TRUE)
  #plotBCV(dge)
  # As expected, we have variance in low expressed regions and low in highly expressed 
  # Let's fit glm
  # First order the contrast based on the design and use the order col of contrast file 
  # to improve contrast specification
  ncont = order_contrast_from_cname(design, tmp.tb)
  ncont = rev(sort(ncont))
  # set the contrast
  ncont = paste0(names(ncont[1]),"-",names(ncont[2]))
  # And fit the model
  fit = glmFit(dge,design = design, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc)
  topTags(gt)
  dt = decideTests(gt, adjust.method="BH",  p.value=cont_fdr, lfc = cont_fc)
  tb_nonorm = table(dt)
  dt_lst_nonorm[[ncont]] = tb_nonorm
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  outToptable = paste0(outPath_1, "other.tmm_topTbl.txt")
  outToptable = sub("other", ncont, outToptable)
  write.table(topTable,outToptable, quote = F, sep = "\t")
  ##################################################################################
  de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 
                      'log-FC'= topTable$logFC, 'FDR'=topTable$FDR, 'cl'= dt)
  colnames(de_tbl) = c('Ave_CPM', 'logFC', 'FDR', 'cl')
  de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0",
                      ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
  
  de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1,
                      ifelse(de_tbl$col == "#333333", 2, 3))
  de_tbl = de_tbl[order(de_tbl$ord),]
  de_tbl$col = factor(de_tbl$col, levels = c("#5E2129", "#C0C0C0", "#333333"))
  de_tbl$exp = sub("-", " vs ", ncont)
  de_tbl$exp = factor(de_tbl$exp)
  
  outDEtbl = paste0(outPath_1, "other.tmm_deTbl.txt")
  outDEtbl = sub("other", ncont, outDEtbl)
  write.table(de_tbl,outDEtbl, quote = F, sep = "\t")
 
  labs_y = function(A, B) bquote(log[2] ~ ' FC ' ~ .(A) ~ '/' ~ .(B))  
 p1 = de_tbl %>% ggplot(aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() + theme_test() +
    xlab(bquote(log[2]~' CPM')) +
    ylab(labs_y(cond_A, cond_B)) +
    scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0", "#333333"), alpha.f = .6), 
                       labels=c(cond_A, 'not-DE', cond_B), name=NULL) +
    geom_hline(aes(yintercept = 0), colour = "red",
               linewidth = 1, linetype="dashed", alpha=0.4) +
    facet_wrap(~ exp) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20), 
          legend.position = "bottom",
          legend.text = element_text(size=16),
          strip.text = element_text(size=18, colour = "black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) +
   guides(color = guide_legend(override.aes = list(size=12)))
  
 ##################################################################################
  feats = c("miRNA_S", "rRNA_S", "tRNA_S", "piRNA_S")
  df_lst = list()
  for (i in feats) {
    df = de_tbl
    df$col = ifelse(grepl(i, rownames(df)), "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df = df[order(df$ord),]
    df$exp = i
    df_lst[[i]] = df
  }
  
  names(df_lst) = NULL
  df2plot = do.call(rbind, df_lst)
  
  p2 = ggplot(df2plot, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) +  xlab(bquote(log[2]~' CPM')) + theme_test() +
    ylab(labs_y(cond_A, cond_B)) +
    scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0"), alpha.f = .6)) +
    geom_hline(aes(yintercept = 0), colour = "red",
               linewidth = 1, linetype="dashed", alpha=0.4) +
    facet_wrap(~ exp) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20), 
          legend.position = "bottom",
          legend.text = element_text(size=16),
          strip.text = element_text(size=18, colour = "black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) +
    guides(color = guide_legend(override.aes = list(size=12))) 
  
  tmm_ma_p = ggarrange(p1, p2, ncol = 1, heights = c(1, .9))
  tmm_ma[[ncont]] = tmm_ma_p
  
  outDGE = paste0(outPath_3, "other.tmm_dge.Rds")
  outDGE = sub("other", ncont, outDGE)
  saveRDS(dge, outDGE)
  
  ##################################################################################
  if (hk_norm == FALSE) {
    print("Only TMM-based DEA was conducted")
    next
  }
  
  # It is important to realize that IP vs Input are not fair contrasts.
  # These libraries are inherently compositionally biased due to the nature of small RNA sequencing.
  # This bias arises because Argonaute IPs selectively enrich for sRNAs with specific properties
  # such as defined length distributions and 5′ nucleotide preferences.
  #
  # As a result, standard normalization methods that assume similar global distributions (like TMM)
  # may fail or require adjustment. To address this, we focus on a subset of regions (e.g., rRNA fragments)
  # that are not expected to be enriched in IPs and thus serve as a more stable reference for scaling.
  
  get_trimmed_category_counts = function(dge_full, category_pattern = "^rRNA_S", groups, 
                                         min_count = 10, logratioTrim = 0.3, sumTrim = 0.05,
                                         stringent = FALSE, use_all_category = FALSE) {
    
    # 1. Get the rows corresponding to the category of interest
    category_counts = dge_full$counts[grepl(category_pattern, rownames(dge_full$counts)), ]
    
    if (nrow(category_counts) == 0L) {
      stop("No rows match category_pattern: ", category_pattern, call. = TRUE)
    }
    
    # In case the you are not interested the TMM-like approach you can use the 
    # whole category
    if (isTRUE(use_all_category)) {
      # Optional guard: warn if any sample has 0 counts in this category
      tmp.counts = colSums(category_counts) == 0
      if (any(tmp.counts)) {
        warning("Some libraries have zero counts for the selected category (",
                paste(colnames(category_counts)[tmp.counts], collapse = ", "),
                ") Norm factors may be unstable.")
      }
      return(category_counts)
    }
    
    
    # 2. Get the corresponding N reads for this category
    lib_sizes_cat = colSums(category_counts)
    
    # 3. Build a new dge object but just for this category
    dge_cat = DGEList(counts = category_counts, group = groups, lib.size = lib_sizes_cat)
    
    # 4. Filter out lowly expressed regions. Since I'm using the libsize only for this category
    #    I do not expect to loose to many of this regions because of this filter. 
    keep = filterByExpr(dge_cat, min.count = min_count)
    dge_cat = dge_cat[keep, , keep.lib.sizes = FALSE]
    message("Features kept after filterByExpr:"); print(table(keep))
    
    # 5. Get the reference sample based on geometric mean
    lib_sizes = dge_cat$samples$lib.size
    geo_mean_lib_size = exp(mean(log(lib_sizes)))
    ref = which.min(abs(log(lib_sizes / geo_mean_lib_size)))
    
    # 6.Now I'll compare each sample against the reference. 
    test_samples = setdiff(seq_len(ncol(dge_cat$counts)), ref)
    
    trimming_masks = lapply(test_samples, function(test) {
      # I need to compare the logFC (M-value) between my reference and each test
      # I will also include a pseudocount to avoid inf with log2(0) 
      logR = log2((dge_cat$counts[, test] + 0.5) / lib_sizes[test]) -
        log2((dge_cat$counts[, ref] + 0.5) / lib_sizes[ref])
      
      # Now the average expression value between ref and each test (A-value)
      # Again, since this is log2 y need pseudocounts
      absE = 0.5 * (log2((dge_cat$counts[, test] + 0.5) / lib_sizes[test]) +
                      log2((dge_cat$counts[, ref] + 0.5) / lib_sizes[ref]))
      
      # This is the key to behave like tmm. I will remove (30%) extreme values in both up and down based on FC 
      keep_M = rank(logR) > length(logR) * logratioTrim & rank(logR) < length(logR) * (1 - logratioTrim)
      # Now the same but for expression where I'll just shrink by 5% of the lower and higer expressed regions
      keep_A = rank(absE) > length(absE) * sumTrim & rank(absE) < length(absE) * (1 - sumTrim)
      # Get the vector per lib
      keep = keep_M & keep_A
      return(keep)
    })
    # 7. Now only keep stable regions in at least one contrast vs ref.
    # Why?
    # If we required regions to pass the trimming in *all* comparisons (i.e., Reduce(&, ...)),
    # we would likely discard most or all features due to natural biological and technical variability
    # This is especially true for sRNA-seq data, where Argonaute IPs often have highly skewed composition
    #
    # By using Reduce(|, ...), we are instead selecting regions that pass trimming in *at least one* comparison
    # This approach is less strict, but more realistic, it allows us to retain regions that are generally stable,
    # even if they appear slightly extreme in some individual comparisons due to noise or enrichment.
    
    if (stringent == TRUE) {
      keep_any = Reduce(`&`, trimming_masks)
    } else{
      keep_any = Reduce(`|`, trimming_masks)
    }
    
    filtered_counts = dge_cat$counts[keep_any, , drop = FALSE]
    
    min_features = 5
    if (nrow(filtered_counts) <  min_features && isTRUE(stringent)) {
      stop("Insufcient feature instances after TMM. 
           Try a different feature set, relax trimming or test use_all_category = TRUE.")
    }
    
    return(filtered_counts)
  }
  
  filtered_counts = get_trimmed_category_counts(
    dge_full = dge,               
    category_pattern = normFeature, 
    groups = groups,
    stringent = strict,
    use_all_category = all_cat_norm
  )
  
  normCounts = filtered_counts
  sumNormRNAs = colSums(normCounts)
  
  ##################################################################################
  counts = read.delim(countsFile, check.names = F)
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.collaps.*", "", colnames(counts))
  if (grepl("-", colnames(counts)) %>% sum() > 0) {
    colnames(counts) = str_replace( colnames(counts) ,"-", ".")
  }
  groups = colnames(counts) %>%  sub("_[^_]*$", "", .) %>% factor()
  lib.size = colSums(counts)
  # How many reads do we have?
  round((lib.size)/1e6, 2)
  # Create DGE and filter by expression
  dge = DGEList(counts=counts, group=groups, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 10)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  
  totalLibsize = colSums(dge$counts)
  normFactors = sumNormRNAs / totalLibsize
  normFactors = normFactors / (prod(normFactors)^(1/length(normFactors)))
  
  dge$samples$norm.factors = normFactors
  
  design = model.matrix(~0+dge$samples$group)
  colnames(design) =  levels(dge$samples$group)
  dge = estimateDisp(dge, design=design, robust=TRUE)
  
  #plotBCV(dge)
  ##################################################################################
  fit = glmFit(dge,design = design,dispersion = dge$common.dispersion)
  #plotQLDisp(fit)
  contrast = makeContrasts(ncont,levels = design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc)
  topTags(gt)
  dt = decideTests(gt, adjust.method="BH", p.value=cont_fdr, lfc = cont_fc)
  normrrna_dt = table(dt)
  dt_lst_normrrna[[ncont]] = normrrna_dt
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  
  outToptable = paste0(outPath_2, "other.norm_topTbl.txt")
  outToptable = sub("other", ncont, outToptable)
  write.table(topTable,outToptable, quote = F, sep = "\t")
  
  
  ##################################################################################
  de_tbl_norm = data.frame('Ave_CPM'=topTable$logCPM, 
                           'log-FC'= topTable$logFC, "FDR"=topTable$FDR,'cl'= dt)
  colnames(de_tbl_norm) = c('Ave_CPM', 'logFC', 'FDR', 'cl')
  de_tbl_norm$col = ifelse(de_tbl_norm$cl == 0, "#C0C0C0", 
                           ifelse(de_tbl_norm$cl == 1, "#15C1FB", "#C0C0C0"))
  
  de_tbl_norm$col = ifelse(grepl(normFeature, rownames(de_tbl_norm)),
                           "#DBBF10", de_tbl_norm$col)
  
  de_tbl_norm$col = factor(de_tbl_norm$col, levels = c("#15C1FB", "#C0C0C0",  "#DBBF10"))
  de_tbl_norm$ord = 1
  de_tbl_norm$ord = ifelse(de_tbl_norm$col == "#15C1FB", 2, de_tbl_norm$ord)
  de_tbl_norm$ord = ifelse(de_tbl_norm$col == "#DBBF10", 3, de_tbl_norm$ord)
  de_tbl_norm$ord = factor(de_tbl_norm$ord)
  
  de_tbl_norm = de_tbl_norm[order(de_tbl_norm$ord),]
  de_tbl_norm[de_tbl_norm$col == "#DBBF10",]
  
  de_tbl_norm$exp = ncont
  
  outDEtbl = paste0(outPath_2, "other.norm_deTbl.txt")
  outDEtbl = sub("other", ncont, outDEtbl)
  write.table(de_tbl_norm,outDEtbl, quote = F, sep = "\t")
  
  normFeaturel = ifelse(grepl("\\|", normFeature), "Norm_regions", normFeature)
  norm_label = sub("^\\^", "", normFeature)
  ma = ggplot(de_tbl_norm, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_test() +
    ylab(labs_y(cond_A, cond_B)) +
    scale_color_manual(values = adjustcolor(c("#155AA3", "#C0C0C0", "red"), alpha.f = .6),
                       labels=factor(c(cond_A, 'non-DE', norm_label), 
                                     levels = c("non-DE",cond_A, norm_label)), 
                       name=NULL) + 
    facet_wrap(~ exp) +
    geom_hline(aes(yintercept = 0), colour = "red", 
               linewidth = 1, linetype="dashed", alpha=0.7)  +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20), 
          legend.position = "bottom",
          legend.text = element_text(size=16),
          strip.text = element_text(size=18, colour = "black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches")) +
  guides(color = guide_legend(override.aes = list(size=12))) 
  ma_norm_lst[[cont]] = ma
  
  feats = c("miRNA_S", "rRNA_S", "tRNA_S", "piRNA_S")
  df_lst = list()
  for (i in feats) {
    df = de_tbl_norm
    df$col = ifelse(grepl(i, rownames(df)), "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df = df[order(df$ord),]
    df$exp = i
    df_lst[[i]] = df
  }
  
  names(df_lst) = NULL
  df2plot = do.call(rbind, df_lst)
  
  p2 = ggplot(df2plot, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) +  xlab(bquote(log[2]~' CPM')) + theme_test() +
    ylab(labs_y(cond_A, cond_B)) +
    theme(axis.title.x = element_text(size = 20),
          axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text( size = 20),
          axis.text.y = element_text( size = 20)) +
    geom_hline(aes(yintercept = 0), colour = "red",
               linewidth = 1, linetype="dashed", alpha=0.4) +
    facet_wrap(~ exp) +
    scale_color_manual(values = adjustcolor(c("#155AA3","#C0C0C0"), alpha.f = .6)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    theme( strip.text = element_text(size = 20, colour = "black")) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches"))
  
  norm_ma_p = ggarrange(ma, p2,  ncol = 1, heights = c(1, .9))
  norm_norm_simp_list[[ncont]] = norm_ma_p
  outDGE = paste0(outPath_3, "other.norm_dge.Rds")
  outDGE = sub("other", ncont, outDGE)
  saveRDS(dge, outDGE)

}




##################################################################################
for (i in names(tmm_ma)) {
  print(i)
  tmp.plot = tmm_ma[[i]]
  fname = paste0("MA.",i, "_tmm.png")
  ggsave(plot = tmp.plot, device = "png", width = 12, height = 14, dpi = 300,
         filename = fname, 
         path = outPath_6)
}


for (i in names(norm_norm_simp_list)) {
  print(i)
  tmp.plot = norm_norm_simp_list[[i]]
  fname = paste0("MA.",i, "_norm.png")
  ggsave(plot = tmp.plot, device = "png", width = 12, height = 14, dpi = 300,
         filename = fname, 
         path = outPath_6)
}

for (i in names(mds_list)) {
  print(i)
  tmp.plot = mds_list[[i]]
  fname = paste0("MDS.",i, ".png")
  ggsave(plot = tmp.plot, device = "png", width = 6, height = 5, dpi = 300,
         filename = fname, 
         path = outPath_6)
}
##################################################################################
