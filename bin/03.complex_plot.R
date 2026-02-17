pacman::p_load(
  ggplot2, ggpubr, pals, reshape2, scales, dplyr, purrr, edgeR, tibble, stringr, grid
)

args = commandArgs(trailingOnly = T)
c_file = args[1]


fn_files  = list.files(pattern = ".expanded.firstnt.tsv")
cls_files = list.files(pattern = ".cls_mtx.tsv")
constrast = read.delim(c_file, sep = "\t", header = T)


constrast = read.delim(c_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

required = c("c_name","group")
if (!all(required %in% colnames(constrast))) {
  stop(
    "Contrast file missing required columns: ",
    paste(setdiff(required, colnames(constrast)), collapse=", "),
    "\nColumns found: ", paste(colnames(constrast), collapse=", "),
    "\nFile used as contrast: ", c_file
  )
}

fn = lapply(fn_files, function(x){
  tmp.x = read.table(x, header = TRUE)
  rownames(tmp.x) = tmp.x$length
  tmp.x$length = NULL
  tmp.x$sample = NULL
  tmp.x$total  = NULL
  tmp.x
})

cls = lapply(cls_files, function(x){
  tmp.x = read.table(x, header = TRUE)
  rownames(tmp.x) = tmp.x$length
  tmp.x$length = NULL
  tmp.x
})

names(fn)  = sub("\\..*", "", fn_files)
names(cls) = sub("\\..*", "", cls_files)

#groups = names(fn) %>%
#  strsplit("_") %>%
#  map(function(x) paste(x[c(1,2,3)], collapse = "_")) %>%
#  unlist() %>%
#  unique()

constrast$c_name = trimws(constrast$c_name)
constrast$group  = trimws(constrast$group)
groups = unique(constrast$group)
names(groups) = groups


pad_matrix_rows = function(mat, len_seq_chr) {
  rn = rownames(mat)
  missing = setdiff(len_seq_chr, rn)
  if (length(missing) > 0) {
    add = matrix(0, nrow = length(missing), ncol = ncol(mat),
                  dimnames = list(missing, colnames(mat)))
    mat = rbind(mat, add)
  }
  mat[len_seq_chr, , drop = FALSE]
}

pick_x_breaks = function(len_seq_num, max_labels = 8, prefer_even = TRUE) {
  len_seq_num = sort(unique(as.numeric(len_seq_num)))
  if (length(len_seq_num) == 0) return(NULL)
  
  # if few values, show all
  if (length(len_seq_num) <= max_labels) return(len_seq_num)
  
  rng  = range(len_seq_num)
  span = diff(rng)
  
  raw_step   = ceiling(span / (max_labels - 1))
  nice_steps = c(1, 2, 3, 4, 5, 6, 8, 10)
  step       = nice_steps[which.min(abs(nice_steps - raw_step))]
  
  start = rng[1]
  if (prefer_even && step %% 2 == 0) start = ceiling(start / 2) * 2
  
  brks = seq(start, rng[2], by = step)
  if (tail(brks, 1) != rng[2]) brks = c(brks, rng[2])
  brks
}


calc_font = function(n_panels) {
  base = if (n_panels <= 4) 16 else if (n_panels <= 9) 10 else if (n_panels <= 16) 7 else 16
  list(
    axis_title   = base,
    axis_text    = max(10, base - 2),
    strip        = base,
    legend_text  = max(10, base - 2),
    legend_title = base
  )
}

get_legend_grob = function(p) {
  get_legend(p + theme(legend.position = "bottom"))
}


n_panels = length(groups)
fs = calc_font(n_panels)

final_ncol = if (n_panels <= 2) 1 else if (n_panels <= 6) 2 else 3


prop_l_cpm      = list()
mean_l_cpm      = list()
mean_tab_cpm    = list()
group_lst_cpm   = list()
mean_fn_cpm_lis = list()

mean_pl = list()
prop_pl = list()

legend_nt  = NULL
legend_cls = NULL



for (i in names(groups)) {
  message(i)
  
  sample_names = constrast$c_name[constrast$group == i]
  libs = fn[sample_names]
  
  if (length(libs) == 0) {
    warning("No fn tables found for group: ", i)
    next
  }
  
  counts_lts = list()
  for (fn_lib in names(libs)) {
    tmp.lib = as.data.frame(libs[[fn_lib]])
    tmp.lib$len = rownames(tmp.lib)
    
    tmp.mlt = melt(tmp.lib)
    rownames(tmp.mlt) = paste0(tmp.mlt$variable, "_", tmp.mlt$len)
    
    tmp.mlt = select(tmp.mlt, value)
    colnames(tmp.mlt) = fn_lib
    
    counts_lts[[fn_lib]] = tmp.mlt
  }
  
  tmp.counts = do.call(cbind, counts_lts)
  tmp.counts[is.na(tmp.counts)] = 0
  mtx_cpm = cpm(tmp.counts)
  
  libs_nt = list()
  for (exp in colnames(mtx_cpm)) {
    x = mtx_cpm[, exp]
    
    cpm_mtx = as.data.frame(x) %>%
      rownames_to_column(var = "rownames") %>%
      mutate(first_char = substr(rownames, 1, 1)) %>%
      group_split(first_char) %>%
      setNames(c("df_A", "df_C", "df_G", "df_T")) %>%
      lapply(function(z) select(z, c(1,2))) %>%
      do.call(cbind, .)
    
    rownames(cpm_mtx) = sub(".*_", "", cpm_mtx$df_A.rownames)
    cpm_mtx = cpm_mtx[, !grepl("rownames", colnames(cpm_mtx)), drop = FALSE]
    colnames(cpm_mtx) = gsub("^df_", "", substr(colnames(cpm_mtx), 4, 4))
    
    libs_nt[[exp]] = cpm_mtx
  }
  
  vals = lapply(libs_nt, rownames) %>% unlist() %>% unique()
  len_all = sort(as.numeric(vals))
  len_seq = seq(min(len_all), max(len_all), by = 1)
  len_seq_chr = as.character(len_seq)
  
  libs_nt = lapply(libs_nt, function(m) pad_matrix_rows(m, len_seq_chr))
  
  mean_fn = Reduce("+", libs_nt) / length(libs_nt)
  
  colnames(mean_fn) = sub("T", "U", colnames(mean_fn))
  rownames(mean_fn) = as.character(as.numeric(rownames(mean_fn)))
  mean_fn[is.na(mean_fn)] = 0
  mean_fn_cpm_lis[[i]] = mean_fn
  
  mean_fn_long = as.data.frame(mean_fn)
  mean_fn_long$Length = as.numeric(rownames(mean_fn_long))
  mean_fn_long = melt(mean_fn_long, id.vars = "Length")
  colnames(mean_fn_long) = c("Length", "First_nt", "Reads")
  
  prop_fn = mean_fn / rowSums(mean_fn) * 100
  prop_fn[is.na(prop_fn)] = 0
  
  prop_fn_long = as.data.frame(prop_fn)
  prop_fn_long$Length = as.numeric(rownames(prop_fn_long))
  prop_fn_long = melt(prop_fn_long, id.vars = "Length")
  colnames(prop_fn_long) = c("Length", "First_nt", "Reads")
  
  mean_fn_long$First_nt = factor(mean_fn_long$First_nt, levels = rev(c("C","A","U","G")))
  prop_fn_long$First_nt = factor(prop_fn_long$First_nt, levels = rev(c("C","A","U","G")))
  
  mean_fn_long$IP = groups[[i]]
  prop_fn_long$IP = groups[[i]]
  
  x_breaks = pick_x_breaks(len_seq)
  
  colors_nt = rev(c("#00008AFF", "#09891A", "#E6110E", "#FFE400"))
  
  mean_p = ggplot(mean_fn_long, aes(x = Length, y = Reads, fill = First_nt)) +
    geom_col() +
    scale_fill_manual(values = colors_nt, name = "5' nt") +  # keep legend definition
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(labels = label_scientific(digits = 2)) +
    ylab("CPM") +
    facet_wrap(~IP) +
    theme_test() +
    theme(
      strip.text   = element_text(size = fs$strip, colour = "black"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = fs$axis_title),
      axis.text.y  = element_text(size = fs$axis_text, angle = 90, hjust = 0.5),
      legend.position = "none"
    )
  
  prop_p = ggplot(prop_fn_long, aes(x = Length, y = Reads, fill = First_nt)) +
    geom_col() +
    scale_fill_manual(values = colors_nt, name = "5' nt") +  # keep legend definition
    scale_x_continuous(breaks = x_breaks) +
    ylab("% of CPM") +
    xlab("Length distribution (nt)") +
    theme_test() +
    theme(
      axis.title.x = element_text(size = fs$axis_title),
      axis.text.x  = element_text(size = fs$axis_text),
      axis.title.y = element_text(size = fs$axis_title),
      axis.text.y  = element_text(size = fs$axis_text, angle = 90, hjust = 0.5),
      legend.position = "none"
    )
  
  if (is.null(legend_nt)) {
    legend_nt = get_legend_grob(
      ggplot(mean_fn_long, aes(x = Length, y = Reads, fill = First_nt)) +
        geom_col() +
        scale_fill_manual(values = colors_nt, name = "5' nt") +
        theme_test() +
        theme(
          legend.text  = element_text(size = fs$legend_text),
          legend.title = element_text(size = fs$legend_title)
        )
    )
  }
  
  mean_l_cpm[[groups[i]]]   = mean_p
  mean_tab_cpm[[groups[i]]] = mean_fn_long
  prop_l_cpm[[groups[i]]]   = prop_p
  
  group_lst_cpm[[groups[i]]] = ggarrange(
    mean_p,
    ggplot() + theme_minimal(),
    prop_p,
    ncol = 1, nrow = 3,
    heights = c(2, 0.05, 1),
    align = c("hv")
  )
}

for (i in names(groups)) {
  message(i)
  
  sample_names = constrast$c_name[constrast$group == i]
  libs = cls[sample_names]
  
  if (length(libs) == 0) {
    warning("No cls tables found for group: ", i)
    next
  }
  
  mean_fn = Reduce("+", libs) / length(libs)
  mean_fn = as.data.frame(mean_fn)
  
  colnames(mean_fn) = sub("^gene", "Protein-coding", colnames(mean_fn))
  colnames(mean_fn) = sub("DNA", "Transposon", colnames(mean_fn))
  colnames(mean_fn) = sub("LTR|LINE", "Retrotransposon", colnames(mean_fn))
  
  len_cls = sort(as.numeric(rownames(mean_fn)))
  len_seq = seq(min(len_cls), max(len_cls), by = 1)
  x_breaks = pick_x_breaks(len_seq)
  
  mean_fn_long = as.data.frame(mean_fn)
  mean_fn_long$Length = as.numeric(rownames(mean_fn_long))
  mean_fn_long = melt(mean_fn_long, id.vars = "Length")
  colnames(mean_fn_long) = c("Length", "Class", "Reads")
  
  prop_fn = mean_fn / rowSums(mean_fn) * 100
  prop_fn[is.na(prop_fn)] = 0
  prop_fn_long = as.data.frame(prop_fn)
  prop_fn_long$Length = as.numeric(rownames(prop_fn_long))
  prop_fn_long = melt(prop_fn_long, id.vars = "Length")
  colnames(prop_fn_long) = c("Length", "Class", "Reads")
  
  mean_fn_long$Class = factor(
    mean_fn_long$Class,
    levels = c(
      "Unannotated","Other_repeat","Other","Other_ncRNA","snRNA","snoRNA",
      "piRNA","lincRNA","rRNA","yRNA","tRNA","miRNA","Unknown",
      "Retrotransposon","pseudogene","Protein-coding"
    )
  )
  prop_fn_long$Class = factor(
    prop_fn_long$Class,
    levels = c(
      "Unannotated","Other","Other_repeat","Other_ncRNA","snRNA","snoRNA",
      "piRNA","lincRNA","rRNA","yRNA","tRNA","miRNA","Unknown",
      "Retrotransposon","Transposon","pseudogene","Protein-coding"
    )
  )
  
  colors_cls = rev(c(
    "#0051FFFF","#0092EFFF","#5FD9BB","#13F240FF","#9AFF16FF","#F5E00A","#FFBA00FF","#FF7F00FF",
    "#FF3E00","#D70B3D","#F55E87","#601573","#3D0B31","#2F2F30","#76767A","#DBDBDB"
  ))
  
  exp = groups[[i]]
  mean_fn_long$exp = exp
  prop_fn_long$exp = exp
  
  mean_p = ggplot(mean_fn_long, aes(x = Length, y = Reads, fill = Class)) +
    geom_col() +
    scale_fill_manual(values = colors_cls, guide = guide_legend(reverse = TRUE), name = "Class") +
    scale_x_continuous(breaks = x_breaks) +
    ylab("Number of aligned reads") +
    facet_wrap(~exp) +
    theme_test() +
    theme(
      strip.text = element_text(size = fs$strip, colour = "black"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = fs$axis_title),
      axis.text.y  = element_text(size = fs$axis_text, angle = 90, hjust = 0.5),
      legend.position = "none"
    )
  
  prop_p = ggplot(prop_fn_long, aes(x = Length, y = Reads, fill = Class)) +
    geom_col() +
    scale_fill_manual(values = colors_cls, guide = guide_legend(reverse = FALSE), name = "Class") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = c(0, 50, 100)) +
    ylab("% of reads") + xlab("sRNA length (nt)") +
    theme_test() +
    theme(
      axis.title.x = element_text(size = fs$axis_title),
      axis.text.x  = element_text(size = fs$axis_text),
      axis.title.y = element_text(size = fs$axis_title),
      axis.text.y  = element_text(size = fs$axis_text, angle = 90, hjust = 0.5),
      legend.position = "none",
      legend.key.size = unit(.9, "cm")
    )
  
 
  if (is.null(legend_cls)) {
    legend_cls = get_legend_grob(
      ggplot(prop_fn_long, aes(x = Length, y = Reads, fill = Class)) +
        geom_col() +
        scale_fill_manual(values = colors_cls, name = "Class") +
        theme_test() +
        theme(
          legend.text  = element_text(size = fs$legend_text),
          legend.title = element_text(size = fs$legend_title)
        )
    )
  }
  
  mean_pl[[exp]] = mean_p
  prop_pl[[exp]] = prop_p
}



choose_ncol = function(n) {
  if (n <= 1) 1
  else if (n <= 4) 2
  else if (n <= 9) 3
  else 4
}
mixed_lst = list()

for (i in names(groups)) {
  g = groups[[i]]
  
  mixed_lst[[g]] = ggarrange(
    mean_l_cpm[[g]],             
    ggplot() + theme_void(),      
    prop_pl[[g]],                 
    ncol = 1, nrow = 3,
    heights = c(1.2, 0.08, .7)
  )
}

n_groups = length(groups)
mixed_ncol = choose_ncol(n_groups)

main_grid = ggarrange(
  plotlist = mixed_lst,
  ncol = mixed_ncol
)

final_plot = ggarrange(
  legend_nt,
  main_grid,
  legend_cls,
  ncol = 1,
  heights = c(0.10, 1, 0.18)
)

ggsave("fn_cls_cpm.png", width = 12)
