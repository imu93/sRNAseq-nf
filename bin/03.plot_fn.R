pacman::p_load(ggplot2, dplyr, stringr, reshape2, tidyr, tibble, purrr, optparse)

# Variables
option_list = list(
  make_option(c("--min"), type="numeric", default=18,
              help= "Min read length", metavar="character"),
  make_option(c("--max"), type="numeric", default=27, 
              help= "Max read length", metavar="character"),
  make_option(c("--reps"), type="logical", default = FALSE,
              help = "Plots per replicate",
              metavar = "logical"),
  make_option(c("--all_lengths"), type="logical", default = FALSE,
              help = "Analyze full observed length range (ignore --min/--max)",
              metavar = "logical")
); 

# Parse args
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
all_lengths = opt$all_lengths
min_length = opt$min
max_length = opt$max
reps = opt$reps

files = list.files(pattern = "expanded.firstnt.tsv")
tabs = lapply(files, function(x){read.table(x, sep = "\t", header = T)})
#tabs = tabs[grepl("Ip", names(tabs))]
names(tabs) = sub(".expanded.firstnt.tsv", "", files) 

tabs = lapply(tabs, function(m) {
  m[is.na(m)] <- 0
  m$sample = NULL
  m$total = NULL
  rownames(m) = m$length
  m$length = NULL
  return(m)
})

tab_counts = imap_dfr(tabs, ~ as.data.frame(.x) %>%
                        rownames_to_column("length") %>%
                        pivot_longer(-length, names_to = "nt", values_to = "count") %>%
                        mutate(
                          sample = .y,
                          length = as.numeric(length)))

tabs  = split(tab_counts, tab_counts$sample)

# Estimate percentage
f_lst = list()
for (i in names(tabs)) {
  tmp.tb = tabs[[i]]
  tmp.tb$percentage = 100*(tmp.tb$count/sum(tmp.tb$count))
  tmp.tb$ID = tmp.tb$sample %>%  sub("_[^_]*$","",.)
  f_lst[[i]] = tmp.tb
}

# Remove useless variables
tabs = NULL

# Bind gorups of libraries based on their condition (ID)
all_df = bind_rows(f_lst, .id = "sample")

if (all_lengths) {
  min_length <- min(all_df$length, na.rm = TRUE)
  max_length <- max(all_df$length, na.rm = TRUE)
}

# Build final tibble using average values per condition 
df = all_df %>%
  group_by(ID, length, nt) %>%
  dplyr::summarise(
    count = mean(count),
    percentage = mean(percentage),
    .groups = "drop"
  )

# keep requested length range first, then format
df = df %>% filter(length >= min_length, length <= max_length)

# Format the final table
df$nt = gsub("T", "U", df$nt)
df$nt = factor(df$nt, levels = c("G", "U", "A", "C"))

# drop inputs (comment out next line if it empties the data)
#df = df %>% filter(!str_detect(tolower(ID), "input"))

# if empty after filtering, exit cleanly (avoids facet error)
if (nrow(df) == 0) {
  message(sprintf("No rows left after filtering (length %d-%d and removing inputs). Skipping plot.", min_length, max_length))
  quit(save = "no", status = 0)
}

# now factor length with desired span
df$length = as.character(df$length)  
df$length = factor(df$length, levels = as.character(min_length:max_length))

# Colors per first nt
colors = c("#FFD600", "#DE0B00", "#05A61C", "#0314A8")
names(colors) = c("G", "U", "A", "C") # This line is not relevant but I'll keep it

tick_lengths = as.character(seq(min_length, max_length, 2))
tick_lengths = intersect(tick_lengths, levels(df$length))

# number of conditions
n_panels = length(unique(df$ID))

# choose ncol so that the facet grid is roughly square
ncol_facets = ceiling(sqrt(n_panels))
nrow_facets = ceiling(n_panels / ncol_facets)

# length span for x axis
len_span = max_length - min_length + 1

# per-panel dimensions (inches)
# base width grows with length span so long x-axes get more room
panel_width  = 2.5 + len_span * 0.06   
panel_height = 4 

# Figure size
plot_width  = panel_width  * ncol_facets
plot_height = panel_height * nrow_facets

# Facet strip size
strip_size = 12
strip_size = strip_size - 2.5 * log10(n_panels)
strip_size = strip_size - 1.0 * (ncol_facets - 2)
strip_size = strip_size - 0.03 * (len_span - 10)
strip_size = max(9, min(18, strip_size))

# Fn plots
# By percentage
percetage_plot =df %>% 
  ggplot(aes(x = length, y = percentage, fill = nt)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ ID, ncol = ncol_facets) +
  scale_fill_manual(values = colors, name = "5' nt") +
  theme_test() + 
  ylab("% of reads") + 
  xlab("Length") +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 16), 
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    strip.text = element_text(size = strip_size), 
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches")
  ) +
  scale_x_discrete(breaks = tick_lengths)

ggsave(
  "length_dit_fn_percentage.png",
  device = "png",
  path = "./",
  plot = percetage_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

################################################################################
# Same but for per-replicate
if (reps == TRUE) {
  
  df_rep = tab_counts %>%
    mutate(ID = sub("_[^_]*$","", sample)) %>%  
    group_by(sample, length, nt) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()
  
  df_rep = df_rep %>% filter(length >= min_length, length <= max_length)
  
  if (nrow(df_rep) == 0) {
    message(sprintf("No rows left for replicate plot (length %d-%d). Skipping.", min_length, max_length))
  } else {
    
    df_rep$nt = gsub("T", "U", df_rep$nt)
    df_rep$nt = factor(df_rep$nt, levels = c("G", "U", "A", "C"))
    
    df_rep$length = as.character(df_rep$length)
    df_rep$length = factor(df_rep$length, levels = as.character(min_length:max_length))
    
    tick_lengths = as.character(seq(min_length, max_length, 2))
    tick_lengths = intersect(tick_lengths, levels(df_rep$length))
    
    n_panels = length(unique(df_rep$sample))
    ncol_facets = ceiling(sqrt(n_panels))
    nrow_facets = ceiling(n_panels / ncol_facets)
    
    len_span = max_length - min_length + 1
    panel_width  = 2.5 + len_span * 0.06
    panel_height = 4
    plot_width  = panel_width  * ncol_facets
    plot_height = panel_height * nrow_facets
    
    strip_size = 16
    strip_size = strip_size - 2.5 * log10(n_panels)
    strip_size = strip_size - 1.0 * (ncol_facets - 2)
    strip_size = strip_size - 0.03 * (len_span - 10)
    strip_size = max(9, min(18, strip_size))
    
    
    
    perc_plot_rep = df_rep %>%
      ggplot(aes(x = length, y = percentage, fill = nt)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ sample, ncol = ncol_facets) +
      scale_fill_manual(values = colors, name = "5' nt") +
      theme_test() +
      ylab("% of reads") +
      xlab("Length") +
      theme(
        axis.text  = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        strip.text = element_text(size = strip_size),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches")
      ) +
      scale_x_discrete(breaks = tick_lengths)
    
    ggsave(
      "length_dit_fn_percentage.reps.png",
      device = "png",
      path = "./",
      plot = perc_plot_rep,
      width = plot_width,
      height = plot_height,
      dpi = 300
    )
  }
  
}
