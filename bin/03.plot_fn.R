pacman::p_load(ggplot2, dplyr, stringr, reshape2, tidyr, tibble, purrr)
args = commandArgs(trailingOnly = TRUE)
min_length = as.integer(args[1])
max_length = as.integer(args[2])

files = list.files(pattern = "expanded.firstnt.tsv")
tabs = lapply(files, function(x){read.table(x, sep = "\t", header = T)})
#tabs = tabs[grepl("Ip", names(tabs))]
names(tabs) = sub(".ps.fastq.collapsed.expanded.firstnt.tsv", "", files) 

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
panel_width  <- 2.5 + len_span * 0.06   
panel_height <- 4 

# Figure size
plot_width  = panel_width  * ncol_facets
plot_height = panel_height * nrow_facets


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
    strip.text = element_text(size = 18), 
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