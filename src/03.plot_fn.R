# This script aims to produce fn plots after ShortStack
# This is script must be polished since I'm just looking for 
# reads in a range of 18-27
# This script aims to produce fn plots after ShortStack
# This is script must be polished since I'm just looking for 
# reads in a range of 18-27
pacman::p_load(ggplot2, dplyr, stringr, reshape2, tidyr, tibble, purrr)

min_length = 18
max_length = 27

files = "first_nt_rlength.Rds"
tabs = readRDS(files)
#tabs = tabs[grepl("Ip", names(tabs))]

names(tabs) = sub(".trim.mapped.Rds", "", names(tabs)) 

tab_counts = imap_dfr(tabs, ~ as.data.frame(.x) %>%
                          rownames_to_column("length") %>%
                          pivot_longer(-length, names_to = "nt", values_to = "count") %>%
                          mutate(
                            sample = .y,
                            length = as.numeric(length)))

tabs  = split(tab_counts, tab_counts$sample)
# Give format and estimate CPM per library


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

# Format the final table
df$nt = gsub("T", "U", df$nt)
df$nt = factor(df$nt, levels = c("G", "U", "A", "C"))
df$length = as.character(df$length)  
df$length = factor(df$length, levels = c(min_length:max_length))
df = df %>% filter(!str_detect(tolower(ID), "input"))

# Colors per first nt
colors = c("#FFD600", "#DE0B00", "#05A61C", "#0314A8")
names(colors) = c("G", "U", "A", "C") # This line is not relevant but I'll keep it

tick_lengths = as.character(seq(min_length, max_length, 2))


# Fn plots
# By percentage
percetage_plot = df %>% ggplot(aes(x=length, y=percentage, fill=nt)) +
  geom_bar(stat="identity") + facet_wrap(~ID) +
  scale_fill_manual(values = colors, name= "5' nt") +
  theme_test() + ylab("% of reads") + xlab("Length") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18), 
        legend.position = "bottom",
        legend.text = element_text(size=14),
        strip.text = element_text(size=18)) +
  scale_x_discrete(breaks = tick_lengths) + ylim(c(0,35))

ggsave("length_dit_fn_percentage.pdf", device = "pdf", width = 5, height =  5,
       path = "./", plot = percetage_plot, dpi = 300)
