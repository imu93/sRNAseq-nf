# This script aims to produce fn plots after ShortStack
# This is script must be polished since I'm just looking for 
# reads in a range of 18-27
 
pacman::p_load(ggplot2, dplyr, stringr)
args = commandArgs(trailingOnly = T)
min_length = args[1]
max_length = args[2]

files = list.files(pattern = ".*length_firstnt.txt")
tabs = lapply(files, read.table)
names(tabs) = files

# Give format and estimate CPM per library
cpm_list = lapply(tabs, function(x){
  tmp.x = x
  colnames(tmp.x) = c("length", "fn")
  tmp.x = tmp.x[tmp.x$fn != "N",]
  total_reads = nrow(tmp.x)
  tmp.df = tmp.x %>%
    count(length, fn) %>%
    group_by(length) %>%
    mutate(CPM = (n / total_reads) * 1e6)
})

# Estimate percentage
cpm_f_lst = list()
for (i in names(cpm_list)) {
  tmp.tb = cpm_list[[i]]
  tmp.tb$percentage =  100*(tmp.tb$n/sum(tmp.tb$n))
  cpm_f_lst[[i]] = tmp.tb
}

# Remove useless variables
tabs = NULL
cpm_list = NULL  
tmp.tb = NULL

# Edit lib IDs to group 
names(cpm_f_lst) = names(cpm_f_lst) %>% 
  str_replace("\\.length_firstnt\\.txt$", "") %>% 
  str_replace("_[^_]+$", "")

# Bind gorups of libraries based on their condition (ID)
all_df = bind_rows(cpm_f_lst, .id = "Condition")

# Build final tibble using average values per condition 
df = all_df %>%
  group_by(Condition, length, fn) %>%
  summarise(
    n = mean(n),
    CPM = mean(CPM),
    percentage = mean(percentage),
    .groups = "drop"
  )

# Format the final table
df$fn = gsub("T", "U", df$fn)
df$fn = factor(df$fn, levels = c("G", "U", "A", "C"))
df$length = as.character(df$length)  
df$length = factor(df$length, levels = c(min_length:max_length))
df = df %>% filter(!str_detect(tolower(Condition), "input"))

# Colors per first nt
colors = c("#FFD600", "#DE0B00", "#05A61C", "#0314A8")
names(colors) = c("G", "U", "A", "C") # This line is not relevant but I'll keep it

tick_lengths = as.character(seq(min_length, max_length, 2))


# Fn plots
# By percentage
percetage_plot = df %>% ggplot(aes(x=length, y=percentage, fill=fn)) +
  geom_bar(stat="identity") + facet_wrap(~Condition) +
  scale_fill_manual(values = colors, name= "5' nt") +
  theme_test() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18), 
        legend.position = "bottom",
        legend.text = element_text(size=14),
        strip.text = element_text(size=18)) +
    scale_x_discrete(breaks = tick_lengths)

# And CPM just in case
cpm_plot = df %>% ggplot(aes(x=length, y=CPM, fill=fn)) +
  geom_bar(stat="identity") + facet_wrap(~Condition) +
  scale_fill_manual(values = colors, name= "5' nt") +
  theme_test() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18), 
        legend.position = "bottom",
        legend.text = element_text(size=14),
        strip.text = element_text(size=18)) +
  scale_x_discrete(breaks = tick_lengths)

# Save pdf files

ggsave("length_dit_fn_percentage.pdf", device = "pdf", width = 8, height =  10,
       path = "./", plot = percetage_plot, dpi = 300)

ggsave("length_dit_fn_cp.pdf", device = "pdf", width = 8, height =  10,
       path = "./", plot = cpm_plot, dpi = 300)
