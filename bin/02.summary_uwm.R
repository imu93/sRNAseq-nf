# Summary of UWM alignment results
# Load packages and siRmap logs

# Note: a general improvement would be to use functions to avoid code repetition
# also I need a dinamic text size based on the number of samples

# Actually at two levels: 
# 1) general size
# 2) label size

# consider this for next version

pacman::p_load(ggplot2, pals, scales, ggpubr, reshape2)
logs = list.files(pattern = ".*.log")
files = lapply(logs, readLines)
names(files) = sub("\\..*", "", logs)

# buold data frame for plotting
df2plot = list()
for (i in names(files)) {
  x = files[[i]]
  # helper function to extract numeric values from log lines
  get_num = function(key) {
    as.numeric(gsub("[ ,]", "", sub(".*:", "", x[grepl(key, x)])))
  }
  uq  = get_num("Expanded unique:")
  mm  = get_num("Expanded multimappers:")
  ave_prob  = get_num("Mean UWM confidence:")
  total = get_num("Reads written")
  #Build data frame
  df = data.frame("ID"=i, "Unique_M"=uq, "Multimappers_M"=mm, "Ave_prob_uwm"=ave_prob, "Total"=total)
  # Estimate percentages
  df$Percentage_unique = round((100 * (df$Unique / unique(total))), 2)
  df$Percentage_multimappers = round((100 * (df$Multimappers / unique(total))), 2)
  df2plot[[i]] = df
}
# Combine and melt data frame for plotting
names(df2plot) = NULL
df2plot = do.call(rbind, df2plot)
df = df2plot
df2plot = melt(df2plot)

# Subset data for plotting
df1 = df2plot[grepl("Unique_M|Multimappers_M", df2plot$variable),]
df1$Class = factor(df1$variable, levels = c("Multimappers_M","Unique_M"))
color = rev(c("#33B7A0", "#1283D4"))
df1$exp = "Numebr of aligned reads"
# This can be improved by using the same theme as in other plots
# I'll do it later
p1 = ggplot(df1, aes(value, ID, fill = Class)) + 
  geom_bar(stat = "identity") +
  theme_test() +  xlab("Number of reads") + ylab("Sample ID") +
  scale_fill_manual(values = color, name = "Aln type") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(labels = unit_format(unit = "M",  scale = 1e-6))  +
  theme(axis.title.x = element_text(face = "bold", size = 18), 
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18)) + facet_wrap(~exp) + 
  theme( strip.text = element_text(size = 18, colour = "black")) +
  theme(legend.text = element_text(size = 16))+
  guides(fill = guide_legend(reverse = TRUE), size=4)

# Same for percentage plot
df2 = df2plot[grepl("Perc", df2plot$variable),]
df2$Class = ifelse(grepl("uniq", df2$variable), "Unique", "Multimappers")
df2$Class = factor(df2$Class, levels = c("Multimappers","Unique"))
df2$exp = "Numebr of aligned reads"
# Again, can be improved later
p2 = ggplot(df2, aes(value, y= ID, fill=Class)) + 
  geom_bar(stat = "identity") +
  theme_test() +  xlab("% of reads") + ylab("Sample ID") +
  scale_fill_manual(values = color, name = "Aln type") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(labels = unit_format(unit = "%")) +
  theme(axis.title.x = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18)) + facet_wrap(~exp) + 
  theme( strip.text = element_text(size = 18, colour = "black")) +
  theme(legend.text = element_text(size = 16))+
  guides(fill = guide_legend(reverse = TRUE), size=4)

p3 = ggarrange(p1, p2, ncol = 1, common.legend = T, legend = "bottom")
# Save summary table
options(sicpen=999)
df = df %>% mutate(df, across(c("Unique_M", "Multimappers_M", "Total"), ~ ./1e6))
df = df %>% mutate(df, across(c("Unique_M", "Multimappers_M", 
                                "Total", "Ave_prob_uwm"), ~ round(., 2)))
df = df[,c(1,2,3,6,7,4,5)]
write.table(df, "Alignment_summary.txt", sep = "\t", quote = F, 
            col.names = T)
# adn save plot
ggsave("Asignment_summary.png", device = "png", width = 12, height = 20,
       dpi = 200, path = "./", plot = p3)
