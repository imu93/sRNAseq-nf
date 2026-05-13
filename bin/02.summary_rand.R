pacman::p_load(ggplot2, pals, scales, ggpubr)
logs = list.files(pattern = ".*.log")
files = lapply(logs, readLines)
names(files) = sub("\\..*", "", logs)

df2plot = list()
for (i in names(files)) {
  
  x = files[[i]]
  get_num = function(key) {
    as.numeric(gsub("[ ,]", "", sub(".*:", "", x[grepl(key, x)])))
  }
  
  uq  = get_num("Unique:")
  mm  = get_num("Multimappers:")
  #rt  = get_num("Rand")
  #mm  = mm-rt
  tot = get_num("Reads written")
  
  
  df = data.frame("ID"=i, "Class"=c("Unique", "Multimappers"),
                  "Counts"=c(uq, mm))
  df$Percentage = round((100 * (df$Counts / tot)), 2)
  df2plot[[i]] = df
}
names(df2plot) = NULL
df2plot = do.call(rbind, df2plot)
df2plot$Class = factor(df2plot$Class, 
                       levels = rev(c("Unique", "Multimappers")))

color = rev(c("#33B7A0", "#1283D4"))
df2plot$exp = "Number of aligned reads"
p1 = ggplot(df2plot, aes(Counts, ID, fill = Class)) + 
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
  guides(fill = guide_legend(reverse = TRUE))


p2 = ggplot(df2plot, aes(Percentage, y= ID, fill=Class)) + 
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
  guides(fill = guide_legend(reverse = TRUE))

p3 = ggarrange(p1, p2, ncol = 1, common.legend = T, legend = "bottom")

ggsave("Asignment_summary.png", device = "png", width = 12, height = 20,
       dpi = 200, path = "./", plot = p3)
