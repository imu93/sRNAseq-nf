pacman::p_load(dplyr, reshape2, ggplot2, pals, scales, ggpubr)
args = commandArgs(trailingOnly = TRUE)
inFile = readLines(args[1])
libs = inFile[grepl("readfile", inFile)]
libs = strsplit(libs, " ") %>% unlist()
libs = libs[2:length(libs)]
libs = sub(".ps.*", "", libs)

total_reads =  inFile[grepl("Unique mappers", inFile)]
total_reads =  total_reads %>% strsplit(" ") %>%
  lapply(function(x){x[5]}) %>%
  unlist() %>% as.numeric()

unique = inFile[grepl("Unique mappers", inFile)]
unique = unique %>% strsplit(" ") %>%
    lapply(function(x){x[3]}) %>%
    unlist() %>% as.numeric()


multi = inFile[grepl("Multi mappers", inFile)]
multi = multi %>% strsplit(" ") %>%
  lapply(function(x){x[3]}) %>%
  unlist() %>% as.numeric()

nomap = inFile[grepl("Non mappers", inFile)]
nomap = nomap %>% strsplit(" ") %>%
  lapply(function(x){x[3]}) %>%
  unlist() %>% as.numeric()




df = data.frame("Total"= total_reads,
                "Uniquemappers" = unique ,
                "Multimappers" = multi, "Nonmappers"= nomap)
rownames(df) = libs

Uper = 100*(df$Uniquemappers/df$Total) %>% round(., 4)
Mper = 100*(df$Multimappers/df$Total) %>% round(., 4)
Nper = 100*(df$Nonmappers/df$Total) %>% round(., 4)
df2 = data.frame("Uniquemappers" = Uper ,"Multimappers" = Mper, "Nonmappers"= Nper)
rownames(df2) = libs
df = df[,c(2:ncol(df))]
df$id = libs
df = melt(df)
df$id = factor(df$id, levels = libs)

color = c("#33B7A0", "#1283D4",  "#352A87")
df$exp = "Numebr of aligned reads"
p1 = df %>% ggplot(aes(x = value, y= id, fill=variable)) + geom_bar(stat = "identity") +
  theme_light() +  xlab("Number of reads") + ylab("Sample ID") +
  scale_fill_manual(values = color, name = "Aln type") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(labels = unit_format(unit = "M",  scale = 1e-6))  +
  theme(axis.title.x = element_text(face = "bold", size = 14), 
        axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10)) + facet_wrap(~exp) + 
  theme( strip.text = element_text(size = 14, colour = "black"))

df2$id = libs
df2 = melt(df2)
df2$id = factor(df2$id, levels = libs) 
df2$exp = "Percentage of aligned reads"
p2 = df2 %>% ggplot(aes(x = value, y= id, fill=variable)) + geom_bar(stat = "identity") +
  theme_light() +  xlab("% of reads") + ylab("Sample ID") +
  scale_fill_manual(values = color, name = "Aln type") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(labels = unit_format(unit = "%")) +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10)) + facet_wrap(~exp) + 
  theme( strip.text = element_text(size = 14, colour = "black"))



p3 = ggarrange(p1, p2, ncol = 1, common.legend = T, legend = "bottom")
  
ggsave("ShortStack_summary.pdf", device = "pdf", width = 12, height = 20,
       dpi = 200, path = "./", plot = p3)
