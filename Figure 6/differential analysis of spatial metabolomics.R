######## differential analysis was performed in pos and neg data separately and then combimed into a single file
library(ggrepel)

volcano <- read.csv("./Figure 6/diff_Y7.csv",row.names = 1)

volcano$p <- (-1)*log10(volcano$p_val_adj)
volcano$p[volcano$p == "Inf"] = 303
volcano$color = ifelse(volcano$p_val_adj > 0.001, "stable",
                       ifelse(volcano$avg_log2FC > 1, "up",
                              ifelse(volcano$avg_log2FC < -1,"down","stable")))
volcano$color2 <- "others"
volcano$color2[(grep("Fatty Acyls",volcano$Class))]="Fatty Acyls"
volcano$color2[(grep("Carboxylic acids and derivatives",volcano$Class))]="Carboxylic acids"
volcano$color2[(grep("Glycerolipids",volcano$Class))]="Glycerolipids"
volcano$color2[volcano$color=="stable"] <- "others"
volcano <- arrange(volcano,desc(volcano$color2))

################# Figure 6M ######################
p <- ggplot(
  volcano, aes(x = avg_log2FC, y = p, colour=color2, shape=group)) +
  geom_point(alpha=0.9, size=2) +
  scale_color_manual(values=c(pal_npg()(4)[c(3,1,4)], "#d2dae2"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

volcano$label <- NA
volcano$label[intersect(grep("Fatty Ac",volcano$Class),grep("up",volcano$color))]=volcano$gene[intersect(grep("Fatty Ac",volcano$Class),grep("up",volcano$color))]
volcano$label[intersect(grep("Triradylcglycerols",volcano$Sub.Class),grep("up",volcano$color))]=volcano$gene[intersect(grep("Triradylcglycerols",volcano$Sub.Class),grep("up",volcano$color))]
volcano$label <- as.numeric(volcano$label)
volcano$label <- signif((volcano$label+0.00001),digits = 7)

p +ylim(0,320)+ geom_label_repel(data = volcano, aes(x = volcano$avg_log2FC, 
                                                     y = p,
                                                     label = label),
                                 #max.overlaps = 100,
                                 size = 3, box.padding = unit(0.5, "lines"),
                                 point.padding = unit(0.8, "lines"),
                                 segment.color = "black",
                                 color = "black",
                                 show.legend = F)
