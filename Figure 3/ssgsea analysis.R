############extract 4 metabolic-associated GO modules #################
#################### Extended Data Figure 5A #######################

metabolic <- fread("./metabolic_selected.txt")

GO2 <- rbind(res_h_RNA,res2_RNA,res5_RNA)[metabolic$GO_terms,tumor]
GO2 <- t(scale(t(GO2)))

top_anno <- HeatmapAnnotation(df =sample_info[,c(1:2,6:19)],
                              col = cl_col,show_legend = T)

col <- brewer.pal(11,"RdBu")[11:1]
heat2 <- Heatmap(GO2[,], 
                 col = colorRamp2(seq(-2,2,length.out = 11),col),
                 cluster_rows = T,
                 cluster_columns = T,
                 show_column_names = F,
                 show_row_names = F,
                 use_raster=F,
                 border = "BLACK",
                 show_row_dend = T,
                 row_gap = unit(1, "mm"),
                 column_split = cluster_info2[tumor,]$class,
                 row_split = 4,
                 column_gap = unit(1, "mm"),
                 top_annotation = top_anno,
                 column_title = NULL)
heat2


#############select pyrimidine realated GO terms ###############
#################### Extended Data Figure 3A #######################

GO3 <- GO2[c(grep("NUCLEOSIDE",metabolic$GO_terms,value = T),grep("PYRIMIDINE",metabolic$GO_terms,value = T),
             grep("CYTIDINE",metabolic$GO_terms,value = T)),tumor]

GO3 <- GO3[!duplicated(rownames(GO3)),]
row.names(GO3) <- tolower(row.names(GO3))
row.names(GO3) <- gsub("_"," ",row.names(GO3))
row.names(GO3) <- gsub("gobp","",row.names(GO3))

heat2 <- Heatmap(GO3, 
                 col = colorRamp2(seq(-2,2,length.out = 11),col),
                 cluster_rows = T,
                 cluster_columns = T,
                 show_column_names = F,
                 show_row_names = T,
                 use_raster=F,
                 show_row_dend = T,
                 clustering_method_rows = "ward.D",
                 clustering_method_columns = "ward.D",
                 row_gap = unit(1, "mm"),
                 column_gap = unit(1, "mm"),
                 top_annotation = top_anno,
                 column_title = NULL)
heat2

############# comparison of RNA-seq and proteomics #############
######################### Figure 3B ############################

ssgsea_IM3_diff <- fread("./IM3_ssgsea_diff.txt") %>% as.data.frame()###load organized data

rownames(dotplot) <- tolower(rownames(dotplot))

pathway <- c(grep("gobp",rownames(dotplot),value = T),grep("kegg",rownames(dotplot),value = T),
             grep("reactome",rownames(dotplot),value = T),grep("hallmark",rownames(dotplot),value = T))
dotplot <- dotplot[pathway,]
dotplot$color <- "others"
dotplot$pathway <- rownames(dotplot)
ssgsea_IM3_diff$color[grep("nucleoside",ssgsea_IM3_diff$pathway)]="nucleoside"
ssgsea_IM3_diff$color[grep("cytidine",ssgsea_IM3_diff$pathway)]="cytidine"
ssgsea_IM3_diff$color[grep("pyrimidine",ssgsea_IM3_diff$pathway)]="pyrimidine"
ssgsea_IM3_diff <- ssgsea_IM3_diff[order(ssgsea_IM3_diff$color,decreasing = T),]
ggplot(ssgsea_IM3_diff,aes(ssgsea_IM3_diff$proteomics,ssgsea_IM3_diff$transcriptome,color=color))+
  geom_point(size=1)+ggthemes::theme_base()+
  scale_color_manual(values = c("#E64B35FF", "#3C5488FF","grey", "#00A087FF"))+ 
  ylim(-11.2,11.2)+ylab("T-value (transcriptome)")+xlab("T-value (proteome)")+
  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)


data_box <- rbind(res_h_RNA,res2_RNA,res5_RNA)[c("GOBP_GLUCOSE_6_PHOSPHATE_METABOLIC_PROCESS",
                              "WP_GLYCEROPHOSPHOLIPID_BIOSYNTHETIC_PATHWAY"),tumor] %>% t() %>% as.data.frame()
data_box$IM_subtype <- sample_info[tumor,]$class


###################### Figure 3F-1 ############################

# on-way ANOVA analysis
aov_res <- aov(GOBP_GLUCOSE_6_PHOSPHATE_METABOLIC_PROCESS ~ IM_subtype, data = data_box)
summary(aov_res)

# If one-way ANOVA p-value < 0.05, then perform pairwise comparison with Tukey HSD test
tukey_res <- TukeyHSD(aov_res)
tukey_df <- data.frame(tukey_res$IM_subtype)

signif_df <- data.frame(
  xstart = as.numeric(gsub("-.*", "", rownames(tukey_df))),
  xend = as.numeric(gsub(".*-", "", rownames(tukey_df))),
  annotations = ifelse(tukey_df$p.adj >= 0.05, "",
                       ifelse(tukey_df$p.adj<0.001,"***",
                              ifelse(tukey_df$p.adj<0.01,"**",
                                     ifelse(tukey_df$p.adj < 0.05,"*",""))))  # text of the label
)
signif_df <- signif_df[signif_df$annotations!="",]

ggplot(data_box,aes(IM_subtype,GOBP_GLUCOSE_6_PHOSPHATE_METABOLIC_PROCESS,fill=IM_subtype,color=IM_subtype))+
  geom_boxplot(outlier.colour = NA,size=1)+geom_jitter(size=1.5,aes(color=IM_subtype))+theme_bw()+
  scale_fill_manual(values = alpha(class,0)) + scale_color_manual(values = class)+ #ylim(0,40) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="Normalized level") + 
  theme(
    panel.border = element_rect(colour = "black",size=1), 
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  ) + ggtitle("Glucose-6-phosphate (G6P) metabolic process")+
  geom_signif(annotations = signif_df$annotations,
              y_position = c(0.24,0.26,0.28,0.3), 
              xmin = signif_df$xstart, 
              xmax = signif_df$xend,color="black")

###################### Figure 3F-2 was drawn with the same code ############################
