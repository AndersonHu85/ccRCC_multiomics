################# Figure 4a ######################

genelist_Figure4a <- read.table("genelist_Figure4a.txt")[,1]###pooled up-regulated and down-regulated genes in IM4

normal <- grep("_N",colnames(mat),value = T)
top_anno <- HeatmapAnnotation(df =sample_info[c(normal,tumor),c(1:2,6:19)],
                              col = cl_col,show_legend = T)

heat <- Heatmap(mat[genelist_Figure4a,c(normal,tumor)], 
                col = colorRamp2(seq(-2,2,length.out = 11),col),
                cluster_rows = T,
                cluster_columns = T,
                show_column_names = F,
                show_row_names = F,
                use_raster=F,
                show_row_dend = T,
                na_col = "grey85",
                clustering_distance_rows = "pearson",
                clustering_method_rows = "ward.D",
                clustering_distance_columns = "pearson",
                clustering_method_columns = "ward.D",
                row_split = 4,
                column_split = c((rep("normal",50)),rep("tumor",100)),
                column_gap = unit(1, "mm"),
                row_gap = unit(1, "mm"),
                top_annotation = top_anno,
                column_title = NULL) # 不需要列标题)
heat

a <- row_order(heat)

library(clusterProfiler)
library(org.Hs.eg.db)

go_UP <- enrichGO(genelist_Figure4a[a[[4]]], #extract genes in module 4,the same process was performed to obtain module 3
                  OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.01, 
                  qvalueCutoff = 0.05,keyType = 'SYMBOL')

########################### Figure 4B ##########################

UP_model3 <- read.table("./Figure 4/module3.txt", sep = "\t", header = T, row.names = NULL)

UP_model3$p.adjust <- -log10(UP_model3$p.adjust)
UP_model3 <- UP_model3[order(UP_model3$p.adjust, decreasing = T),]
UP_model3$order <- c(nrow(UP_model3):1)
p1 <- ggplot(data=UP_model3, aes(x=reorder(Description,order),y=p.adjust, fill=p.adjust)) + 
  geom_bar(position=position_dodge(), stat="identity") +   
  xlab("GO term") + ylab("-log10(adjusted.p-value)") +coord_flip() + theme_bw() 
p1 + scale_fill_gradientn(colours = brewer.pal(9,"Blues")[-1],limit=c(2,10))


UP_model4 <- read.table("./module4.txt", sep = "\t", header = T, row.names = NULL)

UP_model4$p.adjust <- -log10(UP_model4$p.adjust)
UP_model4 <- UP_model4[order(UP_model4$p.adjust, decreasing = T),]
UP_model4$order <- c(nrow(UP_model4):1)
p1 <- ggplot(data=UP_model4, aes(x=reorder(Description,order),y=p.adjust, fill=p.adjust)) + 
  geom_bar(position=position_dodge(), stat="identity") +   
  xlab("GO term") + ylab("-log10(adjusted.p-value)") +coord_flip() + theme_bw() 
p1 + scale_fill_gradientn(colours = brewer.pal(9,"Reds")[-1],limit=c(3,8))


############################# metabolomics ##############################
# calculate mean value of metabolites in each IM subgroups

meta_obj <- CreateSeuratObject(matrix_meta[,c(normal,tumor)])
meta_obj$group <- c(rep("N",50),sample_info[tumor,]$class)
meta_obj <- SetIdent(meta_obj,value = "group")

mean_mat <- as.matrix(AverageExpression(meta_obj, return.seurat = T)@assays$RNA@data) #use seurat to enable rapid visuliztion of metabolites in each group

IM4_meta <- c(
  "Argininosuccinic acid","L-Arginine","Oxidized glutathione","Glutathionate(1-)","Fumaric acid",
  "Aminopropylcadaverine","Spermidine","Spermine","Glycine",
  "L-Ornithine","L-Aspartic acid","Urea","N-Acetylaspartylglutamic acid","N-acetylaspartate","Gamma-aminobutyric acid"
)

data <- matrix_meta[,c(normal,tumor)] %>% as.matrix()
groups <- c(rep("N",50),sample_info$class)  

p_values <- matrix(NA, nrow = nrow(data), ncol = 5) 

for (i in 1:5) {
  p_value <- apply(data, 1, function(x) {
    wilcox.test(x[groups == unique(groups)[i]], x[groups != unique(groups)[i]])$p.value
  })
  p_value <- p.adjust(p_value,"BH")
  p_values[,i] <- p_value
  
}
colnames(p_values) <- unique(groups)
rownames(p_values) <- rownames(data)

b <- p_values[IM4_meta,]

b = ifelse(b >= 0.05, "",
           ifelse(b<0.0001,"****",
                  ifelse(b<0.001,"***",
                         ifelse(b<0.01,"**",
                                ifelse(b < 0.05,"*","")))))

bk=c(seq(-2,0,by=0.01),seq(0,2,by=0.01))

################### Figure 4F ######################
pheatmap(t(scale(t(mean_mat[IM4_meta,]))), border_color = "NA", fontsize = 9,cellheight = 15,cellwidth = 20,
         display_numbers = b,number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))



# calculate mean value of metabolites in each IM subgroups
data_df <- t(matrix[,c(normal,tumor)]) %>% as.data.frame()

data_df$group <- c(rep("N",50),sample_info[tumor,]$class)

group_means <- aggregate(. ~ group, data_df, mean)

result_mat <- data.frame(group_means[-1]) %>% t()
colnames(result_mat) <- group_means$group
rownames(result_mat) <- rownames(matrix_meta)

Ion <- c("SLC7A5","SLC17A9","SLC38A5","SLC18A3","SLC12A8","SLC35F3","SLC38A2",
         "SLC1A5","SLC2A1","SLC2A3","SLC2A2","SLC4A3","SLC6A17")

data <- matrix[,c(normal,tumor)] %>% as.matrix()
groups <- c(rep("N",50),sample_info$class)  

p_values <- matrix(NA, nrow = nrow(data), ncol = 5) 

for (i in 1:5) {
  p_value <- apply(data, 1, function(x) {
    wilcox.test(x[groups == unique(groups)[i]], x[groups != unique(groups)[i]])$p.value
  })
  p_value <- p.adjust(p_value,"BH")
  p_values[,i] <- p_value
  
}
colnames(p_values) <- unique(groups)
rownames(p_values) <- rownames(data)

b <- p_values[Ion,c(2:5,1)]

b = ifelse(b >= 0.05, "",
           ifelse(b<0.0001,"****",
                  ifelse(b<0.001,"***",
                         ifelse(b<0.01,"**",
                                ifelse(b < 0.05,"*","")))))

bk=c(seq(-2,0,by=0.01),seq(0,2,by=0.01))

################### Figure 4F ######################
pheatmap(t(scale(t(result_mat[Ion,]))), border_color = "NA", fontsize = 9,cellheight = 15,cellwidth = 20,
         display_numbers = b[,],number_color="black",fontsize_number=10,
         clustering_distance_rows="correlation",
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))

