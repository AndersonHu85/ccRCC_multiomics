# calculate mean value of metabolites in each IM subgroups
data_df <- t(matrix_meta[,tumor]) %>% as.data.frame()

data_df$group <- sample_info[rownames(data_df),]$class

group_means <- aggregate(. ~ group, data_df, mean)

result_mat <- data.frame(group_means[-1]) %>% t()
colnames(result_mat) <- group_means$group
rownames(result_mat) <- rownames(matrix_meta)

IM3_meta <- c("Glucose 6-phosphate","D-fructose 6-phosphate",
         "Xanthosine","Guanosine","Guanosine 3'-phosphate","Inosine","IDP",
         "dCMP","Hypoxanthine",#"UDP-L-rhamnose","CDP-glucose",
         "LysoPC(20:4(8Z,11Z,14Z,17Z)/0:0)","LysoPC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)",
         "LysoPC(22:5(4Z,7Z,10Z,13Z,16Z)/0:0)","LysoPC(18:2(9Z,12Z))","PC(16:1(9Z)/0:0)",
         "PC(22:4(7Z,10Z,13Z,16Z)/0:0)"
)

row_anno <- data.frame(c(rep("glycolysis",2),rep("pyrimidine derivatives",7),rep("Glycerophospholipids",6)))
colnames(row_anno) <- "class"
rownames(row_anno) <- IM3_meta
cl_right <- list(class=c("glycolysis"="#8DD3C7","pyrimidine derivatives"="#BEBADA","Glycerophospholipids"="#FB8072"))
right_anno <- rowAnnotation(df=row_anno,col = cl_right,show_legend = T)

data <- matrix_meta[IM3_meta,tumor] %>% as.matrix()
groups <- sample_info$class  

p_values <- matrix(NA, nrow = nrow(data), ncol = 4) 

for (i in 1:4) {
  p_value <- apply(data, 1, function(x) {
    wilcox.test(x[groups == unique(groups)[i]], x[groups != unique(groups)[i]])$p.value
  })
  p_value <- p.adjust(p_value,"BH")
  p_values[,i] <- p_value
  
}
colnames(p_values) <- unique(groups)
rownames(p_values) <- rownames(data)

b <- p_values[IM3_meta,]

b = ifelse(b >= 0.05, "",
           ifelse(b<0.0001,"****",
                  ifelse(b<0.001,"***",
                         ifelse(b<0.01,"**",
                                ifelse(b < 0.05,"*","")))))

bk=c(seq(-2,0,by=0.01),seq(0,2,by=0.01))

################### Figure 3E ######################
pheatmap(t(scale(t(result_mat[IM3_meta,]))), border_color = "NA", fontsize = 9,cellheight = 15,cellwidth = 20,
         display_numbers = b,number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         annotation_row = row_anno,annotation_colors = cl_right,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))
