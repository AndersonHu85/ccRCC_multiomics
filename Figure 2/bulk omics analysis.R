library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(ggplot2)

sample_info <- read.csv("./sample_info.csv", row.names = 1, colClasses = "character") ####included in Table S1
tmp_org <- read.csv("./tmp_org.csv",row.names = 1)#####signature gene of four IM subtypes, included in Table S3

tumor <- grep("_T",colnames(coding_combat2_mat),value = T)
all(rownames(sample_info)==tumor)

mat_tumor_only <- coding_combat2_mat[,tumor]%>%t()%>%scale()%>%t()
mat_tumor_only <- na.omit(mat_tumor_only)

cl_col =list(age=c("<40"="#edf8b1","40-60"="#7fcdbb",">60"="#2c7fb8"),
             sex=c("male"="#66C2A5", "female"="#FC8D62","M"="#66C2A5", "F"="#FC8D62"),
             #T=c("T1"="#f1eef6","T2"="#d7b5d8","T3"="#df65b0","T4"="#ce1256"),
             #N=c("N0"="#bcbddc","N1"="#756bb1"#,"pNx"="#9B808F50"),
             grade=c("1"="#fee5d9","2"="#fcae91","3"="#de2d26","4"="#a50f15"),
             stage = c("stage I"="#C6DBEF", "stage II"="#6BAED6", "stage III"="#2171B5", "stage IV"="#08306B"),
             class = c("1"="#8DD3C7", "2"= "#BEBADA", "3"= "#FB8072", "4"= "#80B1D3"),
             X3p = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             X5q = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             X7p = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             X9p = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             X12p = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             X14q = c("WT"="grey90","gain"="#de2d26","loss"="#2171B5","NA"="grey"),
             VHL=c("mut"="#50b9e3","WT"="grey90"),PBRM1=c("mut"="#50b9e3","WT"="grey90"),
             BAP1=c("mut"="#50b9e3","WT"="grey90"),
             SETD2=c("mut"="#50b9e3","WT"="grey90"),KDM5C=c("mut"="#50b9e3","WT"="grey90")
             )

top_anno <- HeatmapAnnotation(df =sample_info[,c(1:2,6:19)],
                              col = cl_col,show_legend = T)


################ Figure 2A, part 1 #####################
heat <- Heatmap(mat_tumor_only[unique(tmp_org$Symbol),tumor], 
                col = colorRamp2(seq(-2,2,length.out = 11),col),
                cluster_rows = T,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = F,
                use_raster=F,
                show_row_dend = T,
                na_col = "grey85",
                #clustering_distance_rows = "pearson",
                clustering_method_rows = "ward.D",
                column_split = sample_info$class,
                column_gap = unit(1, "mm"),
                row_gap = unit(0, "mm"),
                top_annotation = top_anno,
                column_title = NULL) 
heat

################ Figure 2A, part 2 #####################

Fig2_path <- c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","GOMF_COLLAGEN_BINDING",
               "REACTOME_SMOOTH_MUSCLE_CONTRACTION","GOBP_SPROUTING_ANGIOGENESIS","KEGG_FATTY_ACID_METABOLISM",
               "GOBP_LIPID_OXIDATION","GOBP_T_CELL_ACTIVATION","GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
               "GOBP_ACUTE_INFLAMMATORY_RESPONSE","KEGG_COMPLEMENT_AND_COAGULATION_CASCADES","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_ANGIOGENESIS","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
               "GOBP_MITOCHONDRIAL_GENE_EXPRESSION")

ssgsea_RNA_Fi2a_path <- rbind(res_h_RNA,res2_RNA,res5_RNA)[Fig2_path,tumor]
ssgsea_RNA_Fi2a_path_zs <- ssgsea_RNA_Fi2a_path[,tumor] %>% t() %>% scale() %>% t()
rownames(ssgsea_RNA_Fi2a_path_zs) <- tolower(rownames(ssgsea_RNA_Fi2a_path_zs))
rownames(ssgsea_RNA_Fi2a_path_zs) <-  strsplit(rownames(ssgsea_RNA_Fi2a_path_zs),"_")[] %>% 
  lapply(., function(x) x[-1]) %>% 
  sapply(., function(x) paste(x, collapse = " "))

Heatmap(ssgsea_RNA_Fi2a_path_zs[,tumor], 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        clustering_method_rows = "ward.D",
        column_split = cluster_info2[tumor,]$class,
        column_gap = unit(1, "mm"),
        top_annotation = top_anno, 
        column_title = NULL)

################ Figure 2A, part 3 #####################

ssgsea_protein_Fi2a_path <- rbind(res_h_protein,res2_protein,res5_protein)[Fig2_path,tumor]
ssgsea_protein_Fi2a_path_zs <- ssgsea_Fi2a_path[,tumor] %>% t() %>% scale() %>% t()
rownames(ssgsea_protein_Fi2a_path_zs) <- tolower(rownames(ssgsea_protein_Fi2a_path_zs))
rownames(ssgsea_protein_Fi2a_path_zs) <-  strsplit(rownames(ssgsea_protein_Fi2a_path_zs),"_")[] %>% 
  lapply(., function(x) x[-1]) %>% 
  sapply(., function(x) paste(x, collapse = " "))

Heatmap(ssgsea_protein_Fi2a_path_zs[,tumor], 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        clustering_method_rows = "ward.D",
        column_split = cluster_info2[tumor,]$class,
        column_gap = unit(1, "mm"),
        top_annotation = top_anno, 
        column_title = NULL)



################ Figure 2A, part 4 #####################

sig_911 <- fread("911_sig.txt") %>% as.data.frame() ######### 911 gene signature extracted from  single cell sequencing

sig_911 <- lapply(split(sig_911$gene, sig_911$cluster), as.character)

ssgsea_cell_signature_TJRCC <- gsva(as.matrix(coding_combat2_mat), sig_911, method = "ssgsea", 
                                    parallel.sz = 0, verbose = T)

ssgsea_cell_signature_TJRCC_zs <- ssgsea_cell_signature_TJRCC[,tumor] %>% t() %>% scale() %>% t()

Heatmap(ssgsea_cell_signature_TJRCC_zs[,tumor], 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        clustering_method_rows = "ward.D",
        column_split = cluster_info2[tumor,]$class,
        column_gap = unit(1, "mm"),
        top_annotation = top_anno, 
        column_title = NULL) 

################ Figure 2A, part 5 #####################

CIBER <- read.table("./xCell_mRNA_mat_xCell_combat.txt",sep = "\t",header = T, row.names = 1) #####calculated by online tool of xCELL
CIBER <- CIBER[c(3,5,6,7,8,10:12,16,18,19,21,23,31:34,39,42,46:47,52:53,56:58,60:61,65:67),tumor] #####remove cell types not associated with ccRCC

CIBER <- na.omit(t(scale(t(CIBER[,tumor]))))

Heatmap(CIBER[,tumor], 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        clustering_method_rows = "ward.D",
        column_split = cluster_info2[tumor,]$class,
        column_gap = unit(1, "mm"),
        top_annotation = top_anno, 
        column_title = NULL)



arm_level_cnv <- read.table("arm_level_cnv_TJRCC.txt")
cnv_counts <- colSums(arm_level_cnv != 0) %>% data.frame()
colnames(cnv_counts) <- "cnv_count"
cnv_counts$sample <- rownames(cnv_counts)
cnv_counts$IM_group <- sample_info[cnv_counts$sample,]$class %>% as.factor()
cnv_counts$cnv_count  <-  as.numeric(cnv_counts$cnv_count)

# on-way ANOVA analysis
aov_res <- aov(cnv_count ~ IM_group, data = cnv_counts)
summary(aov_res)

# If one-way ANOVA p-value < 0.05, then perform pairwise comparison with Tukey HSD test
tukey_res <- TukeyHSD(aov_res)
tukey_df <- data.frame(tukey_res$IM_group)

signif_df <- data.frame(
  xstart = as.numeric(gsub("-.*", "", rownames(tukey_df))),
  xend = as.numeric(gsub(".*-", "", rownames(tukey_df))),
  annotations = ifelse(tukey_df$p.adj >= 0.05, "",
                      ifelse(tukey_df$p.adj<0.001,"***",
                             ifelse(tukey_df$p.adj<0.01,"**",
                                    ifelse(tukey_df$p.adj < 0.05,"*",""))))  # text of the label
)
signif_df <- signif_df[signif_df$annotations!="",]

###################### Figure 2B ############################
ggplot(cnv_counts,aes(IM_group,cnv_count,fill=IM_group,color=IM_group))+
  geom_boxplot(outlier.colour = NA,size=1)+geom_jitter(size=1.5,aes(color=IM_group))+theme_bw()+
  scale_fill_manual(values = alpha(class,0)) + scale_color_manual(values = class)+ ylim(0,40) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="Normalized level") + 
  theme(
    panel.border = element_rect(colour = "black",size=1), 
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  ) + ggtitle("CNV burden")+
  geom_signif(annotations = signif_df$annotations,
              y_position = c(33,35,37,39), 
              xmin = signif_df$xstart, 
              xmax = signif_df$xend,color="black")



SNA_number <- annovar.laml_demo@variants.per.sample %>% 
  as.data.frame() %>%
  transform(Tumor_Sample_Barcode = toupper(sub("T$", "_T", Tumor_Sample_Barcode)))
row.names(SNA_number) <- SNA_number$Tumor_Sample_Barcode
SNA_number <- SNA_number[tumor,]
SNA_number$Tumor_Sample_Barcode <- tumor
SNA_number$Variants <- ifelse(is.na(SNA_number$Variants), 0, SNA_number$Variants)
SNA_number$IM_group <- sample_info[tumor,]$class %>% as.factor()

# on-way ANOVA analysis
aov_res <- aov(Variants ~ IM_group, data = SNA_number)
summary(aov_res)

###################### Figure 2C ############################
ggplot(SNA_number,aes(IM_group,Variants,fill=IM_group,color=IM_group))+
  geom_boxplot(outlier.colour = NA,size=1)+geom_jitter(size=1.5,aes(color=IM_group))+theme_bw()+
  scale_fill_manual(values = alpha(class,0)) + scale_color_manual(values = class)+ ylim(0,300) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="Normalized level") + 
  theme(
    panel.border = element_rect(colour = "black",size=1), 
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  ) + ggtitle("TMB")



########################### Figure 2D-1 ###############################
column <- as.data.frame(table(sample_info$class, sample_info$grade))

colnames(column) <- c("IM_group","grade","Freq")

grade=c("1"="#fee5d9","2"="#fcae91","3"="#de2d26","4"="#a50f15")
p1 <- ggplot(column,aes(IM_group,weight=Freq,fill=grade))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = grade) #[c(1,5,17,14,21)]
p1 + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################### Figure 2D-2 ###############################
column <- as.data.frame(table(sample_info$class, sample_info$stage))

colnames(column) <- c("IM_group","stage","Freq")

stage = c("stage I"="#C6DBEF", "stage II"="#6BAED6", "stage III"="#2171B5", "stage IV"="#08306B")
p1 <- ggplot(column,aes(IM_group,weight=Freq,fill=stage))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = stage) #[c(1,5,17,14,21)]
p1 + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################# Extended Data Figure 2F-2 ##################

a <- unique(as.character(sample_info[tumor,]$class))

##### differential analysis of pathway level ssgsea value between IM subtypes #####
results <- lapply(a, function(m) {
  group <- ifelse(sample_info[tumor,]$class == m, "select", "ref")
  group_mat <- model.matrix(~factor(group))
  colnames(group_mat) <- c("ref", "select")
  fit <- lmFit(rbind(res_h_RNA, res2_RNA, res5_RNA)[, tumor], group_mat)
  fit2 <- eBayes(fit)
  allDiff <- topTable(fit2, adjust = 'fdr', coef = 2, number = 200000)
  allDiff$X <- row.names(allDiff)
  ##### only select T-value from the result ############
  allDiff <- allDiff[, c("X", "t")]
  colnames(allDiff)[2] <- paste0("IM",m)
  return(allDiff)
})

b <- Reduce(function(x, y) {merge(x, y, by = "X")}, results)
b <- data.frame(b[-1],row.names = b$X)

gsva_T_RNA <- b[Fig2_path,]
rownames(gsva_T_RNA) <- tolower(rownames(gsva_T_RNA))
rownames(gsva_T_RNA) <-  strsplit(rownames(gsva_T_RNA),"_")[] %>% 
  lapply(., function(x) x[-1]) %>% 
  sapply(., function(x) paste(x, collapse = " "))

bk=c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
pheatmap(gsva_T_RNA, border_color = NA, fontsize = 9,cellheight = 15,cellwidth = 20,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = brewer.pal(11,"RdBu")[11:6])(length(bk)/2),
                   colorRampPalette(colors = brewer.pal(11,"RdBu")[6:1])(length(bk)/2)))


################# Extended Data Figure 2F-2 ##################

a <- unique(as.character(sample_info[tumor,]$class))

results <- lapply(a, function(m) {
  group <- ifelse(sample_info[tumor,]$class == m, "select", "ref")
  group_mat <- model.matrix(~factor(group))
  colnames(group_mat) <- c("ref", "select")
  fit <- lmFit(rbind(res_h_protein, res2_protein, res5_protein)[, tumor], group_mat)
  fit2 <- eBayes(fit)
  allDiff <- topTable(fit2, adjust = 'fdr', coef = 2, number = 200000)
  allDiff$X <- row.names(allDiff)
  ##### only select T-value from the result ############
  allDiff <- allDiff[, c("X", "t")]
  colnames(allDiff)[2] <- paste0("IM",m)
  return(allDiff)
})

b <- Reduce(function(x, y) {merge(x, y, by = "X")}, results)
b <- data.frame(b[-1],row.names = b$X)

gsva_T_protein <- b[Fig2_path,]
rownames(gsva_T_protein) <- tolower(rownames(gsva_T_protein))
rownames(gsva_T_protein) <-  strsplit(rownames(gsva_T_protein),"_")[] %>% 
  lapply(., function(x) x[-1]) %>% 
  sapply(., function(x) paste(x, collapse = " "))

bk=c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
pheatmap(gsva_T_protein, border_color = NA, fontsize = 9,cellheight = 15,cellwidth = 20,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = brewer.pal(11,"RdBu")[11:6])(length(bk)/2),
                   colorRampPalette(colors = brewer.pal(11,"RdBu")[6:1])(length(bk)/2)))


############# relationship between AA and IM subgroups #############
#################### Extended Data Figure 2G #######################

TJ_AA_sample <- SNA$Tumor_Sample_Barcode[SNA$group=="TJ"&SNA$SBS=="AA"] %>% sub("T$", "_T", .)
sample_info$AA <- "non-AA"
sample_info$AA[rownames(sample_info)%in%TJ_AA_sample]="AA"
column <- as.data.frame(table(sample_info$class, sample_info$AA))

network <- column
colnames(network) <- c("class", "AA", "num")
factor_list <- sort(unique(c(levels(network$class), levels(network$AA))))
num_list <- 0:(length(factor_list)-1)
levels(network$class) <- num_list[factor_list %in% levels(network$class)]
levels(network$AA) <- num_list[factor_list %in% levels(network$AA)]
network$class <- as.numeric(as.character(network$class))
network$AA <- as.numeric(as.character(network$AA))
attribute <- data.frame(name=c(factor_list))

color_scale <- 'd3.scaleOrdinal() 
.domain(["1","2","3","4","AA","non-AA"]) 
.range(["#8DD3C7","#BEBADA", "#FB8072", "#80B1D3","#D24B27B2", "#6E4B9EB2"])'

library(networkD3)
sn <- sankeyNetwork(Links = network, Nodes = attribute,
                    Source = "class", Target = "AA",
                    Value = "num", NodeID = "name", units = "TWh",colourScale = color_scale,
                    fontSize= 0, nodeWidth = 100)
sn

########### correlation between signature of cell subgroups in TJ-RCC ###############
specfic_gene_sig <- fread("./specific_gene_sig.txt")%>%as.data.frame() ###Table S4

pathway_list <- vector("list",length(specfic_gene_sig))

for (i in seq_along(specfic_gene_sig)) {
  pathway_list[[i]] <- unique(na.omit(specfic_gene_sig[,i]))
}

names(pathway_list) <- colnames(specfic_gene_sig)


pathway_list <- lapply(specfic_gene_sig, function(x) {
  unique(na.omit(x)) 
})

ssgsea_specfic_gene_sig_TJRCC <- gsva(as.matrix(coding_combat2_mat), pathway_list, method = "ssgsea", 
                                    parallel.sz = 0, verbose = T)

COR <- cor(t(ssgsea_specfic_gene_sig_TJRCC[,tumor]),method = "spearman")


bk=c(seq(-1,-0.01,by=0.01),seq(0,1,by=0.01))
pheatmap(COR, border_color = "NA", fontsize = 9,cellheight = 15,cellwidth = 15,
         clustering_distance_rows="correlation", clustering_distance_cols ="correlation",
         cluster_col=T, cluster_rows=T, border= NULL, breaks=bk, 
         treeheight_row = 20,treeheight_col = 20, #display_numbers =T,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))
