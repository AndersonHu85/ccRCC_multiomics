Epi_snRNA <- snRNA[,snRNA$celltype2=="Malignant"]

Epi_snRNA <- NormalizeData(Epi_snRNA)

Epi_snRNA <- FindVariableFeatures(Epi_snRNA, nfeatures = 2000)
VariableFeaturePlot(Epi_snRNA)
Epi_snRNA <- ScaleData(Epi_snRNA, verbose = T)
Epi_snRNA <- RunPCA(Epi_snRNA, npcs = 50, verbose = FALSE)
PCAPlot(Epi_snRNA)
ElbowPlot(Epi_snRNA, ndims = 50)

Epi_snRNA <- RunHarmony(Epi_snRNA, c("sample"),max.iter.harmony = 20)


Epi_snRNA <- CellCycleScoring(Epi_snRNA, g2m.features = cc.genes$g2m.genes,
                          s.features = cc.genes$s.genes)


########################### Figure 2E #############################
Epi_snRNA  <- RunUMAP(Epi_snRNA, reduction = "harmony", dims = 1:40, n.components = 2)

DimPlot(Epi_snRNA, reduction = "umap",  label = F, repel = TRUE, raster= F,group.by = "sample",
        cols = celltype
) ##### UMAP dimension in Extended Data Figure 9C & 9D
celltype <- c(brewer.pal(11,"Set3")[c(1,3:8,10:11)],pal_npg(alpha = 0.7)(9),#color,
              brewer.pal(11,"Paired")[11:1],brewer.pal(8,"Set1")[8:1],
              brewer.pal(8,"Set2"))[-15]

Epi_snRNA  <- RunUMAP(Epi_snRNA, reduction = "pca", dims = 1:40, n.components = 2)

DimPlot(Epi_snRNA, reduction = "umap",  label = F, repel = TRUE, raster= F,group.by = "sample",
        cols = celltype
) ##### UMAP dimension in Figure 6D


Epi_snRNA <- AddModuleScore(Epi_snRNA,DCCD_score)
Epi_snRNA$predicted_class <- "IM4-like"
Epi_snRNA$predicted_class[Epi_snRNA$Cluster1>Epi_snRNA$Cluster2] <- "IM2-like" ##classify malignant cells into two groups

DimPlot(Epi_snRNA, reduction = "umap",  label = F, repel = TRUE, raster= F,group.by = "predicted_class",
        cols = celltype
)


################# Figure 6F ####################
column <- as.data.frame(table(Epi_snRNA$sample, Epi_snRNA$predicted_class))

column$group2 <- plyr::mapvalues(x = column$Var1, 
                                 from = rownames(sample_info), 
                                 to = sample_info$class)

colnames(column) <- c("Patient","cluster","Freq","group2")
column$group2 <- as.character(column$group2)

p1 <- ggplot(column,aes(Patient,weight=Freq,fill=cluster))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("IM4-like"="#E64B35FF", "IM2-like"="#4DBBD5FF")) 
p1 +  facet_grid(~group2,scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(.2, "cm"))



################# Extended Data Figure 9E ################
a <- as.data.frame(table(Epi_snRNA$sample, Epi_snRNA$predicted_class))
a <- spread(a,key = "Var2", value = "Freq")
a <- data.frame(a[2:ncol(a)], row.names = a$Var1)
a <- a/rowSums(a)
a <- as.data.frame(a)

a$group <- plyr::mapvalues(x = rownames(a), 
                                 from = rownames(sample_info), 
                                 to = sample_info$class)

a$group<- as.character(a$group)


class <- c("1"="#8DD3C7", "2"="#BEBADA", "3"="#FB8072", "4"="#80B1D3" )

ggplot(a,aes(factor(group,levels =c("1","2","3","4")),
             a$IM4.like,color=group))+
  geom_boxplot(outlier.colour = NA,size=1)+geom_jitter(size=1.5,aes(color=group))+theme_classic()+
  #scale_fill_manual(values = alpha(class,0)) +
  scale_color_manual(values = class)+
  geom_signif(comparisons = list(c("3","4"),c("2","4"),c("4","1")),
              map_signif_level=T,#size = 1,
              textsize=5,test=wilcox.test,step_increase=0.1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ labs(y="Proportion", x="IM group") + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.7), 
    axis.ticks = element_line(size = 0.7, color="black") ,
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12) 
  )+ scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.00))


############ CNV analysis of snRNA-seq data #####################

control <- snRNA_NAT[,snRNA_NAT$celltype2 == "PT"] ##use cell of origination of ccRCC to generate the reference of infercnv
Epi_sub$label <- Epi_sub$sample
control$label <- "control"

Epi_cnv <- merge(control, Epi_sub)

matrix <- Epi_cnv@assays$RNA@counts
meta <- Epi_cnv@meta.data
meta2 <- as.data.frame(meta$label)
table(meta$label)
rownames(meta2) <- rownames(meta)
colnames(meta2) <- "label"
class(meta2$label)

library(infercnv)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=matrix,
                                    annotations_file=meta2,
                                    delim="\t",
                                    gene_order_file=gtf_df,
                                    ref_group_names=c("control")) 

dir="./inferCNV"

gtf1 <- rtracklayer::import('dir/Homo_sapiens.GRCh38.87.chr.gtf') #load the reference genome downloaded from ensemble
gtf_df <- as.data.frame(gtf1)
gtf_df <- gtf_df[,c(1:3,12)]
seq <- as.character(gtf_df$seqnames)
seq <- paste("chr",seq,sep = "")
gtf_df$seqnames <- seq
gtf_df <- na.omit(gtf_df)
gtf_df <- gtf_df[!duplicated(gtf_df[,4]),]
gene <- rownames(Epi_cnv)
gtf_df <- gtf_df[gtf_df$gene_name %in% gene,]
rownames(gtf_df) <- gtf_df$gene_name
gtf_df <- gtf_df[,1:3]
gtf_df <- gtf_df[gtf_df$seqnames != "chrMT",] #only use chr 1-22


gtf_df$seqnames <- as.factor(gtf_df$seqnames)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir= dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             num_threads = 20,
                             HMM=F,
                             no_plot = T)

##### cnv matrix could be obtained from the run.final.infercnv_obj in the output folder




