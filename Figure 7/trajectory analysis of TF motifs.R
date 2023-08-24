Epi_trajectory <- Epi_snRNA[,Epi_snRNA$nUMI>median(Epi_snRNA$nUMI)]#only use top 50% to minimize the influence of sequencing depth

library(monocle3)
featureData <- as.data.frame(Epi_trajectory@assays$RNA@counts@Dimnames[[1]])
colnames(featureData) <- "gene_short_name"

rownames(featureData) <- featureData[,1]
cds <- new_cell_data_set((Epi_trajectory@assays$RNA@data),
                         cell_metadata = Epi_trajectory@meta.data,
                         gene_metadata = featureData)
cds


#cds <- preprocess_cds(cds, num_dim = 50)
#cds <- align_cds(cds, alignment_group = "sample")

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- as.matrix(Epi_trajectory@reductions[["umap"]]@cell.embeddings)

cds@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- as.matrix(Epi_trajectory@reductions[["pca"]]@cell.embeddings)
cds@int_colData@listData[["reducedDims"]]@listData[["Aligned"]] <- as.matrix(Epi_trajectory@reductions[["harmony"]]@cell.embeddings)

library(destiny)

dm <- DiffusionMap(Epi_trajectory@reductions$harmony@cell.embeddings[,1:15],
                   k = 3,n_pcs = NA, verbose = TRUE, n_eigs=2)#obtain diffusion map dimension

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,1] <- dm@eigenvectors$DC1
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,2] <- dm@eigenvectors$DC2

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "predicted_class")

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = c("partition"))

cds <- learn_graph(cds,close_loop=F)
plot_cells(cds,show_trajectory_graph=T,
           group_label_size = 4,
           color_cells_by = "sample",
           cell_size = 0.7,
           label_cell_groups=FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE) + scale_color_manual(values = celltype)
cds <- order_cells(cds)

################# Figure 7A ##################
plot_cells(cds,show_trajectory_graph=F,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           cell_size = 0.7,
           label_branch_points=FALSE,
           graph_label_size=1.5)

Epi_trajectory$monocle3_pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
meta <- Epi_trajectory@meta.data

################# Figure 7B ##################
ggplot(meta, aes(monocle3,Cluster1)) + #Cluster1 represents IM2 score from AddModuleScore
  geom_point(aes(color = predicted_class, alpha=0.9),size=0.05,meta) +
  geom_smooth(color="black",se = T) + scale_fill_npg()+
  theme_classic() + 
  scale_color_manual(
    values = c("#4DBBD5","#E64B35")
  )


##################### draw pseudoheatmap with monocle 2 ####################
library(monocle)

seq <- seq(from=1,to=ncol(Epi_trajectory),by=7)
metadata <- Epi_trajectory@meta.data[seq,]
featureData <- as.data.frame(Epi_trajectory@assays$RNA@counts@Dimnames[[1]])
rownames(featureData) <- featureData[,1]
pd <- new("AnnotatedDataFrame", data = metadata)
fd <- new("AnnotatedDataFrame", data = featureData)

monocle <- new("CellDataSet", exprs=as.matrix(Epi_trajectory@assays$RNA@data[,seq]), phenoData = pd,
               expressionFamily = negbinomial(),lowerDetectionLimit=0.1)
monocle

monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)

monocle@phenoData@data[["Pseudotime"]]=monocle@phenoData@data$monocle3_pseudotime

newdata <- data.frame(Pseudotime = seq(min(monocle2$pseudotime), 
                                       max(monocle2$pseudotime), length.out = 1000))

TF_final <- read.csv("Figure 6/TF_final.csv")
cds_Epi_subset <- monocle[c(unique(TF_final$symbol)),]
genSmoothCurves_mat<-genSmoothCurves(cds_Epi_subset,
                                     new_data = newdata,
                                     trend_formula = "~sm.ns(pseudotime, df = 3)",
                                     cores = 1)

order <- order(monocle2$pseudotime,decreasing = F)

mark_gene <- c("ZBTB7C","CEBPB","HOXD10","HIF1A","TFE3","PPARA","TFEC","HNF1B",
               "HNF4A","HNF1A","AR","ZNF711","EMX2","NEF2L2","FOSL1","JUNB","FOSL2")
mark_gene <- TF_final$symbol[TF_final$symbol%in%mark_gene]
gene_pos <- which(TF_final$symbol %in% mark_gene)

mark_gene <- rownames(genSmoothCurves_mat)[rownames(genSmoothCurves_mat)%in%mark_gene]
gene_pos <- which(rownames(genSmoothCurves_mat) %in% mark_gene)
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = mark_gene))

cols <- brewer.pal(11,"RdBu")[11:1]

####################### Figure 7C ############################
heat <- Heatmap(t(scale(t(genSmoothCurves_mat[,order]))), 
                col = colorRamp2(seq(-2,2,length.out = 11),cols),
                cluster_rows = T,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = T,
                use_raster=F,
                border = NA,
                show_row_dend = T,
                na_col = "grey85",
                clustering_method_rows = "ward.D",
                column_gap = unit(1, "mm"),
                row_gap = unit(0, "mm"),
                row_dend_width = unit(20, "mm"),
                right_annotation = row_anno,
                column_title = NULL) 
heat
dev.off()

#### Figure 7G was drawn using the same pipeline


