library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
options(stringsAsFactors = F)

##########count matrix was obtained from mendeley################

merge@meta.data$nUMI <- merge@meta.data$nCount_RNA
merge@meta.data$nGene <- merge@meta.data$nFeature_RNA


mito.genes <- grep(
  pattern = "^MT-",
  x = rownames(x = merge@assays$RNA@data),
  value = TRUE)

percent.mito <- Matrix::colSums(merge@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(merge@assays$RNA@counts)

merge <- AddMetaData(
  object = merge,
  metadata = percent.mito,
  col.name = "percent.mito")
VlnPlot(merge, c("nUMI","nGene","percent.mito"), group.by = "sample",pt.size = 0)


merge <- subset(merge, subset = percent.mito < 0.3)
merge <- subset(merge, subset = nGene < 7500)
merge <- subset(merge, subset = nUMI > 1000)
merge <- subset(merge, subset = nUMI < 30000)

merge <- NormalizeData(merge)

merge <- FindVariableFeatures(merge)
VariableFeaturePlot(merge)
merge <- ScaleData(merge, verbose = T)
merge <- RunPCA(merge, npcs = 100, verbose = FALSE)
PCAPlot(merge)
ElbowPlot(merge, ndims = 100)
merge <- RunHarmony(merge, c("sample"))
merge <- FindNeighbors(merge, reduction = "harmony", dims = 1:75, do.plot = T)
merge <- FindClusters(merge, resolution = 1)

merge <- RunUMAP(merge, reduction = "harmony", dims = 1:75)

merge <- CellCycleScoring(merge, g2m.features = cc.genes$g2m.genes,
                          s.features = cc.genes$s.genes)

#########################Extended Data Figure 2-1##############################
celltype <- c(brewer.pal(11,"Set3")[c(1,3:8,10:11)],pal_aaas(alpha = 0.7)(10)[-2],brewer.pal(8,"Set2"),brewer.pal(11,"Paired"))#[13:1]

DimPlot(merge, reduction = "umap", group.by = "celltype2", label = F, repel = TRUE,#pt.size = 1,
        cols = celltype
)


merge <- SetIdent(merge, cells = NULL, 
                  value = merge@meta.data$celltype2)

all.markers_celltype <- FindAllMarkers(merge, verbose = T)

marker <- all.markers_celltype[all.markers_celltype$avg_logFC > 0.8,]
marker <- marker[!duplicated(marker$gene),]

marker %>% group_by(cluster) %>% top_n(60, avg_logFC) -> top_marker

ave <- AverageExpression(merge,return.seurat = T)

cell_exp <- ave@assays$RNA@data[top_marker$gene,]

cell_exp <- (scale(t(cell_exp)))
cell_exp <- na.omit(cell_exp)

#########################Extended Data Figure 2-2##############################
Heatmap(cell_exp, 
        col = colorRamp2(seq(-4,4,length.out = 11),col),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        use_raster=F,
        column_gap = unit(1, "mm"),
        column_title = NULL) 


#####subtyping with consensusclusterplus##################


coding_matrix <- read.csv("./coding_mRNA.csv", sep = ",", header = T, row.names = 1)

coding_matrix <- log2(coding_matrix+1)

tumor <- grep(colnames(d), pattern = "_T", value = T)

d=coding_matrix[top_marker$gene,tumor]

d = sweep(d,1, apply(d,1,median,na.rm=T))
d <- as.matrix(d)

library(ConsensusClusterPlus)
title="./ConsensusClusterPlus/"
results = ConsensusClusterPlus(d,maxK=20,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title=title,
                               clusterAlg="hc",
                               distance="pearson",
                               seed=1262118388.71279,
                               plot="png",
                               verbose = T)
####k=4 was choosen according to the CDF plot#################