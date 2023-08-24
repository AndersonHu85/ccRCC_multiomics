library(dplyr)
library(ggplot2)
library(cowplot)
library(harmony)
library(Seurat)
library(ggsci)
library(Rcpp)
library(RColorBrewer)
options(stringsAsFactors = F)

snRNA_NAT <- readRDS("./normal_snRNA.rds")

snRNA_NAT@meta.data$nUMI <- snRNA_NAT@meta.data$nCount_RNA
snRNA_NAT@meta.data$nGene <- snRNA_NAT@meta.data$nFeature_RNA


mito.genes <- grep(
  pattern = "^MT-",
  x = rownames(x = snRNA_NAT@assays$RNA@data),
  value = TRUE)

percent.mito <- Matrix::colSums(snRNA_NAT@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(snRNA_NAT@assays$RNA@counts)

snRNA_NAT <- AddMetaData(
  object = snRNA_NAT,
  metadata = percent.mito,
  col.name = "percent.mito")
VlnPlot(snRNA_NAT, c("nUMI","nGene","percent.mito"), group.by = "sample",pt.size = 0)

snRNA_NAT <- subset(snRNA_NAT, subset = percent.mito < 0.05)
snRNA_NAT <- subset(snRNA_NAT, subset = n_genes > 500)
snRNA_NAT <- subset(snRNA_NAT, subset = nUMI > 1000)

snRNA_NAT <- snRNA_NATizeData(snRNA_NAT)

snRNA_NAT <- FindVariableFeatures(snRNA_NAT, nfeatures = 3000)
VariableFeaturePlot(snRNA_NAT)
snRNA_NAT <- ScaleData(snRNA_NAT, verbose = T)
snRNA_NAT <- RunPCA(snRNA_NAT, npcs = 50, verbose = FALSE)
PCAPlot(snRNA_NAT)
ElbowPlot(snRNA_NAT, ndims = 50)
snRNA_NAT <- RunHarmony(snRNA_NAT, c("sample"))
snRNA_NAT <- FindNeighbors(snRNA_NAT, reduction = "harmony", dims = 1:50, do.plot = T)
snRNA_NAT <- FindClusters(snRNA_NAT, resolution = 1)

snRNA_NAT <- RunUMAP(snRNA_NAT, reduction = "harmony", dims = 1:50, n.components = 2)

p1 <- DimPlot(snRNA_NAT, reduction = "umap", group.by = "sample",cols = celltype, label = TRUE)+ NoLegend()
p2 <- DimPlot(snRNA_NAT, reduction = "umap", group.by = "RNA_snn_res.1", cols = celltype, pt.size = 0.4,
              label = TRUE,raster = F) + NoLegend()
plot_grid(p1,p2)

snRNA_NAT <- CellCycleScoring(snRNA_NAT, g2m.features = cc.genes$g2m.genes,
                           s.features = cc.genes$s.genes)

celltype <- c(brewer.pal(11,"Set3")[c(1,3:8,10:11)],pal_npg(alpha = 0.7)(9),
              brewer.pal(11,"Paired")[11:1],brewer.pal(8,"Set1")[8:1],
              brewer.pal(8,"Set2"))#[13:1]

############ Extended Data Figure 6B ##################
DimPlot(snRNA_NAT, reduction = "umap", group.by = "celltype2", label = T, repel = TRUE,
        cols = celltype, pt.size =0.4,
)

############# scoring 4 modules in snRNA-seq data ##################

genelist <- list(module1=genelist_Figure4a[a[[1]]],
                 module2=genelist_Figure4a[a[[2]]],
                 module3=genelist_Figure4a[a[[3]]],
                 module4=genelist_Figure4a[a[[4]]],)#a is the row cluster tree of Figure 4a

snRNA_NAT <- AddModuleScore(snRNA_NAT,genelist)

############ Extended Data Figure 6D ##################
FeaturePlot(snRNA_NAT, c("Cluster1","Cluster2","Cluster3","Cluster4"),
            raster=F,max.cutoff = "q99", min.cutoff = "q0",
            cols = c("grey85",brewer.pal(9,"YlOrRd")[-1]),reduction = "umap", slot = "data")


library(COSG)

snRNA_NAT <- SetIdent(snRNA_NAT, cells = NULL, 
                  value = snRNA_NAT@meta.data$celltype2)

marker_cosg <- cosg(
  snRNA_NAT,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=5)

head(marker_cosg$names)
marker=marker_cosg$names
marker <- marker[,order(colnames(marker))]

############ Extended Data Figure 6C ##################

col <- viridis::viridis(11,option = "A")
A <- DotPlot(object = snRNA_NAT, features = unique(as.character(as.matrix(marker))), 
             group.by = "celltype3")+scale_color_gradientn(colours = col)
A + ggthemes::theme_base() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
