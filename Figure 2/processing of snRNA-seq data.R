library(dplyr)
library(ggplot2)
library(cowplot)
library(harmony)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggsci)
library(Rcpp)
library(RColorBrewer)
library(ggpointdensity)
library(symphony)
library(viridis)
options(stringsAsFactors = F)

snRNA <- readRDS("./snRNA.rds")

VlnPlot(snRNA, c("nUMI","nGene","percent.mito"), group.by = "sample",pt.size = 0)

snRNA <- SetIdent(snRNA, cells = NULL, 
                  value = snRNA@meta.data$celltype2)

snRNA <- NormalizeData(snRNA)

snRNA <- FindVariableFeatures(snRNA, nfeatures = 2000)
VariableFeaturePlot(snRNA)
snRNA <- ScaleData(snRNA, verbose = T)
snRNA <- RunPCA(snRNA, npcs = 50, verbose = FALSE)
PCAPlot(snRNA,cols=celltype)
ElbowPlot(snRNA, ndims = 50)
snRNA <- RunHarmony(snRNA, c("sample"),max.iter.harmony = 20)
snRNA <- FindNeighbors(snRNA, reduction = "harmony", dims = 1:50, do.plot = T)
snRNA <- FindClusters(snRNA, resolution = 1)

snRNA  <- RunUMAP(snRNA, reduction = "harmony", dims = 1:50, n.components = 2)

snRNA <- CellCycleScoring(snRNA, g2m.features = cc.genes$g2m.genes,
                        s.features = cc.genes$s.genes)


########################### Figure 2E #############################

DimPlot(snRNA, reduction = "umap",  label = F, repel = TRUE, raster= F,group.by = "celltype",
        cols = c("#fc772d", "#c6d5a8", "#b48f60", "#a75565", "#e7553e")
)
DimPlot(snRNA, reduction = "umap",  label = F, repel = TRUE, raster= F,group.by = "class",
        cols = brewer.pal(11,"Set3")[c(1,3:5)]
)

#######Sub-clustering of different cell lineages were performed using the same pipeline


########################### Figure 2F #############################

column <- as.data.frame(table(snRNA$sample, snRNA$celltype))
group <- as.data.frame(table(snRNA$sample, snRNA$class))
group <- group[group$Freq>0,]

column$IM_subtype <- plyr::mapvalues(x = column$Var1, 
                                 from = as.character(group$Var1), 
                                 to = as.character(group$Var2))

colnames(column) <- c("Patient","cluster","Freq","IM_subtype")

column$IM_subtype <- as.character(column$IM_subtype)

p1 <- ggplot(column,aes(Patient,weight=Freq,fill=cluster))+geom_bar(position="fill")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
  scale_fill_manual(values = c(celltype))
p1 +  facet_grid(~IM_subtype,scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(.2, "cm")) + ylab("Proportion")



