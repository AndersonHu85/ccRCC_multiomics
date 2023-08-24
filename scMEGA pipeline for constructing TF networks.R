###############scMEGA pipeline################
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(scMEGA))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(JASPAR2020))
suppressMessages(library(TFBSTools))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))

inputdata.10x <- Read10X_h5("./filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# filter peaks by chromosome
atac_counts <- atac_counts[grep("chr", rownames(atac_counts)), ]

obj.rna <- CreateSeuratObject(counts = rna_counts)
obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")

barcode <- rownames(meta_epi)

obj.rna <- obj.rna[,barcode]

# create seurat object
Epi_aggr <- CreateSeuratObject(
  counts = obj.rna@assays$RNA@counts,
  assay = "RNA",
  meta.data = meta_epi
)

Epi_aggr[["ATAC"]] <- CreateChromatinAssay(
  counts = obj.atac@assays$ATAC@counts,
  sep = c(":", "-"),
    min.cells = 1,
    genome = 'hg38',
    fragments = './atac_fragments.tsv.gz'
)


# extract gene annotations from EnsDb
annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE))

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(Epi_aggr) <- annotations


Epi_aggr@active.assay <- "RNA"


Epi_aggr <- NormalizeData(Epi_aggr)
Epi_aggr <- FindVariableFeatures(Epi_aggr)
Epi_aggr <- ScaleData(Epi_aggr)
Epi_aggr <- RunPCA(Epi_aggr)

Epi_aggr <- RunHarmony(
  Epi_aggr,
  group.by.vars = c("sample"),
  reduction = "pca",
  max.iter.harmony = 30,
  dims.use = 1:30,
  project.dim = FALSE,
  plot_convergence = FALSE
)

epi_sub <- subset(Epi_tra,cells=colnames(Epi_aggr))
Epi_aggr@reductions$harmony <- epi_sub@reductions$harmony

Epi_aggr <- RunDiffusionMap(Epi_aggr, reduction = "harmony",dims = 1:15)

Epi_aggr <- AddTrajectory(object = Epi_aggr, 
                          trajectory = c("IM2_like", "IM4_like"),
                          group.by = "predicted_class", 
                          reduction = "dm",
                          dims = 1:2, 
                          use.all = FALSE)
p <- TrajectoryPlot(object = Epi_aggr, 
                    reduction = "dm",
                    continuousSet = "blueYellow",
                    size = 1,
                    addArrow = FALSE)

p


############SFigure 7F####################
res_TF <- SelectTFs(object = Epi_aggr, return.heatmap = TRUE)
df.cor <- res_TF$tfs
ht_TF <- res_TF$heatmap
draw(ht_TF)



res <- SelectGenes(object = Epi_aggr,
                   labelTop1 = 0,
                   labelTop2 = 0)
df.p2g <- res$p2g

tf.gene.cor <- GetTFGeneCorrelation(object = Epi_aggr, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")


motif.matching <- Epi_aggr@assays$ATAC@motifs@data
colnames(motif.matching) <- Epi_aggr@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
  
df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)
				 
# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
  subset(correlation > 0.5) %>%
  select(c(tf, gene, correlation)) %>%
  rename(weights = correlation)

p <- GRNPlot(df.grn2, 
             tfs.use = c("HNF1B","HNF1A","HNF4G","HNF4A","RORC","TFEC","EMX2","EMX1","AR",
                         "ZBTB7C","BACH1" ,"CEBPB",
                         "FOSL2" ,"FOSL1",
                         "NFKB2","NFKB1","RELA",
                         "EHF","HOXD10","HOXD11","E2F7"
                         ), #####previously defined IM4 or IM2 related TFs
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             #genes.use = c(IM4_feature$IM4_up,IM4_feature$IM2_up),
             genes.highlight = IM4_feature$IM2_up,
             seed = 42, 
             plot.importance = T,
             min.importance = 2,
             remove.isolated = F)

options(repr.plot.height = 40, repr.plot.width = 40)


Epi_aggr <- AddTargetAssay(object = Epi_aggr, df.grn = df.grn2)
p1 <- PseudotimePlot(object = Epi_aggr, tf.use = "HNF1A")+scale_color_npg()
p2 <- PseudotimePlot(object = Epi_aggr, tf.use = "CEBPB")+scale_color_npg()

p1 + p2
