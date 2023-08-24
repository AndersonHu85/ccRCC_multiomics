library(ArchR)
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
library(parallel)
options(stringsAsFactors = F)

addArchRThreads(threads = 35) 
memory.limit(1000000000)
addArchRGenome("hg38")

inputFiles <- "./atac_fragments.tsv.gz" 
names(inputFiles) <- "aggr"


############load ATAC data################
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ccRCC",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
df

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  size =0.5,
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)


##############select single nucleic with paired snRNA-seq data##############

id <- projHeme1$cellNames
id <- gsub("aggr#","",id)
id2 <- colnames(annotated_snRNA)[which(annotated_snRNA$cohort=="10X_ARC")]
idxPass <- which(id %in% id2)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme2 <- projHeme1[cellsPass, ]

meta <- annotated_snRNA@meta.data[cellNames,]
rownames(meta) <- paste("aggr#",rownames(meta),sep = "")
meta <- meta[projHeme2$cellNames,]

projHeme2 <- addCellColData(ArchRProj = projHeme2, data = meta$class,
                            cells = cellsPass, name = "class",force = TRUE)
projHeme2 <- addCellColData(ArchRProj = projHeme2, data = meta$sample,
                            cells = cellsPass, name = "sample",force = TRUE)
projHeme2 <- addCellColData(ArchRProj = projHeme2, data = meta$celltype,
                            cells = cellsPass, name = "celltype",force = TRUE)


########dimensional reduction###############
projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.8), 
    sampleCells = 20000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

projHeme2 <- addUMAP(
  ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "UMAP", 
  nNeighbors = 30, minDist = 0.5, metric = "cosine"
)

######add snRNA-seq data##############
seRNA <- import10xFeatureMatrix(
  input = c("./filtered_feature_bc_matrix.h5"),
  names = c("aggr")
)

seRNA <- subset(seRNA, cells=projHeme2@cellColData@rownames)

projHeme2 <- addGeneExpressionMatrix(
  input = projHeme2,
  seRNA = seRNA,
  chromSizes = getChromSizes(projHeme2),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)




meta_EC <- snRNA_EC@meta.data ####a subset of snRNA-seq data with only endothelial cells
rownames(meta_EC) <- paste("aggr#",rownames(meta_EC),sep = "")
meta_EC <- meta_EC[rownames(meta_EC)%in%projHeme2$cellNames,]

EC <- projHeme2[rownames(meta_EC), ]


EC <- addUMAP(
  ArchRProj = EC, reducedDims = "IterativeLSI", name = "UMAP", force = TRUE,
  nNeighbors = 30, minDist = 0.5, metric = "cosine"
)

EC <- addImputeWeights(EC)


##########SFigure 7c##########
pathToMacs2 <- ("/home/usr/miniconda3/bin/macs2")

EC <- addReproduciblePeakSet(
  ArchRProj = EC, 
  groupBy = "class", 
  pathToMacs2 = pathToMacs2
)

EC <- addPeakMatrix(EC)

EC <- addPeak2GeneLinks(
  ArchRProj = EC,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneExpressionMatrix"
)

################## Extended Data Figure 4D ####################

p <- plotBrowserTrack(
  ArchRProj = EC, 
  groupBy = "class", 
  geneSymbol = "FOS", 
  upstream = 40000,
  downstream = 50000,
  loops = getPeak2GeneLinks(EC)
)
grid::grid.newpage()
grid::grid.draw(p[[1]])

p <- plotBrowserTrack(
  ArchRProj = EC, 
  groupBy = "class", 
  geneSymbol = "FOS", 
  upstream = 40000,
  downstream = 50000,
  loops = getPeak2GeneLinks(EC)
)
grid::grid.newpage()
grid::grid.draw(p[[1]])

