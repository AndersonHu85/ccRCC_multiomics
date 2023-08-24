########### following the process of Figure 2

meta_epi <- Epi_snRNA@meta.data
rownames(meta_epi) <- paste("aggr#",rownames(meta_epi),sep = "")
meta_epi <- meta_epi[rownames(meta_epi)%in%Epi$cellNames,]
length(Epi$cellNames)
Epi <- projHeme2[rownames(meta_epi), ]

Epi <- addCellColData(ArchRProj = Epi, data = meta_epi$predicted_class,
                      cells = rownames(meta_epi), name = "predicted_class",force = TRUE)


Epi <- addUMAP(
  ArchRProj = Epi, reducedDims = "IterativeLSI", name = "UMAP", force = TRUE,
  nNeighbors = 30, minDist = 0.5, metric = "cosine"
)

Epi <- addImputeWeights(Epi)


pathToMacs2 <- ("dir/macs2") #path to MACS2

Epi <- addReproduciblePeakSet(
  ArchRProj = Epi, 
  groupBy = "sample", 
  pathToMacs2 = pathToMacs2
)

Epi <- addPeakMatrix(Epi)

Epi <- addPeak2GeneLinks(
  ArchRProj = Epi,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneExpressionMatrix"
)

############ Extended Data Figure 10C ################
p <- plotPeak2GeneHeatmap(ArchRProj = Epi, groupBy = "sample")
p
plotPDF(plotList = list(p), 
        name = "P2L_heat.pdf", 
        ArchRProj = Epi, 
        addDOC = FALSE, width = 9, height = 6)


############Figure 7L#############
Epi <- addMotifAnnotations(ArchRProj = Epi, motifSet = "cisbp", name = "Motif")

motifPositions <- getPositions(Epi)
motifs <- c("AR")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

seFoot <- getFootprints(
  ArchRProj = Epi, 
  positions = motifPositions["AR_689"], 
  groupBy = "sample"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = Epi, 
  normMethod = "Subtract",
  plotName = "AR_motif",
  addDOC = FALSE,
  plot = T,
  smoothWindow = 5
)

#############Figure 7N#############
p <- plotBrowserTrack(
  ArchRProj = Epi, 
  groupBy = "sample", 
  geneSymbol = "CEBPB", 
  upstream = 100000,
  downstream = 50000,
  loops = getPeak2GeneLinks(Epi)
)
grid::grid.newpage()
grid::grid.draw(p[[1]])


plotPDF(plotList = p, 
        name = "Plot-Tracks-CEBPB.pdf", 
        ArchRProj = Epi, 
        addDOC = FALSE, width = 9, height = 6)


############Figure 7H#############

Epi <- addReproduciblePeakSet(
  ArchRProj = Epi, 
  groupBy = "predicted_class", 
  pathToMacs2 = pathToMacs2
)

Epi <- addPeakMatrix(Epi)

Epi <- addPeak2GeneLinks(
  ArchRProj = Epi,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneExpressionMatrix"
)

Epi <- addMotifAnnotations(ArchRProj = Epi, motifSet = "cisbp", name = "Motif")

Epi <- addBgdPeaks(Epi,force = TRUE)

Epi <- addDeviationsMatrix(
  ArchRProj = Epi, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(Epi, name = "MotifMatrix", plot = TRUE,n = 15)

plotVarDev


#############Figure 7I#############
motifPositions <- getPositions(Epi)
motifs <- c("HNF1A","HNF4A","FOSL1","JUNB")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

seFoot <- getFootprints(
  ArchRProj = Epi, 
  positions = motifPositions[markerMotifs], 
  groupBy = "predicted_class"
)

col <- c("IM2_like"="#4DBBD5FF","IM4_like"="#E64B35FF")

plotFootprints(
  seFoot = seFoot,
  ArchRProj = Epi, 
  pal = col,
  normMethod = "Subtract",
  plotName = "HNF_AP1",
  addDOC = FALSE,
  plot = T,
  #height = 6.5,
  smoothWindow = 5
)