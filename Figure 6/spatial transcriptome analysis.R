############ ST data was processed per sample, use X98_T as the example ##############

library(SPATA2)

X98 <-
  initiateSpataObject_10X(
    directory_10X = 'dir/Y7/filtered_feature_bc_matrix/', # the directory from which to load the data. It has been deposited at Mendeley and Zenodo
    sample_name = "X98"
  )

plotSurface(object = X98,
            color_by = c("SFTPA2"), display_image=T,#alpha_by = "SCGB3A1",
            smooth = F, normalize =F,#smooth_span = 0.5,
            pt_size = 1.2) + 
  ggplot2::labs(title = "Denoised")

################ Figure 6K-1 ########################
plotSurface(object = X98, color_by = "seurat_clusters", pt_clrp = "npg")


ref_cor <- Read10X("/R_cor/filtered_feature_bc_matrix/") ####use ST data of renal cortex as the control
ref_cor <- as.matrix(as.matrix(ref_cor))
colnames(ref_cor) <- paste(colnames(ref_cor),"ref",sep = "_")

meta_ref <- data.frame(rep("ref",ncol(ref_cor)),row.names = colnames(ref_cor))
colnames(meta_ref) <- "sample"

gtf1 <- rtracklayer::import('./Homo_sapiens.GRCh38.87.chr.gtf')
gtf_df <- as.data.frame(gtf1)
gtf_df <- gtf_df[,c(10,12,1:3)]
gtf_df <- gtf_df[!duplicated(gtf_df[,2]),]
colnames(gtf_df) <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                      "start_position" , "end_position") ###make sure colnames are in this order
gene <- rownames(ref_cor)
gtf_df <- gtf_df[gtf_df$hgnc_symbol %in% gene,]
gtf_df <- gtf_df[gtf_df$chromosome_name != "MT",]
gtf_df$chromosome_name <- as.character(gtf_df$chromosome_name)


X98 <- runCnvAnalysis(
  object = X98,
  directory_cnv_folder = "./X98_T/CNV/", # example directory to 
  ref_annotation = meta_ref, 
  ref_mtr = ref_cor, 
  ref_regions = cnv_ref$regions,
  gene_pos_df = gtf_df
)

################# Figure 6K-2 ##################
infercnv_obj <- readRDS("./X98_T/CNV/run.final.infercnv_obj")
tumor_mat <- infercnv_obj@expr.data
tumor_mat <- as.matrix(tumor_mat)
gene_order <- infercnv_obj@gene_order
gene_order <- as.data.frame(gene_order)
gene_order <- gene_order[gene_order[,1]%in%c(1:22),]
gene_order$chr <- as.numeric(as.character(gene_order$chr))
gene_order <- gene_order[order(gene_order$chr,decreasing = F),]
tumor_mat <- tumor_mat[rownames(gene_order),]
tumor_mat <- t(tumor_mat)
tumor_mat[1:5,1:5]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(ggplot2)

seurat_clusters <- pal_npg()(8)
names(seurat_clusters) <- c(0,1,2,3,4,5,6,7)
cl_col =list(seurat_clusters=seurat_clusters
)

right_anno <- rowAnnotation(df = annotation,col = cl_col,show_legend = T)
col <- c(brewer.pal(11, "RdBu"))[c(11:1)]

heat <- Heatmap(tumor_mat[X98_barcode,], #X98_barcode is the barcode of the SPATA2 object
                col = colorRamp2(seq(0.7,1.3,length.out = 11),col),
                cluster_rows = T,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = F,
                use_raster=T,
                border = "BLACK",
                show_row_dend = T,
                clustering_method_rows = "ward.D",
                row_gap = unit(0.5, "mm"),
                column_split = gene_order$chr,
                column_gap = unit(0, "mm"),
                row_split = annotation$seurat_clusters,
                row_gap = unit(0, "mm"),
                top_annotation = top_anno,
                right_annotation = right_anno,
                column_title = NULL)
heat

################# Figure 6L #####################
plotSurface(object = X98, color_by = "Chr9", pt_clrsp = "Blues 3", rev = F, c1 = 1,smooth = T)


X98 <-
  runAutoencoderDenoising(
    object = X98, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

X98 <- setActiveExpressionMatrix(object = X98, mtr_name = "denoised")

IM4_feature <- read.table("./IM4_feature.txt",header = T) ######part of the tmp_org file used in bulk RNA-seq analysis
pathway <- IM4_feature
pathway_list <- vector("list",length(pathway))

for (i in seq_along(pathway)) {
  pathway_list[[i]] <- unique(na.omit(pathway[,i]))
}

names(pathway_list) <- colnames(pathway)


pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) 
})

################# Figure 6J ####################
plotSurfaceAverage(object = X98, 
                   color_by = c(pathway_list),
                   smooth = F, 
                   pt_size = 2) + 
  ggplot2::labs(title = "Denoised")

################# Extended Data Figure 10A ##############
X98 <- createTrajectories(object = X98)

plotTrajectory(object = X98, 
               trajectory_name = "DCCD",
               color_by = "seurat_clusters",
               pt_clrp = "npg",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               display_image = FALSE) +
  legendTop()

################# Extended Data Figure 10B ##############
TF <- read.csv("./TF_final.csv") ####the final version of TF list

hm_colors <- viridis::inferno(n = 256)

tra_spata <- SPATA2::plotTrajectoryHeatmap(object = X98, 
                                           trajectory_name = "DCCD",
                                           variables = TF$symbol, 
                                           arrange_rows = "maxima",
                                           colors = hm_colors,
                                           cluster_rows = T,
                                           show_rownames = T, 
                                           split_columns = F, # splits the heatmap to highlight the trajectory parts
                                           smooth_span = 0.5)


