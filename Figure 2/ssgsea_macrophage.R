Myeloid <- snRNA[,snRNA$celltype%in%c("Myeloid")]

matrix <- Myeloid@assays$RNA@data[,grep("Macro",Myeloid$celltype2)]
res_h <- gsva(as.matrix(matrix), geneset, method = "ssgsea", parallel.sz = 0, verbose = T)
res2 <- gsva(as.matrix(matrix), geneset2, method = "ssgsea", parallel.sz = 0, verbose = T)
res5 <- gsva(as.matrix(matrix), geneset5, method = "ssgsea", parallel.sz = 0, verbose = T)

a <- unique(Myeloid$celltype2[grep("Macro",Myeloid$celltype2)])

results <- lapply(a, function(m) {
  group <- ifelse(Myeloid$celltype2 == m, "select", "ref")
  group_mat <- model.matrix(~factor(group))
  colnames(group_mat) <- c("ref", "select")
  fit <- lmFit(rbind(res_h, res2, res5), group_mat)
  fit2 <- eBayes(fit)
  allDiff <- topTable(fit2, adjust = 'fdr', coef = 2, number = 200000)
  allDiff$X <- row.names(allDiff)
  ##### only select T-value from the result ############
  allDiff <- allDiff[, c("X", "t")]
  colnames(allDiff)[2] <- m
  return(allDiff)
})

b <- Reduce(function(x, y) {merge(x, y, by = "X")}, results)
b <- data.frame(b[-1],row.names = b$X)

ssgsea_macro <- c("GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY", "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX",                      
                  "GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE", "GOBP_INTERLEUKIN_18_PRODUCTION",                         
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GOMF_SCAVENGER_RECEPTOR_ACTIVITY",                       
                  "GOBP_ENDOTHELIAL_CELL_FATE_COMMITMENT", "GOCC_GROWTH_FACTOR_COMPLEX" )

ssgsea_macro <- b[ssgsea_macro,]
rownames(ssgsea_macro) <- tolower(rownames(ssgsea_macro))
rownames(ssgsea_macro) <-  strsplit(rownames(ssgsea_macro),"_")[] %>% 
  lapply(., function(x) x[-1]) %>% 
  sapply(., function(x) paste(x, collapse = " "))

####################### Extended Data Figure 3J ##############################  

bk=c(seq(-20,-0.1,by=0.01),seq(0,20,by=0.01))
pheatmap(ssgsea_macro, border_color = NA, fontsize = 9,cellheight = 12,cellwidth = 15,
         cluster_col=T, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = brewer.pal(11,"RdBu")[11:6])(length(bk)/2),
                   colorRampPalette(colors = brewer.pal(11,"RdBu")[6:1])(length(bk)/2)))
