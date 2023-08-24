library(KEGGREST)
library(plyr)


LC_matrix <- fread("./LC_matrix.csv") 
GC_matrix <- fread("./GC_matrix.csv")
anno_meta <- list(LC_MS=LC_matrix[,1:15],
                  GC_MS=GC_matrix[,1:10])
LC_matrix <- data.frame(LC_matrix[,16:165],row.names = LC_matrix$Metabolites)
GC_matrix <- data.frame(GC_matrix[,11:160],row.names = GC_matrix$`Metabolite name`)


matrix_meta <- rbind(LC_matrix,GC_matrix)
matrix_meta <- log2(matrix_meta)

mean_N <- rowMeans(matrix_meta[,1:50])
mean_T <- rowMeans(matrix_meta[,51:150])
diff_meta <- data.frame(mean_N,mean_T)
diff_meta$log2FC <- diff_meta$mean_T-diff_meta$mean_N
p_value <- c()
for (i in 1:nrow(matrix_meta)) {
  p <- wilcox.test(t(matrix_meta[i,1:50]),t(matrix_meta[i,51:150]))
  p_value <- c(p_value,p[["p.value"]])
}
diff_meta$p_value <- p_value
diff_meta$adj_p <- p.adjust(p_value,"fdr")
diff_meta$metabollites <- rownames(diff_meta)
diff_meta$group <- c(rep("LC",nrow(LC_matrix)),rep("GC",nrow(GC_matrix)))
diff_meta$KEGG_ID <- c(anno_meta$LC_MS$kegg,anno_meta$GC_MS$KEGG)
diff_meta$HMDB <- c(anno_meta$LC_MS$`Compound ID`,anno_meta$GC_MS$HMDB)


##### load kegg annotation file #############
# source("./get.kegg.all.R")
# source("./get.kegg.byId.R")
# 
# keggAll = get.kegg.all()
# save(keggAll, file = "keggAll.Rdata")

load("keggAll.Rdata")

KEGG_anno <- keggAll[["compound"]]
KEGG_anno <- KEGG_anno[KEGG_anno$ENTRY %in%  c(anno_meta$LC_MS$kegg,anno_meta$GC_MS$KEGG),]
KEGG_anno <- KEGG_anno[,c(1,8)]


test <- KEGG_anno[, 2]
test <- strsplit(test, "///")
test <- unlist(test)
path_total <- unique(test)

path <- matrix(NA, nrow = length(path_total), ncol = nrow(KEGG_anno))
colnames(path) <- KEGG_anno[, 1]

path_list <- vector("list", nrow(KEGG_anno))

for (i in 1:nrow(KEGG_anno)) {
  test <- strsplit(KEGG_anno[i, 2], "///")[[1]]
  test <- as.data.frame(x = test, row.names = test)[path_total, , drop = FALSE]
  test <- as.matrix(test)
  colnames(test) <- KEGG_anno[i, 1]
  path_list[[i]] <- test
}

path <- do.call(cbind, path_list)

b <- rowSums(!is.na(path))
path <- data.frame(b,path)

############# extract differetilly expressed metabolites ###############

up <- diff_meta$KEGG_ID[diff_meta$log2FC > 1 & diff_meta$adj_p < 0.05]
up <- up[up %in% colnames(path)] %>% unique()
down <- diff_meta$KEGG_ID[diff_meta$log2FC < -1 & diff_meta$adj_p < 0.05]
down <- down[down %in% colnames(path)] %>% unique()


############# calculate number of differetilly expressed metabolites in each pathway ###############

up <- as.matrix(path[, up])
b <- apply(up, 1, function(x) sum(!is.na(x)))
path$up <- b

down <- as.matrix(path[, down])
b <- apply(down, 1, function(x) sum(!is.na(x)))
path$down <- b

path <- path[,c(ncol(path),ncol(path)-1,1:(ncol(path)-2))]
colnames(path)[1:3] <- c("down","up","total")
path$DA <- (path$up-path$down)/path$total
path <- path[,c(ncol(path),1:(ncol(path)-1))]
rownames(path) <- path_total
path <- path[path$total >= 3,]

PATH_use <- c("map01210: 2-Oxocarboxylic acid metabolism",             "map02010: ABC transporters" ,                          
              "map00250: Alanine, aspartate and glutamate metabolism", "map00592: alpha-Linolenic acid metabolism" ,           
              "map00520: Amino sugar and nucleotide sugar metabolism", "map00970: Aminoacyl-tRNA biosynthesis" ,               
              "map00590: Arachidonic acid metabolism",                 "map00330: Arginine and proline metabolism" ,           
              "map00220: Arginine biosynthesis",                       "map00053: Ascorbate and aldarate metabolism" ,         
              "map00410: beta-Alanine metabolism",                     "map01230: Biosynthesis of amino acids",                
              "map04973: Carbohydrate digestion and absorption",       "map04979: Cholesterol metabolism",                     
              "map00020: Citrate cycle (TCA cycle)" ,                  "map00460: Cyanoamino acid metabolism"  ,               
              "map00270: Cysteine and methionine metabolism",          "map00472: D-Arginine and D-ornithine metabolism" ,     
              "map00471: D-Glutamine and D-glutamate metabolism",      "map00051: Fructose and mannose metabolism" ,           
              "map00052: Galactose metabolism",                        "map00480: Glutathione metabolism" ,                    
              "map00561: Glycerolipid metabolism" ,                    "map00564: Glycerophospholipid metabolism",             
              "map00260: Glycine, serine and threonine metabolism" ,   "map00010: Glycolysis / Gluconeogenesis" ,              
              "map00300: Lysine biosynthesis" ,                        "map00310: Lysine degradation",                         
              "map00760: Nicotinate and nicotinamide metabolism",      "map00190: Oxidative phosphorylation" ,                 
              "map00770: Pantothenate and CoA biosynthesis" ,          "map00030: Pentose phosphate pathway",                  
              "map00360: Phenylalanine metabolism",                    "map03320: PPAR signaling pathway",                     
              "map00640: Propanoate metabolism",                       "map04974: Protein digestion and absorption" ,          
              "map00230: Purine metabolism"   ,                        "map00240: Pyrimidine metabolism" ,                     
              "map00620: Pyruvate metabolism" ,                        "map04923: Regulation of lipolysis in adipocytes",      
              "map00430: Taurine and hypotaurine metabolism",          "map00380: Tryptophan metabolism",                      
              "map00350: Tyrosine metabolism" ,                        "map00290: Valine, leucine and isoleucine biosynthesis",
              "map00280: Valine, leucine and isoleucine degradation",  "map04270: Vascular smooth muscle contraction")

DA2 <- path[PATH_use,1:4]
DA2$color = ifelse(DA2$DA == 0, "grey75",
                   ifelse(DA2$DA > 0, "red",
                          ifelse(DA2$DA < -0,"blue","grey75")))
DA2$X <- sapply(strsplit(rownames(DA2),": "),"[",2)

library(ggplot2)
library(ggsci)
library(RColorBrewer)

DA2$total[DA2$total>20]=20 ###balance the size of each dot
DA2 <- DA2[order(DA2$X,decreasing = T),]
DA2$X <- factor(DA2$X, levels = DA2$X)

col <- brewer.pal(11,"RdBu")[11:1]

#####################Figure 1I#########################

ggplot(DA2,aes(DA,X),color = DA, fill = color)+
  geom_bar(aes(DA,X),data=DA2,position=position_dodge(),width=0.15, stat="identity") +
  geom_point(aes(DA,X, size=total, color = DA),data = DA2)+
  theme_bw() + scale_color_gradientn(colors = col,limits=c(-0.85,0.85)) + 
  xlim(-0.85,0.85) + geom_vline(xintercept = 0.2,linetype=2) + 
  geom_vline(xintercept = -0.2,linetype=2) + 
  labs(x = "DA score", y = "Pathways")+ 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
ggsave("DA_sup_IM1_3.pdf",width = 10,height = 3)


