library(ggplot2)
library(ggsci)


coding_combat2_mat <- read.csv("./bulkRNA_matrix_TPM.csv",row.names = 1) ###matrix of bulk RNA-seq
coding_combat2_mat[1:5,1:5]

pca <- prcomp(t(coding_combat2_mat), scale=F) 
pca <- pca[["x"]]
pca <- as.data.frame(pca)


#############################Figure 1G_1###################################

NT <- sapply(strsplit(rownames(pca),"_"),"[",2)
group <- ArchR::ArchRPalettes[2][[1]][c(14,16)]
names(group) <- c("N","T")

ggplot(pca, aes(PC1,PC2, color = NT)) + geom_point(size=2,alpha=1) + 
  scale_color_manual(values = group)+
  stat_ellipse(data=pca,aes(x=PC1,y=PC2,fill=NT,color=NT),
               geom = "polygon",alpha=0.2,level=0.96,type="t",linetype = 0,show.legend = F)+
  scale_fill_manual(values = group)+theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15, face = 'plain'),
    axis.title = element_text(color = 'black',size = 15, face = 'plain'),
    axis.ticks = element_line(color = 'black')) + #coord_fixed(ratio = 1)+
  geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")


################processing global proteomics data###############

library(data.table)
library(NormalyzerDE)

protein <- fread("../../蛋白组/protein_gene_level.csv")
protein <- data.frame(protein[,-1], row.names = protein$V1)

NT <- sapply(strsplit(colnames(protein),"_"),"[",2)

design <- data.frame(sample=colnames(protein),group=NT)

write.table(design,"design_protein.txt",col.names = T, row.names = F,sep = "\t",quote = F)
write.table(protein,"matrix_protein.txt",col.names = T, row.names = F,sep = "\t",quote = F)


jobName <- "vignette_run"

experimentObj <- setupRawDataObject("matrix_protein.txt", "design_protein.txt", 
                                    "default", TRUE, "sample", "group")
normObj <- getVerifiedNormalyzerObject(jobName, experimentObj)
normResults <- normMethods(normObj)

protein_VSN <- normResults@normalizations$VSN
protein_VSN[is.na(protein_VSN)==TRUE]=0
rownames(protein_VSN) <- rownames(protein)

pca <- prcomp(t(protein_VSN), scale=F) 
pca <- pca[["x"]]
pca <- as.data.frame(pca)


#############################Figure 1G_2###################################

NT <- sapply(strsplit(rownames(pca),"_"),"[",2)
group <- ArchR::ArchRPalettes[2][[1]][c(14,16)]
names(group) <- c("N","T")

ggplot(pca, aes(PC1,PC2, color = NT)) + geom_point(size=2,alpha=1) + 
  scale_color_manual(values = group)+
  stat_ellipse(data=pca,aes(x=PC1,y=PC2,fill=NT,color=NT),
               geom = "polygon",alpha=0.2,level=0.96,type="t",linetype = 0,show.legend = F)+
  scale_fill_manual(values = group)+theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15, face = 'plain'),
    axis.title = element_text(color = 'black',size = 15, face = 'plain'),
    axis.ticks = element_line(color = 'black')) + 
  geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")





################processing metabolism data#################

library(data.table)
library(NormalyzerDE)
library(dplyr)

LC_matrix <- fread("./LC_matrix.csv") 
LC_matrix <- data.frame(LC_matrix[,16:165],row.names = LC_matrix$Metabolites)

GC_matrix <- fread("./GC_matrix.csv")
GC_matrix <- data.frame(GC_matrix[,11:160],row.names = GC_matrix$`Metabolite name`)

all(colnames(LC_matrix)==colnames(GC_matrix))

matrix_meta <- rbind(LC_matrix,GC_matrix)
matrix_meta <- log2(matrix_meta)

NT <- sapply(strsplit(colnames(matrix_meta),"_"),"[",2)


pca <- prcomp(t(matrix_meta), scale=F) 
pca <- pca[["x"]]
pca <- as.data.frame(pca)


#############################Figure 1G_2###################################

NT <- sapply(strsplit(rownames(pca),"_"),"[",2)
group <- ArchR::ArchRPalettes[2][[1]][c(14,16)]
names(group) <- c("N","T")

ggplot(pca, aes(PC1,PC2, color = NT)) + geom_point(size=2,alpha=1) + 
  scale_color_manual(values = group)+
  stat_ellipse(data=pca,aes(x=PC1,y=PC2,fill=NT,color=NT),
               geom = "polygon",alpha=0.2,level=0.96,type="t",linetype = 0,show.legend = F)+
  scale_fill_manual(values = group)+theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15, face = 'plain'),
    axis.title = element_text(color = 'black',size = 15, face = 'plain'),
    axis.ticks = element_line(color = 'black')) + #coord_fixed(ratio = 1)+
  geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")



############correlation between CNV, mRNA and protein##################

cnv_gene <- read.table("../../SNV/gistic/broad_data_by_genes.txt",sep = "\t",header = T, row.names = 1)

cnv_gene <- cnv_gene[,grep("_T",colnames(cnv_gene),invert = T)]
meta_loc <- cnv_gene[,1:2]
cnv_gene <- cnv_gene[,3:97]
colnames(cnv_gene) <- gsub(".sorted.rmdup","",colnames(cnv_gene))
colnames(cnv_gene) <- gsub("Tp","T",colnames(cnv_gene))
colnames(cnv_gene) <- gsub("TP","T",colnames(cnv_gene))
colnames(cnv_gene)[colnames(cnv_gene)=="w62T"]="W62T"
colnames(cnv_gene) <- sub("(.)$", "_\\1", colnames(cnv_gene))

library(multiOmicsViz)

##################Extended data figure 1G-1##################
gene_select <- intersect(rownames(cnv_gene),rownames(coding_combat2_mat))
WTS_cnv <- coding_combat2_mat[gene_select,colnames(cnv_gene)]

targetOmicsList <- list()
targetOmicsList[[1]] <- WTS_cnv

chr <- as.character(c(1:22))
outputfile <- paste("./","/heatmap_mRNA_CNV",sep="")
multiOmicsViz(cnv_gene[,],sourceOmicsName="CNV","All",targetOmicsList,
              "mRNA","All",0.01,outputfile=outputfile)

##################Extended data figure 1G-2##################
gene_select <- intersect(rownames(cnv_gene),rownames(protein_VSN))
protein_cnv <- protein_VSN[gene_select,colnames(cnv_gene)]

targetOmicsList <- list()
targetOmicsList[[1]] <- protein_cnv

outputfile <- paste("./","/heatmap_protein_CNV",sep="")
multiOmicsViz(cnv_gene, sourceOmicsName="CNV","All",targetOmicsList,
              "protein","All",0.01,outputfile=outputfile)



##############ssgsea analysis##################
########## obtained from MSigDB ###############
library(GSVA)
library(GSEABase)
library(limma)

geneset <- getGmt("./h.all.v7.2.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category="h"))
geneset2 <- getGmt("./c2.all.v7.5.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category="c2"))
geneset5 <- getGmt("./c5.all.v7.5.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category="c2"))

matrix <- protein_VSN
res_h_protein <- gsva(as.matrix(matrix), geneset, method = "ssgsea", parallel.sz = 0, verbose = T)
res2_protein <- gsva(as.matrix(matrix), geneset2, method = "ssgsea", parallel.sz = 0, verbose = T)
res5_protein <- gsva(as.matrix(matrix), geneset5, method = "ssgsea", parallel.sz = 0, verbose = T)

NT <- sapply(strsplit(colnames(matrix),"_"),"[",2)

################limma##################

group<-model.matrix(~factor(NT))
colnames(group) <- c("normal","tumor")
fit<-lmFit(rbind(res_h,res2,res5),group)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
allDiff$X <- row.names(allDiff)
allDiff$X <- tolower(allDiff$X)

gsva_use <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",  "HALLMARK_HYPOXIA",                                
              "HALLMARK_INFLAMMATORY_RESPONSE",   "HALLMARK_GLYCOLYSIS",                            
              "HALLMARK_MTORC1_SIGNALING",  "HALLMARK_ANGIOGENESIS",                         
              "HALLMARK_ADIPOGENESIS",  "HALLMARK_FATTY_ACID_METABOLISM",              
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION",   "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION", 
              "KEGG_ARGININE_AND_PROLINE_METABOLISM",    "KEGG_TAURINE_AND_HYPOTAURINE_METABOLISM",         
              "KEGG_RETINOL_METABOLISM",     "KEGG_NITROGEN_METABOLISM",                        
              "KEGG_PROPANOATE_METABOLISM",     "KEGG_FATTY_ACID_METABOLISM",                      
              "KEGG_CITRATE_CYCLE_TCA_CYCLE",     "KEGG_GLUTATHIONE_METABOLISM",                     
              "KEGG_GALACTOSE_METABOLISM",     "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",           
              "KEGG_CYSTEINE_AND_METHIONINE_METABOLISM",    "KEGG_GLYCEROLIPID_METABOLISM",                  
              "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM",  "GOMF_FATTY_ACID_LIGASE_ACTIVITY",                
              "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",   "GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS")

data <- allDiff[gsva_use,]
data$X <- gsub("_"," ",data$X)

data <- data[order(data$t,decreasing = F),]
data$order <- 1:nrow(data)
data$color <- ifelse(data$t > 0, "up", "down")
data$X <- factor(data$X, levels = data$X)

colors <- c("up" = "#6E4B9E", "down" = "#D24B27")

###################Figure 1H-2#####################
ggplot(data,aes(x=X, y=t, fill = color)) + geom_bar(stat="identity",alpha=0.7) + 
  theme_classic() + coord_flip() + scale_fill_manual(values = colors) + 
  ylab("T-value")+xlab("Pathways")


#Figure 1H-1 could be drawn with matrix=coding_combat2_mat

matrix <- coding_combat2_mat
res_h_RNA <- gsva(as.matrix(matrix), geneset, method = "ssgsea", parallel.sz = 0, verbose = T)
res2_RNA <- gsva(as.matrix(matrix), geneset2, method = "ssgsea", parallel.sz = 0, verbose = T)
res5_RNA <- gsva(as.matrix(matrix), geneset5, method = "ssgsea", parallel.sz = 0, verbose = T)

NT <- sapply(strsplit(colnames(matrix),"_"),"[",2)

group<-model.matrix(~factor(NT))
colnames(group) <- c("normal","tumor")
fit<-lmFit(rbind(res_h,res2,res5),group)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
allDiff$X <- row.names(allDiff)
allDiff$X <- tolower(allDiff$X)

data <- allDiff[gsva_use,]
data$X <- gsub("_"," ",data$X)

data <- data[order(data$t,decreasing = F),]
data$order <- 1:nrow(data)
data$color <- ifelse(data$t > 0, "up", "down")
data$X <- factor(data$X, levels = data$X)

colors <- c("up" = "#6E4B9E", "down" = "#D24B27")

###################Figure 1H-1#####################
ggplot(data,aes(x=X, y=t, fill = color)) + geom_bar(stat="identity",alpha=0.7) + 
  theme_classic() + coord_flip() + scale_fill_manual(values = colors) + 
  ylab("T-value")+xlab("Pathways")



