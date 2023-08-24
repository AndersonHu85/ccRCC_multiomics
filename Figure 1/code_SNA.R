library(maftools)
library(sigminer)
library(NMF)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(ggplot2)
library(agricolae)
library(ggsignif)
library(dplyr)
options(stringsAsFactors = F)

annovar.laml_demo <- annovarToMaf(annovar = "./mutation.maf", #####MAF file of TJ-RCC
                                               refBuild = 'hg38',
                                               tsbCol = 'Tumor_Sample_Barcode', 
                                               MAFobj = T)



#######################Figure 1B#################################
###########Global Mutation landscape of TJ-RCC###################

color <- ArchR::ArchRPalettes[5][[1]][6:1]#[-3]
names(color) <- c("C>T","C>G", "C>A", "T>A", "T>C", "T>G")

vc_cols = ArchR::ArchRPalettes[1][[1]][8:1]#[-4]
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

oncoplot(maf = annovar.laml_demo,draw_titv=T,
         genes = c("VHL","PBRM1","SETD2","BAP1","KDM5C","PCLO","CSMD3","SYNE1","XIRP2","MTOR","TP53",
                   "PTEN","NFE2L2","KEAP1","PKHD1L1","RB1"),
         sepwd_samples = 0,
         fontSize = 0.8,
         keepGeneOrder = F,
         colors = vc_cols,
         bgCol = "grey90",
         titv_col=color,
         altered = T,
         drawBox = T,
         removeNonMutated=F,
         borderCol = NA,
         showTumorSampleBarcodes = F)


##############Extended data Figure 1B#################
annovar.laml <- annovarToMaf(annovar = "./merge_490.txt", ###merge_490.txt is the file containing SNV and INDEL information of 3 Chinese cohort 
                             refBuild = 'hg38',
                             tsbCol = 'Tumor_Sample_Barcode', 
                             MAFobj = T)

somaticInteractions(maf = annovar.laml, 
                    genes = c("VHL","PBRM1","TTN","BAP1","SETD2","KDM5C","MUC16","PCLO","MTOR","XIRP2",
                              "PKHD1L1","CSMD3","ABCA13","MALRD1","SYNE1","TP53",
                              "PTEN")
)


#######################Figure 1B#################################
###########Global Mutation landscape of TJ-RCC###################

dir="../../SNV/merge/"   ####dir is the direction deposited vcf files
list <- list.files(dir,pattern = "hg38_multianno.vcf", 
                   full.names = TRUE) %>% grep("indel",.,invert = T,value = T) 

maf <- read_vcf(list)
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE
)
mt_tally$nmf_matrix[1:5, 1:5]

mt_sig <- sig_auto_extract(mt_tally$nmf_matrix,
                           K0 = 10, nrun = 10,
                           strategy = "stable"
)

##############Extended data Figure 1C#################
###############six mutation signature#################

sim_v3 <- get_sig_similarity(mt_sig, sig_db = "SBS")
show_sig_profile(mt_sig, mode = "SBS", style = "cosmic",x_label_angle = 90,palette = color)+scale_color_npg()


##############Extended data Figure 1D#################
###########NMF cluster of mutation signature#################

grp <- get_groups(mt_sig,method = "exposure")
grp_label <- grp$group  ###extract NMF cluster
names(grp_label) <- grp$sample


show_sig_exposure(mt_sig, style = "cosmic", groups = grp_label,
                  palette = pal_npg()(7), rm_grid_line = T,rm_panel_border=F,rm_space=T)



##############Extended data Figure 1E#################
###########annotion of mutation signature#################
mark_sig <- colnames(sim_v3$similarity)[rowMax(t(sim_v3$similarity))>0.7]
mark_sig <- colnames(sim_v3$similarity)[colnames(sim_v3$similarity)%in%mark_sig]
gene_pos <- which(colnames(sim_v3$similarity) %in% mark_sig)
row_anno <-  rowAnnotation(mark_sig = anno_mark(at = gene_pos, 
                                                 labels = mark_sig))

col <- ArchR::ArchRPalettes[19][[1]]
Heatmap(t(sim_v3$similarity), 
        col = colorRamp2(seq(0,1,length.out = 9),col),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = F,
        use_raster=T,
        border = "BLACK",
        show_row_dend = T,
        na_col = "grey85",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D",
        column_gap = unit(1, "mm"),
        right_annotation = row_anno,
        row_gap = unit(0, "mm"),
        column_title = NULL) # 不需要列标题)


##################Extended data Figure 1F######################

library(BSgenome.Hsapiens.UCSC.hg38) 

df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # chromosome name
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# chromosome length
)
df$chromNum <- 1:length(df$chromName) 

df <- df[1:22,]###chr1-22

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength))
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

scores <- read.table("./gistic/scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)

chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]

range(scores$G.score)

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1

range(scores$G.score)


df$ypos <- rep(c(0.7,0.8),11)

ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-1,1)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )



###############Figure 1C###################

SNA <- annovar.laml@variants.per.sample
SNA <- as.data.frame(SNA)
row.names(SNA) <- SNA$Tumor_Sample_Barcode

group <- as.data.frame(grp_label)
group$sample <- rownames(group)
group$sample <- gsub("_annovar.hg38_multianno","",group$sample)
group$sample <- gsub("Tp","T",group$sample)
group$sample <- gsub("_snp","",group$sample)

rownames(group) <- group$sample
SNA$Tumor_Sample_Barcode <- as.character(SNA$Tumor_Sample_Barcode)
SNA <- SNA[row.names(group),]
rownames(SNA) <- rownames(group)
SNA$Tumor_Sample_Barcode <- rownames(SNA)
SNA[is.na(SNA)=="TRUE"]=0

group <- group[(SNA$Tumor_Sample_Barcode),]
all(group$sample==SNA$Tumor_Sample_Barcode)
SNA$SBS="non-AA"
SNA$SBS[which(group$grp_label=="2")]="AA"

SNA$group <- "TJ"
SNA$group[grep("_T",SNA$Tumor_Sample_Barcode)]="Pek"
SNA$group[grep("ccRCC",SNA$Tumor_Sample_Barcode)]="FUD"

SBS <- ArchR::ArchRPalettes[2][[1]][c(14,16)]
names(SBS) <- c("AA","non-AA")

ggplot(SNA,aes(SBS,Variants,color=SBS)) + 
  geom_boxplot(outlier.colour = NA) +
  theme_bw() +
  geom_signif(comparisons = list(c("AA","non-AA"))) + 
  geom_jitter(data = SNA,aes(shape=group,alpha=0.7)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  )  + scale_color_manual(values = SBS)



#########################Figure 1E#########################
AA_sample <- SNA$Tumor_Sample_Barcode[SNA$SBS=="AA"]
vcf_AA <- annovar.laml@data[annovar.laml@data$Tumor_Sample_Barcode%in%AA_sample,]
vcf_non_AA <- annovar.laml@data[!annovar.laml@data$Tumor_Sample_Barcode%in%AA_sample,]

mutation_map <- c("T>A" = "T>A", "A>T" = "T>A",
                  "C>T" = "C>T", "G>A" = "C>T",
                  "C>G" = "C>G", "G>C" = "C>G",
                  "C>A" = "C>A", "G>T" = "C>A",
                  "T>C" = "T>C", "A>G" = "T>C",
                  "T>G" = "T>G", "A>C" = "T>G")

vcf_AA$AA <- mutation_map[paste(vcf_AA$Reference_Allele, vcf_AA$Tumor_Seq_Allele2, sep = ">")]
vcf_non_AA$AA <- mutation_map[paste(vcf_non_AA$Reference_Allele, vcf_non_AA$Tumor_Seq_Allele2, sep = ">")]


gene_hub <- c("VHL","PBRM1","SETD2","BAP1","KDM5C")
genes = c("TTN","MUC16","PCLO","MTOR","XIRP2",
          "PKHD1L1","CSMD3","ABCA13","MALRD1","SYNE1")
vcf_AA_hub <- vcf_AA[vcf_AA$Hugo_Symbol%in%gene_hub&is.na(vcf_AA$AA)=="FALSE",]
vcf_AA_genes <- vcf_AA[vcf_AA$Hugo_Symbol%in%genes&is.na(vcf_AA$AA)=="FALSE",]
pie_data <- table(vcf_AA_hub$AA)%>%data.frame()
pie_data[4,2]/sum(pie_data$Freq)

vcf_non_AA_hub <- vcf_non_AA[vcf_non_AA$Hugo_Symbol%in%gene_hub&is.na(vcf_non_AA$AA)=="FALSE",]
vcf_non_AA_genes <- vcf_non_AA[vcf_non_AA$Hugo_Symbol%in%genes&is.na(vcf_non_AA$AA)=="FALSE",]

pie_data <- table(vcf_non_AA_genes$AA)%>%data.frame()
ggplot(data=pie_data, mapping=aes(x="Freq",y=Freq,fill=Var1))+
  geom_bar(stat="identity",width=0.5,position='stack')+
  coord_polar("y", start=0)+
  scale_fill_manual(values=color)+theme_bw()

pie_data <- table(vcf_AA_genes$AA)%>%data.frame()
ggplot(data=pie_data, mapping=aes(x="Freq",y=Freq,fill=Var1))+
  geom_bar(stat="identity",width=0.5,position='stack')+
  coord_polar("y", start=0)+
  scale_fill_manual(values=color)+theme_bw()

pie_data <- table(vcf_non_AA_hub$AA)%>%data.frame()
ggplot(data=pie_data, mapping=aes(x="Freq",y=Freq,fill=Var1))+
  geom_bar(stat="identity",width=0.5,position='stack')+
  coord_polar("y", start=0)+
  scale_fill_manual(values=color)+theme_bw()

pie_data <- table(vcf_AA_hub$AA)%>%data.frame()
ggplot(data=pie_data, mapping=aes(x="Freq",y=Freq,fill=Var1))+
  geom_bar(stat="identity",width=0.5,position='stack')+
  coord_polar("y", start=0)+
  scale_fill_manual(values=color)+theme_bw()


#######################Figure 1F#########################
vcf_AA_hub$Tumor_Sample_Barcode <- as.character(vcf_AA_hub$Tumor_Sample_Barcode)
vcf_non_AA_hub$Tumor_Sample_Barcode <- as.character(vcf_non_AA_hub$Tumor_Sample_Barcode)
vcf_AA_genes$Tumor_Sample_Barcode <- as.character(vcf_AA_genes$Tumor_Sample_Barcode)
vcf_non_AA_genes$Tumor_Sample_Barcode <- as.character(vcf_non_AA_genes$Tumor_Sample_Barcode)

pie_data <- table(vcf_non_AA_hub$AA,vcf_non_AA_hub$Tumor_Sample_Barcode)%>%data.frame()
pie_data2 <- table(vcf_AA_hub$AA,vcf_AA_hub$Tumor_Sample_Barcode)%>%data.frame()
number <- table(as.character(vcf_non_AA_hub$Tumor_Sample_Barcode))%>%data.frame()
number2 <- table(as.character(vcf_AA_hub$Tumor_Sample_Barcode))%>%data.frame()

pie_data <- pie_data[pie_data$Var1=="T>A",]
pie_data2 <- pie_data2[pie_data2$Var1=="T>A",]
pie_data$Freq <- pie_data$Freq/number$Freq
pie_data2$Freq <- pie_data2$Freq/number2$Freq
pie_data$group <- "non-AA"
pie_data2$group <- "AA"
pie_merge <- rbind(pie_data,pie_data2)

ggplot(pie_merge,aes(group,Freq,color=group)) + geom_boxplot(outlier.colour = NA) +theme_bw()+
  geom_signif(comparisons = list(c("AA","non-AA")),test = "t.test",
  ) + geom_jitter(data = pie_merge,aes(alpha=0.7),width = 0.2,height = 0.05) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  )  + scale_color_manual(values = SBS)


pie_data <- table(vcf_non_AA_genes$AA,vcf_non_AA_genes$Tumor_Sample_Barcode)%>%data.frame()
pie_data2 <- table(vcf_AA_genes$AA,vcf_AA_genes$Tumor_Sample_Barcode)%>%data.frame()
number <- table(as.character(vcf_non_AA_genes$Tumor_Sample_Barcode))%>%data.frame()
number2 <- table(as.character(vcf_AA_genes$Tumor_Sample_Barcode))%>%data.frame()

pie_data <- pie_data[pie_data$Var1=="T>A",]
pie_data2 <- pie_data2[pie_data2$Var1=="T>A",]
pie_data$Freq <- pie_data$Freq/number$Freq
pie_data2$Freq <- pie_data2$Freq/number2$Freq
pie_data$group <- "non-AA"
pie_data2$group <- "AA"
pie_merge <- rbind(pie_data,pie_data2)

ggplot(pie_merge,aes(group,Freq,color=group)) + geom_boxplot(outlier.colour = NA) +theme_bw()+
  geom_signif(comparisons = list(c("AA","non-AA")),test = "t.test",
  ) + geom_jitter(data = pie_merge,aes(alpha=0.7),width = 0.2,height = 0.05) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black",size=0.5), 
    axis.title.y=element_text(size=14),axis.text.y=element_text(size=14),
    axis.title.x=element_blank(),axis.text.x=element_text(size=12)
  )  + scale_color_manual(values = SBS)

