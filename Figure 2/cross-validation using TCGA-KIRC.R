#####cross-validation using TCGA-KIRC, related to Figure 2#############

TCGA <- fread("./TCGA_KIRC_pre_processed.txt")
TCGA <- data.frame(TCGA[,-1],row.names = TCGA$V1)

tmp_org2 <- tmp_org[tmp_org$Symbol%in%rownames(TCGA),]

library(Biobase)
library(CMScaller)

colnames(tmp_org2) <- c("probe","class")

emat <- ematAdjust(TCGA[tmp_org2$probe,], normMethod = 'RLE')
res <- ntp(emat, tmp_org2, doPlot=TRUE, nPerm=1000, seed = 42,nCores =30) #######predict subtype using ntp algorithm

IM_group_TCGA <- res$prediction %>% data.frame()
rownames(IM_group_TCGA) <- colnames(TCGA)
colnames(IM_group_TCGA) <- "predicted_subtype"

top_anno <- HeatmapAnnotation(df = IM_group_TCGA, col = 
                                list(predicted_subtype=c("IM1"="#E64B35FF", "IM2"="#4DBBD5FF", "IM3"="#00A087FF", "IM4"="#3C5488FF" )),
                              show_legend = T)

col <- brewer.pal(11,"RdBu")[11:1]
heat2 <- Heatmap(t(scale(t(TCGA)))[tmp_org2$probe, ], 
                 col = colorRamp2(seq(-2,2,length.out = 11),col),
                 cluster_rows = F,
                 cluster_columns = F,   
                 row_names_gp = gpar(fontsize = 10),
                 show_column_names = F,
                 show_row_names = F,
                 border = "black",
                 row_split = tmp_org2$class,
                 column_split = res$prediction,
                 column_gap = unit(0, "mm"),
                 row_gap = unit(0, "mm"),
                 top_annotation = top_anno
)
heat2


clinical_TCGA <- read.table("TCGA_clinical.txt") ########load organized clinical data
clinical_TCGA <- clinical_TCGA[colnames(TCGA),] #### keep in the same order
all(rownames(clinical_TCGA)==colnames(TCGA))
#[1] TRUE

table(res$prediction,clinical_TCGA$prediction) 

#     IM1 IM2 IM3 IM4
# IM1  93   0   0   0
# IM2   0 159   0   0
# IM3   0   0  74   0
# IM4   0   0   0 116

########################### Extended Data Figure 3B ###############################
column <- as.data.frame(table(clinical_TCGA$prediction, clinical_TCGA$grade))

colnames(column) <- c("IM_group","grade","Freq")

grade=c("G1"="#fee5d9","G2"="#fcae91","G3"="#de2d26","G4"="#a50f15","GX"="grey85")
p1 <- ggplot(column,aes(IM_group,weight=Freq,fill=grade))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = grade) #[c(1,5,17,14,21)]
p1 + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################### Extended Data Figure 3C ###############################
column <- as.data.frame(table(clinical_TCGA$prediction, clinical_TCGA$stage))

colnames(column) <- c("IM_group","stage","Freq")

stage = c("stage i"="#C6DBEF", "stage ii"="#6BAED6", "stage iii"="#2171B5", "stage iv"="#08306B","not reported"="grey85")
p1 <- ggplot(column,aes(IM_group,weight=Freq,fill=stage))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = stage) #[c(1,5,17,14,21)]
p1 + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################### Extended Data Figure 3D ###############################

library(survminer)
library(survival)
library(cowplot)

os_sub <- clinical_TCGA[!clinical_TCGA$stage%in%c("stage i"),]

fit <- survfit(Surv(OS.time/30, OS) ~ prediction, data = clinical_TCGA)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(base_line_size=1, base_rect_size=1), # Change ggplot2 theme
                palette = pal_npg()(4)
)
p


fit <- survfit(Surv(PFI.time/30, PFI) ~ prediction, data = clinical_TCGA)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(base_line_size=1, base_rect_size=1), # Change ggplot2 theme
                palette = pal_npg()(4)
)
p


