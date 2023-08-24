DCCD_score <- list(IM2_sig=tmp_org$probe[tmp_org$class=="IM2"],
                   IM4_sig=tmp_org$probe[tmp_org$class=="IM4"])


TCGA_DCCD_score <- gsva(as.matrix(TCGA), DCCD_score, method = "ssgsea", parallel.sz = 0, verbose = T)

IM1_3_sample <- rownames(clinical_TCGA)[clinical_TCGA$prediction%in%c("IM1","IM3")]

clinical_TCGA_IM1_3 <- clinical_TCGA[IM1_3_sample,]
clinical_TCGA_IM1_3$IM2_score <- TCGA_DCCD_score["IM2_sig",IM1_3_sample]
clinical_TCGA_IM1_3$IM4_score <- TCGA_DCCD_score["IM4_sig",IM1_3_sample]

clinical_TCGA_IM1_3$identity_score <- scale(clinical_TCGA_IM1_3$IM4_score)-scale(clinical_TCGA_IM1_3$IM2_score)
clinical_TCGA_IM1_3$group_IM <- "IM4-like"
clinical_TCGA_IM1_3$group_IM[clinical_TCGA_IM1_3$identity_score<0] <- "IM2-like"

table(clinical_TCGA_IM1_3$group_IM)

clinical_TCGA_IM1_3 <- clinical_TCGA_IM1_3[order(clinical_TCGA_IM1_3$identity_score,decreasing = T),]

mat_use <- TCGA[,rownames(clinical_TCGA_IM1_3)]%>%t()%>%scale()%>%t()

sig <- intersect(c(DCCD_score$IM4_sig,DCCD_score$IM2_sig),rownames(TCGA))

Heatmap(mat_use[sig,rownames(clinical_TCGA_IM1_3)], 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        row_split = c(rep("1_IM4_sig",285),rep("2_IM2_sig",177)),
        column_split = c(rep("1_IM4-like",77),rep("2_IM2-like",90)),
        column_gap = unit(0, "mm"), 
        row_gap = unit(0, "mm"), 
        border = "black",
        column_title = NULL)

clinical_TCGA$group_IM2 <- clinical_TCGA$prediction
clinical_TCGA[rownames(clinical_TCGA_IM1_3),]$group_IM2 = clinical_TCGA_IM1_3$group_IM


################## Figure 6C #################
fit <- survfit(Surv(OS.time/30, OS) ~ group_IM, data = clinical_TCGA)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(base_line_size=1, base_rect_size=1), # Change ggplot2 theme
                palette = pal_npg()(4)[c(1,2,4,3)]
)
p

########################### Extended Data Figure 7C ###############################

fit <- survfit(Surv(PFI.time/30, PFI) ~ group_IM, data = clinical_TCGA)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(base_line_size=1, base_rect_size=1), # Change ggplot2 theme
                palette = pal_npg()(4)[c(1,2,4,3)]
)
p

##### adjusting IM subgroup of all the cohorts was performed following the same pipeline


#################### Figure 6G ##################

column <- table(clinical_TCGA$group_IM,clinical_TCGA$X9p)%>%data.frame()
column$group <- "09p"
column2 <- table(clinical_TCGA$group_IM,clinical_TCGA$X12p)%>%data.frame()
column2$group <- "12p"
column3 <- table(clinical_TCGA$group_IM,clinical_TCGA$X20p)%>%data.frame()
column3$group <- "20p"
column4 <- table(clinical_TCGA$group_IM,clinical_TCGA$X14q)%>%data.frame()
column4$group <- "14q"
column <- rbind(column,column2,column3,column4)


p1 <- ggplot(column,aes(factor(Var1,levels = c("IM2","IM2-like","IM4-like","IM4")
),
weight=Freq,fill=Var2))+
  geom_bar(position="fill")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = 
                      c("grey85","#de2d26","#2171B5")
  )
p1 +  facet_grid(~group,scales = "free", space = "free") + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       axis.line = element_line(colour = "black",size=1), 
       axis.ticks = element_line(size = 1, color="black"),
       axis.ticks.length = unit(.2, "cm")) + xlab("Adjusted IM group") + ylab("Proportion")






