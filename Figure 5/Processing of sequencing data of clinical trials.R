############ datasets of clinical trial are not provided in source data due to policy limitation. Data deposited in EGA needs DAC review
############relating to figure 5#####################3

library(survminer)
library(survival)
library(data.table)
library(CMScaller)
library(tidyverse)
library(forestplot)
library(meta)


matrix_IMmotion <- fread("./EGAF00004859518/IMmotion151.expression.data.TPM.anon.20201106.csv")
matrix_IMmotion <- matrix_IMmotion[is.na(matrix_IMmotion$symbol)=="FALSE",]
matrix_IMmotion <- matrix_IMmotion[!duplicated((matrix_IMmotion$symbol)),]
matrix_IMmotion <- data.frame(matrix_IMmotion[,4:ncol(matrix_IMmotion)],row.names = matrix_IMmotion$symbol)

clinical_IMmotion <- fread("./EGAF00004859519/IMmotion151_clinical_anon_20201106.csv")
clinical_IMmotion <- data.frame(clinical_IMmotion)
colnames(clinical_IMmotion) <- clinical_IMmotion[4,]
clinical_IMmotion <- data.frame(clinical_IMmotion[5:nrow(clinical_IMmotion),])
clinical_IMmotion$RNASEQ_SAMPLE_ID <- gsub("-",".",clinical_IMmotion$RNASEQ_SAMPLE_ID)

clinical_IMmotion <- clinical_IMmotion[order(clinical_IMmotion$ARM),]
sample <- intersect(clinical_IMmotion$RNASEQ_SAMPLE_ID,colnames(matrix_IMmotion))
clinical_IMmotion <- clinical_IMmotion[clinical_IMmotion$RNASEQ_SAMPLE_ID%in%sample,]
clinical_IMmotion$PFS_CENSOR[clinical_IMmotion$PFS_CENSOR=="TRUE"]=1
clinical_IMmotion$PFS_CENSOR[clinical_IMmotion$PFS_CENSOR=="FALSE"]=0
clinical_IMmotion$PFS_CENSOR <- as.numeric(clinical_IMmotion$PFS_CENSOR)
matrix_IMmotion <- matrix_IMmotion[,sample]

clinical_IMmotion$PFS_MONTHS <- as.numeric(clinical_IMmotion$PFS_MONTHS)

#####annotation using ntp algorithm######################
tmp_org <- read.csv("tmp_org.csv",row.names = 1)
colnames(tmp_org) <- c("probe","class")
tmp_org_use <- tmp_org[tmp_org$probe %in%rownames(matrix_IMmotion),]

emat <- ematAdjust(matrix_IMmotion, normMethod = 'RLE')
res <- ntp(emat, tmp_org_use, doPlot=TRUE, nPerm=1000, seed = 42,nCores =30)
clinical_IMmotion$prediction <- res[sample,]$prediction

clinical_IMmotion_sub <- clinical_IMmotion[clinical_IMmotion$ARM!="atezo_bev",]

fit <- survfit(Surv(clinical_IMmotion_sub$PFS_MONTHS, 
                    clinical_IMmotion_sub$PFS_CENSOR) ~ prediction, data = clinical_IMmotion_sub)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                size=1,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(), # Change ggplot2 theme
                palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF" )
)
p


clinical_IMmotion_sub <- clinical_IMmotion[clinical_IMmotion$ARM=="atezo_bev",]

fit <- survfit(Surv(clinical_IMmotion_sub$PFS_MONTHS, 
                    clinical_IMmotion_sub$PFS_CENSOR) ~ prediction, data = clinical_IMmotion_sub)
p <- ggsurvplot(fit,
                pval = T ,conf.int = F,censor = T, censor.size = 3,
                axes.offset = T,
                size=1,
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_classic(), # Change ggplot2 theme
                palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF" )
)
p

top_anno <- HeatmapAnnotation(df = clinical_IMmotion[,c(9,92)], col = 
                                list(ARM=c("atezo_bev"="#FB8072", "sunitinib"="#80B1D3"),
                                     prediction=c("IM1"="#E64B35FF", "IM2"="#4DBBD5FF", "IM3"="#00A087FF", "IM4"="#3C5488FF" )),
                              show_legend = T)


Heatmap(t(scale(t(matrix_IMmotion[tmp_org_use$probe,sample]))), 
        col = colorRamp2(seq(-2,2,length.out = 11),col),
        cluster_rows = F,
        cluster_columns = F,   
        row_names_gp = gpar(fontsize = 10),
        clustering_distance_rows = "pearson",
        clustering_method_rows = "ward.D",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "ward.D",
        show_column_names = F,
        show_row_names = F,
        border = "black",
        row_split = tmp_org2$class,
        column_split = clinical_IMmotion$prediction,
        column_gap = unit(0, "mm"),
        row_gap = unit(0, "mm"),
        top_annotation = top_anno, # 在热图上边增加注释
        #column_title = NULL
) # 不需要列标题)


clinical_IMmotion_sub <- lapply(c("IM1", "IM2", "IM3", "IM4"), function(i) {
  subset <- clinical_IMmotion[clinical_IMmotion$prediction == i, ]
  subset$ARM[subset$ARM == "atezo_bev"] <- "2_Combo"
  subset$ARM[subset$ARM == "sunitinib"] <- "1_Sunitinib"
  fit <- coxph(Surv(subset$PFS_MONTHS, subset$PFS_CENSOR) ~ ARM, data = subset)
  fit <- summary(fit)
  data.frame(fit$coefficients, fit$conf.int[], row.names = i)
})

HR_List <- do.call(rbind, clinical_IMmotion_sub)
HR_List$group <- rownames(HR_List)

tabletext <- data.frame(
  predicted_cluster =  c("IM1", "IM2", "IM3", "IM4"),
  atezo_bev = c("66", "134", "68", "139"),
  sunitinib =  c("68", "139", "78", "131"),
  mean  = HR_List$exp.coef., 
  lower = HR_List$lower..95,
  upper = HR_List$upper..95
)

tabletext$`HR (95% CI)` <- ifelse(is.na(tabletext$mean), "",
                                  sprintf("%.2f (%.2f to %.2f)",
                                          tabletext$mean, tabletext$lower, tabletext$upper))

tabletext$` ` <- paste(rep(" ", 8), collapse = " ")

library(forestplot)
library(forestploter)

p <- forest(data = tabletext[,c(1:3,8,7)], 
            lower = tabletext$lower, 
            upper = tabletext$upper, 
            est = tabletext$mean, 
            ci_column = 4 ,
            sizes = tabletext$mean, 
            ref_line = 1, 
            xlim = c(0,2), 
            ticks_at = c(0,1,2,3),
            
)

print(p)

##################### three cohorts were processed following the same workflow ################# 
