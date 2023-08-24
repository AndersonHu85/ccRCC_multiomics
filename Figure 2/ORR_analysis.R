library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
library(tidyr)

out.prefix <- "./ORR"

do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}

################### Figure 2G ######################

Immune <- snRNA[,snRNA$celltype%in%c("Myeloid","T & NK")]

meta <- Immune@meta.data
meta$loc <- meta$class

meta$meta.cluster <- meta$celltype2

OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                out.prefix=sprintf("%s.Immune_cell",out.prefix),
                                pdf.width=4,pdf.height=8,verbose=1
)

a=OR.immune.list[["OR.dist.tb"]]
a <- as.data.frame(a)
rownames(a) <- a$rid
a <- a[,-1]
a <- na.omit(a)

col <- viridis(11,option = "D")

b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
b <- spread(b,key = "cid", value = "adj.p.value")
b <- data.frame(b[,-1],row.names = b$rid)
b <- b[rownames(a),]

b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
           ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                  ifelse(b<0.001&(a>1.5|a<0.5),"***",
                         ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))

bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))

pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
         display_numbers = b,number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))


################### Extended Data Figure 3G ######################

EC <- snRNA[,snRNA$celltype%in%c("EC")]

meta <- EC@meta.data
meta$loc <- meta$class

meta$meta.cluster <- meta$celltype2

OR.EC.list <- do.tissueDist(cellInfo.tb=meta,
                                out.prefix=sprintf("%s.EC",out.prefix),
                                pdf.width=4,pdf.height=8,verbose=1
)

a=OR.EC.list[["OR.dist.tb"]]
a <- as.data.frame(a)
rownames(a) <- a$rid
a <- a[,-1]
a <- na.omit(a)

col <- viridis(11,option = "D")

b <- OR.EC.list$count.dist.melt.ext.tb[,c(1,2,6)]
b <- spread(b,key = "cid", value = "adj.p.value")
b <- data.frame(b[,-1],row.names = b$rid)
b <- b[rownames(a),]

b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
           ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                  ifelse(b<0.001&(a>1.5|a<0.5),"***",
                         ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))

bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))

pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
         display_numbers = b[,],number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))

