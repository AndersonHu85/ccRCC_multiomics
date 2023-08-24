library(Cardinal)
library(ggpubr)
library(ggsci)
library(ggrepel)

############## Figure 3H ###############
#Spatial metabolomics visualization was conducted for each individual sample. We illustrate this process using the example of R29_T
##ImzML files are all accessible at zenodo and Mendeley

R29_T_neg <- readImzML(name="R29-T-neg", folder = "./")
image(R29_T_neg,mz=132.0302, normalize.image="none",colorscale=c(rep("black",2),col.map("jet",n=99)),
      contrast.enhance="histogram",smooth.image = "gaussian")
dev.off()

R29_T_pos <- readImzML(name="R29-T-pos", folder = "./")

image(R29_T_pos,mz=137.04581, normalize.image="none",colorscale=c(rep("black",2),col.map("jet",n=99)),
      contrast.enhance="histogram",smooth.image = "gaussian")
