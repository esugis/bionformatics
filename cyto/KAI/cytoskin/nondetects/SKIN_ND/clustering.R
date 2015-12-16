#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
library(plyr)
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

library(amap)
genes=data.skin
hc.data=hcluster(genes, method = "correlation", diag = FALSE, upper = FALSE,
                 link = "complete", members = NULL, nbproc = 2,
                 doubleprecision = TRUE)

plot(hc.data)
ct=cutree(hc.data, k=6)

save(hc.data,ct,file="hclust.RData")
#mydata= na.omit(data.skin)
#dim(mydata)
#wss <- nrow(mydata-1)*sum(apply(mydata,2,var))
#for (i in 1:) wss[i] <- sum(kmeans(mydata, 
#                                     centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#     ylab="Within groups sum of squares")
