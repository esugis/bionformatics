###Script finds differentially expressed genes in the filtered dataset
#with removed cols and rows containing >50% NAs. The data is logtransformed

library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/updated_names_PL_PLN/")
###########################Put Controlls, Pso, Vit in one matrix##########
load(file="data.pca.Rdata")
data=data.pca
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="All data log transformed")

#### PL vs PNL
up=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSLvsPSNLup_filt.txt",sep="\t",header=T)
down=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSLvsPSNLdown_filt.txt",sep="\t",header=T)
genes=c(as.vector(up[,1]), as.vector(down[,1]))
ph_pt_heatm=data[rownames(data)%in%genes,25:93]
pheatmap(ph_pt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in psoriasis lesional skin (PL) vs nonlesional skin (PNL)")

##### Vh vs VT using all the data
####no diff expressed genes

########### PL vs C

up=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSLvsCup_filt.txt",sep="\t",header=T)
down=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSLvsCdown_filt.txt",sep="\t",header=T)
genes=c(as.vector(up[,1]), as.vector(down[,1]))
t=data[,1:58]
data_heatm=t[rownames(t)%in%genes,]
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in lesional skin in psoriatic patients (PL) vs controls (C)")

############# PNL vs C
up=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSNLvsCup_filt.txt",sep="\t",header=T)
down=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/PSNLvsCdown_filt.txt",sep="\t",header=T)
genes=c(as.vector(up[,1]), as.vector(down[,1]))
t=data[, c(1:24,59:93)]
data_heatm=t[rownames(t)%in%genes,]
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in nonlesional skin in psoriatic parients (PNL) vs controls (C)")


############# VL vs C
up=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/VLvsCup_filt.txt",sep="\t",header=T)
down=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp/VLvsCdown_filt.txt",sep="\t",header=T)
genes=c(as.vector(up[,1]), as.vector(down[,1]))
t=data[, c(1:24,94:109)]
data_heatm=t[rownames(t)%in%genes,]
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in lesional skin in vitiligo parients (VL) vs controls (C)")


