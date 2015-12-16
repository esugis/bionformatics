#This script reads in the files computed in DiffExpPVC_removedC024.R 
#Script

#setwd("~/Documents/KAI/cytoskin/SKIN_filter_relaxed/")
library(stringr)
#controlls
load(file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.control_FR.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.psoriasis_FR.RData")     
###Vitiliigo#
load(file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.vitiliigo_FR.RData")
###PCA plot of the data



ctrls=DataControllsSkin
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:43]
ctrls <- apply(ctrls[,2:43],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisSkin)# psoriasis
rows=rownames(pso[2:43,])
cols=as.character(pso[1,])
ps <- apply(pso[2:43,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiliigoSkin)# vitiliigo
rows=rownames(vit[2:43,])
cols=as.character(vit[1,])
vi <- apply(vit[2:43,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

data.all=cbind(ctrls, ps,vi)
summary(data.all)
data.all[rownames(data.all)%in%c("IL17F","IL17A","IL10"),c(1:23,56:86)]=0.00000000000000001
data.log=log2(data.all)
data.skin=data.log

colnames(data.skin)=str_replace(colnames(data.skin),"PH","PL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"PT","PNL")#healthy
colnames(data.skin)=str_replace(colnames(data.skin),"VH","VL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"VT","VNL")#healthy
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))
save(data.skin, file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.skin.Rdata")

pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL, PNL, VL and VNL groups",
         filename="data_log_NAs.png")

#diff expression skin pso


