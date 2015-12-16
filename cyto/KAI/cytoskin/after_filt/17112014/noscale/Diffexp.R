
#get the structures DataControllsSkin DataPsoriasisSkin DataVitiliigoSkin
#source("~/Documents/KAI/cytoskin/CT_files.R")
#source("~/Documents/KAI/cytoskin/CreateMetadata.R")
library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp")
load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")
data_filt=data_log2_filt
##Controlls
ctrls=data_filt[,1:24]

##PSO.h
ph=data_filt[,25:56]

##PSO.t
pt=data_filt[,57:87]

##VIT.h
vh=data_filt[,88:103]

##VIT.t
vt=data_filt[,104:118]

###########################Put Controlls, Pso, Vit in one matrix##########
data=cbind(ctrls,ph,pt,vh,vt)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "row", main="All data after standardization")
summary(data)

#### PH vs PT
p.test=data[,25:87]
p=p.test
cols=colnames(p)
colnames(p)=cols
#rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(p)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PH")
cont.vals =c("PT")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
print(thr.p.value)
tt 
ph_pt=tt
genes=rownames(ph_pt)
ph_pt_heatm=p[rownames(p)%in%genes,]
colnames(ph_pt_heatm)=cols
ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(ph_pt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriasis.Differentially expressed in PHvsPT",
filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/PHvsPT.png")
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PHvsPTup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PHvsPTdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PhvsPt=Gene
##### Vh vs VT using all the data
v.test=data[,88:118]
v=v.test
rows=rownames(v)
cols=colnames(v)
colnames(v)=cols
rownames(v)=rows
attribute=colnames(v)

attribute=str_replace(attribute,"0.*","")
coldata=colnames(v)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VH")
cont.vals =c("VT")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(v)=coldata
head(v)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(v, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)

if(nrow(tt) > 0){
  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
}

print(thr.p.value)
tt
vh_vt=tt
genes=rownames(vh_vt)
vh_vt_heatm=v[rownames(v)%in%genes,]
colnames(vh_vt_heatm)=cols
vh_vt_heatm=vh_vt_heatm[ , ! apply( vh_vt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(vh_vt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed genes in VHvsVT",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/VHvsVT.png")
####no diff expressed genes

########### PH vs C
t=data[,1:56]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PH")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(t, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1,, lfc = thr.lfc, adjust.method = adj.method, number =10000000)
if(nrow(tt) > 0){
  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
}

print(thr.p.value)
print(tt) 
ph_c=tt
genes=rownames(ph_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriasis. Differentially expressed in PHvsC",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/PHvsC.png")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PHvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PHvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

PHvsC=Gene
############# PT vs C

t=data[, c(1:24,57:87)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PT")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(t, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
#if(nrow(tt) > 0){
#  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
#}

print(thr.p.value)
print(tt) 
pt_c=tt
genes=rownames(pt_c)
data_heatm=t[rownames(t)%in%genes,]
data_heatm=as.data.frame(t(data_heatm)
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriasis.Differentially expressed in PTvsC",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/PTvsC.png")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PTvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PTvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PTvsC=Gene

############# VH vs C
t=data[, c(1:24,88:103)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VH")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(t, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
#if(nrow(tt) > 0){
#  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
#}

print(thr.p.value)
print(tt) 
vh_c=tt
genes=rownames(vh_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed in VHvsC",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/VHvsC.png")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="VHvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="VHvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

VHvsC=Gene
#no diff expressed genes
##### VT vs C
t=data[, c(1:24,104:118)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VT")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(t, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
#if(nrow(tt) > 0){
#  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
#}

print(thr.p.value)
print(tt) 
vt_c=tt
genes=rownames(vt_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed in VTvsC",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/VTvsC.png")

tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="VTvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="VTvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
VTvsC=Gene



#See phenoLimma.R for the next analysis.


