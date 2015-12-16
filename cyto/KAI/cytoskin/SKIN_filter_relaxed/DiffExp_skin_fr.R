###Script finds differentially expressed genes in the filtered dataset

library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/SKIN_filter_relaxed/")
###########################Put Controlls, Pso, Vit in one matrix##########
load(file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.skin.Rdata")
data=data.skin

pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="All data after standardization")

#### PSick vs PHealthy
p.test=data[,24:86] #
p=p.test
cols=colnames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(p)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL")
cont.vals =c("PNL")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
ph_pt=tt
genes=rownames(ph_pt)
ph_pt_heatm=p[rownames(p)%in%genes,]
colnames(ph_pt_heatm)=cols
ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(ph_pt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in psoriasis lesional skin (PL) vs nonlesional skin (PNL)")
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PLvsPNLup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PLvsPNLdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PhvsPt=Gene
#[1] "IL17A"   "PI3"     "S100A8"  "S100A9"  "CCL27"   "FOXP3"   "IL1F6"   "IL1b"    "CXCL2"   "IL8b"    "CCL20"  
#[12] "AIM2"    "CCL2"    "CXCL10"  "IFNG"    "PYCARD"  "IL1RN"   "IL22RA1" "TNF"     "IFNGR"   "IL22RA2" "IFIH1"  
#[23] "NLRP1"


########### PL vs C
t=data[,1:55]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(t, mm)
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
print(tt) 
ph_c=tt
genes=rownames(ph_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in lesional skin in psoriatic patients (PS) vs controls (C)")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PLvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PlvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

PLvsC=Gene

#[1] "IL10"    "IL17F"   "IL17A"   "PI3"     "S100A9"  "S100A8"  "IL1F6"   "LCN2"    "CCL20"   "CXCL10"  "IFNG"   
#[12] "CXCL1"   "CCL27"   "CXCL8"   "CTLA4"   "IL1RN"   "CXCL2"   "CCL2"    "IL1b"    "FOXP3"   "AIM2"    "IL22RA1"
#[23] "TNF"     "DEFB1"   "CCL5"    "IFIH1"   "IL26"    "EOMES"  
############# PNL vs C

t=data[, c(1:23,56:86)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PNL")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
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
#data_heatm=as.data.frame(t(data_heatm))
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in nonlesional skin in psoriatic parients (PH) vs controls (C)")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PNLvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PNLvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PNLvsC=Gene
#[1] "IL1F6"   "CXCL10"  "IFNG"    "IL22RA1" "CCL20"   "DEFB1" 

PLvsPNL=PhvsPt
PNLvsC
PLvsC

PLvsPNL[PLvsPNL%in%PLvsC]
PLvsC[!PLvsC%in%PLvsPNL]
