###Script finds differentially expressed genes in the filtered dataset
#with removed cols and rows containing >50% NAs. The data is logtransformed

library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/diffexp")
###########################Put Controlls, Pso, Vit in one matrix##########
load(file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data_log2_filt.RData")
data=data_log2_filt
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
case.vals=c("PS")
cont.vals =c("PH")
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
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
ph_pt=tt
genes=rownames(ph_pt)
ph_pt_heatm=p[rownames(p)%in%genes,]
colnames(ph_pt_heatm)=cols
ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(ph_pt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in psoriasis lesional skin (PS) vs nonlesional skin (PH)")
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PSvsPHup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PSvsPHdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PhvsPt=Gene

##### Vh vs VT using all the data
v.test=data[,88:117]
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
case.vals=c("VS")
cont.vals =c("VH")
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
head(tt) 
vh_vt=tt
genes=rownames(vh_vt)
vh_vt_heatm=v[rownames(v)%in%genes,]
colnames(vh_vt_heatm)=cols
vh_vt_heatm=vh_vt_heatm[ , ! apply( vh_vt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(vh_vt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed genes in VHvsVT")
####no diff expressed genes

########### PS vs C
t=data[,1:55]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PS")
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
write.table(up,file="PSvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PSvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

PHvsC=Gene
############# PHealthy vs C

t=data[, c(1:23,56:86)]
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
write.table(up,file="PHvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PHvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PTvsC=Gene

############# VSick vs C
t=data[, c(1:23,87:102)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VS")
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
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in lesional skin in vitiligo patients (VS) vs controls (C)")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="VSvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="VSvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

VHvsC=Gene
##### VHealthy vs C

############ PH&VH vs C ######

t=data[, c(1:59,95:110)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VH","PH")
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
ph.vh_c=tt
genes=rownames(ph.vh_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Vitiligo.Differentially expressed in PH&VHvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PH_VH_vs_C_up_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PH_VH_vs_C_down_filt.txt",sep="\t",row.names=F,quote=FALSE)
PH_VH_vs_C=Gene

########### PT&VT vs C ######

t=data[, c(1:25,60:94,111:126)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VT","PT")
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
pt.vt_c=tt
genes=rownames(pt.vt_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Vitiligo. Differentially expressed in PT&VTvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PT_VT_vs_C_up_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PT_VT_vs_C_down_filt.txt",sep="\t",row.names=F,quote=FALSE)

PT_VT_vs_C=Gene
###### PT vs VT

t=data[, c(60:94,111:126)]

attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PT")
cont.vals =c("VT")
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

print(tt) 
pt_vt=tt
genes=rownames(pt_vt)
pt_vt_heatm=t[rownames(t)%in%genes,]
pheatmap(pt_vt_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Vitiligo.Differentially expressed in PTvsVT")
######no differentially expressed genes between PT and VT


###### PH vs VH

t=data[, c(26:59,95:110)]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PH")
cont.vals =c("VH")
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

print(tt) 
ph_vh=tt
genes=rownames(ph_vh)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Vitiligo. Differentially expressed in PHvsVH")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PHvsVHup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PHvsVTdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PHvsVH=Gene

######venn diagram

require(gplots)
venn(list(PHvsPT=PhvsPt,PHvsC=PHvsC,PTvsC=PTvsC,VHvsC=VHvsC,PHvsVH=PHvsVH ))

#v <- venneuler(c(PHvsPT=length(PhvsPt), PHvsC=length(PHvsC),VHvsC=length(VHvsC), "PHvsPT&PHvsC"=length(!PhvsPt%in%PHvsC), "PHvsPT&VHvsC"=length(!PhvsPt%in%VHvsC),"PHvsC&VHvsC"=length(!PHvsC%in%VHvsC) ))
#plot(v)


#See phenoLimma.R for the next analysis.


