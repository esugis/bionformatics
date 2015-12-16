
#get the structures DataControllsSkin DataPsoriasisSkin DataVitiliigoSkin
source("~/Documents/KAI/cytoskin/CT_files.R")
source("~/Documents/KAI/cytoskin/CreateMetadata.R")
library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/after_filt/11112014")
#load("~/Documents/KAI/cytoskin/after_filt/31102014/data.RData")
#Controlls
ctrls=DataControllsSkin
rows=as.character(ctrls[,1])
cols=colnames(ctrls[,2:43])
ctrls <- apply(ctrls,2,as.numeric)
ctrls=ctrls[,2:43]
colnames(ctrls)=cols
rownames(ctrls)=rows
ctrls=t(ctrls)


##PSO.h
pso.h=t(DataPsoriasisSkin[1:34,])#haige psoriasis
rows=rownames(pso.h[2:43,])
cols=as.character(pso.h[1,])
ph <- apply(pso.h[2:43,],2,as.numeric)
colnames(ph)=cols
rownames(ph)=rows


#PSO.t
pso.t=t(DataPsoriasisSkin[35:69,])#terve psoriasis  
rows=rownames(pso.t[2:43,])
cols=as.character(pso.t[1,])
pt <- apply(pso.t[2:43,],2,as.numeric)
colnames(pt)=cols
rownames(pt)=rows


#VIT.h
vit.h=t(DataVitiliigoSkin[1:16,])#haige vitiliigo
rows=rownames(vit.h[2:43,])
cols=as.character(vit.h[1,])
vh <- apply(vit.h[2:43,],2,as.numeric)
colnames(vh)=cols
rownames(vh)=rows


#VIT.t
vit.t=t(DataVitiliigoSkin[17:32,])#terve vitiliigo
rows=rownames(vit.t[2:43,])
cols=as.character(vit.t[1,])
vt <- apply(vit.t[2:43,],2,as.numeric)
colnames(vt)=cols
rownames(vt)=rows   


############PSO H vs T

#p=t(DataPsoriasisSkin)#haige psoriasis
#rows=rownames(p[2:43,])
#cols=as.character(p[1,])
#p <- apply(p[2:43,],2,as.numeric)
##standardisation
#scale_rows = function(x){
#  m = apply(x, 1, mean, na.rm = T)
#  s = apply(x, 1, sd, na.rm = T)
#  return((x - m) / s)
#}
#p.st=scale_rows(p)
#p=p.st#

#colnames(p)=cols
#rownames(p)=rows
#attribute=colnames(p)

#attribute=str_replace(attribute,"0.*","")
#coldata=colnames(p)
#coldata=str_replace(coldata,"0.*","")
#attribute = factor(attribute)
#attribute2 = rep("Muu", length(attribute))
## Case-control variables and corresponding group values 
#case.vals=c("PH")
#cont.vals =c("PT")
#attribute2[attribute %in% c(case.vals)] = "Group1"
#attribute2[attribute %in% c(cont.vals)] = "Group2"
#attribute = factor(attribute2)

#colnames(p)=coldata
#head(p)

#mm = model.matrix(~ attribute - 1)
#colnames(mm) = make.names(colnames(mm))
#colnames(mm)
#contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
#fit = lmFit(p, mm)
#fit = contrasts.fit(fit, contrast)
#fit = eBayes(fit)
#adj.method="fdr" 
#thr.p.value=0.05
#thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)

#if(nrow(tt) > 0){
#  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
#}

#print(thr.p.value)
#head(tt) 
#ph_pt=tt
#genes=rownames(ph_pt)
#ph_pt_heatm=p[rownames(p)%in%genes,]
#colnames(ph_pt_heatm)=cols

#ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
#pheatmap(ph_pt_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Differentially expressed in PHvsPT")


#ph_pt_heatm <- na.omit(ph_pt_heatm)#removes rows with NAs
##########try out pheatmaps with annotations
#annotation=t(p.meta)
#annotationh=annotation
#annotationt=annotation

#rownames(annotationh)=str_replace(rownames(annotation),"P","PH")
#rownames(annotationt)=str_replace(rownames(annotation),"P","PT")

#annotation.h=annotationh[rownames(annotationh)%in%colnames(ph_pt_heatm)[1:34],colnames(annotationh)%in%c("Gender","Skin_phototype")]
#annotation.t=annotationt[rownames(annotationt)%in%colnames(ph_pt_heatm)[34:69],colnames(annotationt)%in%c("Gender","Skin_phototype")]
#annotation.pso.h.t=rbind(annotation.h, annotation.t)

#Sample=rownames(annotation.pso.h.t)
#Sample=str_replace(Sample,"0.*","")
#annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

#annotation.pso.h.t=annotation.pso.h.t[order(annotation.pso.h.t[,3]),]#order by second colum
#target=rownames( annotation.pso.h.t)#vector according to which we will order ph_pt_heatm
#ph_pt_heatm=ph_pt_heatm[,match(target, colnames(ph_pt_heatm))]#order columns by vector

#pheatmap(ph_pt_heatm, cluster_rows = T, cluster_cols = F, scale = "row",annotation = annotation.pso.h.t)


###Vitiliigo only standardization
#v=cbind(vh,vt)
#v.st=scale_rows(v)
#v=v.st
#rows=rownames(v)
#cols=colnames(v)
#colnames(v)=cols
#rownames(v)=rows
#attribute=colnames(v)

#attribute=str_replace(attribute,"0.*","")
#coldata=colnames(v)
#coldata=str_replace(coldata,"0.*","")
#attribute = factor(attribute)
#attribute2 = rep("Muu", length(attribute))
## Case-control variables and corresponding group values 
#case.vals=c("VH")
#cont.vals =c("VT")
#attribute2[attribute %in% c(case.vals)] = "Group1"
#attribute2[attribute %in% c(cont.vals)] = "Group2"
#attribute = factor(attribute2)
#colnames(v)=coldata
#head(v)
#mm = model.matrix(~ attribute - 1)
#colnames(mm) = make.names(colnames(mm))
#colnames(mm)
#contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
#fit = lmFit(v, mm)
#fit = contrasts.fit(fit, contrast)
#fit = eBayes(fit)
#adj.method="fdr" 
#thr.p.value=0.05
#thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)

#if(nrow(tt) > 0){
  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
}

#print(thr.p.value)
#head(tt) 
#vh_vt=tt
#genes=rownames(vh_vt)
#vh_vt_heatm=v[rownames(v)%in%genes,]
#pheatmap(vh_vt_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Vitiligo.Differentially expressed in PHvsPT")

#no diff expr genes were found for vitiligo neither in VT VH standartized nor when all the data was standartized together

###########################Put Controlls, Pso, Vit in one matrix##########
data=cbind(ctrls,ph,pt,vh,vt)
#data.st=scale(data)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "row", main="All data after standardization")
summary(data)
#dim(data.st) #[1]  42 126
#data=data.st

#data is not normaly-distributed, therefore it will  b elog-transformed
data.log=log(data)
library(mvnormtest)
#multivariate shapiro normality test
mshapiro.test((data))



#### PH vs PT
p.test=data[,26:94] #dim [1] 42 69
##emove cols and rows with >50%NAs and compare the results of diffexp analysis
# remove rows with more than 50 NAs
numNAs <- apply(p.test, 1, function(z) sum(is.na(z)))
p.test.new=p.test[!(numNAs > 0.5*length(colnames(p.test))),]
pheatmap(p.test.new, cluster_rows = F, cluster_cols = F, scale = "row", main="data")
dim(p.test.new)#[1] 35 69
#remove cols with more than 50% NAs
p.test.new1=t(p.test.new)
numNAs <- apply(p.test.new1, 1, function(z) sum(is.na(z)))
p.test.new2=p.test.new1[!(numNAs >0.5*length(colnames(p.test.new1))),]
p.test.new3=t(p.test.new2)
pheatmap(p.test.new3, cluster_rows = F, cluster_cols = F,scale = "row",main="All data scaled, contains NA")

p.test=p.test.new3 #[1] 35 63
##########
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
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)

print(thr.p.value)
tt 
ph_pt=tt
genes=rownames(ph_pt)
ph_pt_heatm=p[rownames(p)%in%genes,]
colnames(ph_pt_heatm)=cols
ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(ph_pt_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis.Differentially expressed in PHvsPT")
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
#write.table(up,file="PHvsPTup_filt.txt",sep="\t",row.names=F,quote=FALSE)
#write.table(down,file="PHvsPTdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PhvsPt=Gene
##### Vh vs VT using all the data
v.test=data[,95:126]
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
head(tt) 
vh_vt=tt
genes=rownames(vh_vt)
vh_vt_heatm=v[rownames(v)%in%genes,]
colnames(vh_vt_heatm)=cols
vh_vt_heatm=vh_vt_heatm[ , ! apply( vh_vt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(vh_vt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed genes in VHvsVT")
####no diff expressed genes

########### PH vs C
t=data[,1:59]
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
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis. Differentially expressed in PHvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PHvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PHvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

PHvsC=Gene
############# PT vs C

t=data[, c(1:25,60:94)]
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
pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriasis.Differentially expressed in PTvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PTvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PTvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
PTvsC=Gene

############# VH vs C
t=data[, c(1:25,95:110)]
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
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed in VHvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="VHvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="VHvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)

VHvsC=Gene
##### VT vs C
t=data[, c(1:25,111:126)]
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
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Vitiligo. Differentially expressed in VTvsC")

tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="VTvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="VTvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
VTvsC=Gene
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


