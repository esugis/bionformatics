setwd("~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/anova")
#Load annotations for controls and patients
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/DDCTBloodB2mgF.RData")

###############Differential expression
library("limma")
library("pheatmap")
data=t(DDCTBloodB2mgF)
genenames=colnames(DDCTBloodB2mgF)
data=data.frame(data[complete.cases(data),])
data=sapply(data,function(x) as.numeric(as.character(x)))
data=data*(-1)
head(data)
rownames(data)=genenames
pheatmap(data_tmp_or, cluster_rows = F, cluster_cols = F,scale ="row", main="Gene expression in blood samples (B2mg) in controls, psoriasis (72h, PMA, SEB)")

##diff expressed selected genes Pus vs Cus,Pseb vs Cseb, Ppma vs Cpma


data=rbind(data, colnames(data))
data=rbind(data, colnames(data))
data=rbind(data, colnames(data))
data_tmp=t(data)
colnames(data_tmp)[c(15,16,17)]=c("Group","ID","Stimulation")
library(stringr)
data_tmp=data.frame(data_tmp)
#extract sample group
data_tmp$Group=str_replace(data_tmp$Group,"C.*","C")
data_tmp$Group=str_replace(data_tmp$Group,"P.*","P")
#extract sample number ID
data_tmp$ID=str_replace(data_tmp$ID,"C","")
data_tmp$ID=str_replace(data_tmp$ID,"P0","")#some of the psoriasis numeric IDs have 0 in the beginning
data_tmp$ID=str_replace(data_tmp$ID,"P","")
data_tmp$ID=str_replace(data_tmp$ID,"_.*","")
data_tmp$ID=as.numeric(as.character(data_tmp$ID))
#extract sample stimulation
data_tmp$Stimulation=str_replace(data_tmp$Stimulation,"C.*_","")
data_tmp$Stimulation=str_replace(data_tmp$Stimulation,"P.*_","")

#order based on Group, Treatment, ID
data_tmp=data_tmp[order(data_tmp$Group,data_tmp$Stimulation, data_tmp$ID),]
data_tmp_or=sapply(data_tmp,function(x) as.numeric(as.character(x)))
data_tmp_or=data_tmp_or[,1:14]
rownames(data_tmp_or)=rownames(data_tmp)
data_tmp_or=t(data_tmp_or)

#PpmavsCpma
data_PpmavsCpma=data_tmp_or[,c(24:46,106:138)]
library(limma)
p=data_PpmavsCpma
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"C.*","C")
attribute=str_replace(attribute,"P.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C.*","C")
coldata=str_replace(coldata,"P.*","P")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
#contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
#
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
PpmavsCpma_tt =tt 
PpmavsCpma=rownames(tt)
save(PpmavsCpma,PpmavsCpma_tt,file="PpmavsCpma.RData")

#PsebvsCseb
data_PsebvsCseb=data_tmp_or[,c(47:70,139:173)]
library(limma)
p=data_PsebvsCseb
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"C.*","C")
attribute=str_replace(attribute,"P.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C.*","")
coldata=str_replace(coldata,"P.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
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
tt

PsebvsCseb_tt =tt 
PsebvsCseb=rownames(tt)
save(PsebvsCseb,PsebvsCseb_tt,file="PsebvsCseb.RData")


#PunvsCun

data_PsebvsCseb=data_tmp_or[,c(1:23,71:105)]
p=data_PsebvsCseb
cols=colnames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"C.*","C")
attribute=str_replace(attribute,"P.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C.*","C")
coldata=str_replace(coldata,"P.*","P")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
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
tt

PunvsCun_tt =tt 
PunvsCun=rownames(tt)
save(PunvsCun,PunvsCun_tt,file="PunvsCun.RData")

selected_genes_blood=unique(c(PunvsCun,PsebvsCseb,PpmavsCpma))
######create separate annotations for each group
annotation_un=annotation
annotation_seb=annotation
annotation_pma=annotation
rownames(annotation_un)=paste(rownames(annotation_un), "72h", sep = "_")
rownames(annotation_seb)=paste(rownames(annotation_seb), "SEB", sep = "_")
rownames(annotation_pma)=paste(rownames(annotation_pma), "PMA", sep = "_")

annotation_un=annotation_un[1:35,]
annotation_seb=annotation_seb[1:35,]
annotation_pma=annotation_pma[1:35,]
p=data_tmp_or[,71:173]
dim(p) #14,103
################################# Density plots of gene expression in selected genes in PL and PNL  ##############
annotation.pso=rbind(annotation_un,annotation_seb,annotation_pma)
rownames(annotation.pso)=str_replace(rownames(annotation.pso),"P0","P")
rownames(annotation.pso)=str_replace(rownames(annotation.pso),"P0","P")
colnames(p)=str_replace(colnames(p),"P0","P")
annotation.pso=as.data.frame(annotation.pso)
Sample=rownames(annotation.pso)
Sample=str_replace(Sample,".*_72h","72h")
Sample=str_replace(Sample,".*_SEB","SEB")
Sample=str_replace(Sample,".*_PMA","PMA")
annotation.pso=data.frame(cbind(Sample,annotation.pso))
annotation.do=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Sample","Gender","Age_of_disease_onset")]
annotation.do$Age_of_disease_onset=as.numeric(as.character(annotation.do$Age_of_disease_onset))
annotation.do$Gender=as.character(annotation.do$Gender)
str(annotation.do)
annotation.do[annotation.do[,3]>=40, 3]<-"late"
annotation.do[annotation.do[,3]<40,3]<-"early"
annotation.do[annotation.do[,3]==5,3]<-"early"
annotation.do=annotation.do[order(annotation.do[,3]),]

annotation.do=cbind(annotation.do, interactions=interaction(annotation.do$Sample, annotation.do$Age_of_disease_onset))
ann=cbind(annotation.do, sample=rownames(annotation.do))

data=t(DDCTBloodB2mgF)
genenames=colnames(DDCTBloodB2mgF)
data=data.frame(data[complete.cases(data),])
data=sapply(data,function(x) as.numeric(as.character(x)))
data=data*(-1)
rownames(data)=genenames
head(data)

dt=t(data)
dt=dt[71:173,]
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[15]="sample"
dt[,15]=as.character(dt[,15])
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"

tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
save(dt,file="dt.RData")
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed blood data after imputation of non-detects.")
ggsave("SEB_PMA_72h_distribution_plot.png", width=10, height=10)


save(selected_genes_blood,dt,data_tmp_or,annotation.pso,file="data_for_anova.RData")
