#This script reads in the files computed in DiffExpPVC_removedC024.R 
#Script
#-Reads the files,
#-Proceeses the data(filtering based on 2 or more of the parralell experiments are undefined, 
#refenece gene is not stable or had very low expression,SEM should be < 0,2)
#-Assembles data set from the filtered data.
#-Scales the data(convert to Z scores)
#-FIlteres out rows and cols with >50%missing
#-Imputes missing values
#-Does K-means clustering and finds GO annotations for each cluster.
#-Does herarchical clustering and GO annotations for selected groups.
# Additionally using pvclust() finds the clusters that are strongly supported by the data.
setwd("~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/clustering")
library(stringr)
#controlls
load(file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/datacontrol_QR.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/datapsoriasis_QR.RData")     
###Vitiliigo#
load(file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/datavitiligo.RData")
###PCA plot of the data

ctrls=DataControllsBlood_QR
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:16]
ctrls <- apply(ctrls[,2:16],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBlood_QR)# psoriasis
rows=rownames(pso[2:16,])
cols=as.character(pso[1,])
ps <- apply(pso[2:16,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiligoBlood_QR)# vitiliigo
rows=rownames(vit[2:16,])
cols=as.character(vit[1,])
vi <- apply(vit[2:16,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[16,]=gsub("C.*_","",ctrlstst[16,])
ctrlstst=ctrlstst[,order(ctrlstst[16,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[17,]=gsub("C","",ctrlstst[17,])
ctrlstst[17,]=gsub("_.*","",ctrlstst[17,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,17]=as.numeric(as.character(ctrlstst[,17]))
C_72h=ctrlstst[ctrlstst[,16]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,17]),]
C_PMA=ctrlstst[ctrlstst[,16]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,17]),]
C_SEB=ctrlstst[ctrlstst[,16]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,17]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:15,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
ctrls.log=log2(ctrls.data)
rownames(ctrls.log)=genes
rownames(ctrls.log)
pheatmap(ctrls.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in controls (72h, PMA, SEB)")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[16,]=gsub("P.*_","",pstst[16,])
pstst=pstst[,order(pstst[16,])]
pstst=rbind(pstst, colnames(pstst))
pstst[17,]=gsub("P","",pstst[17,])
pstst[17,]=gsub("_.*","",pstst[17,])
pstst=as.data.frame(t(pstst))
pstst[,17]=as.numeric(as.character(pstst[,17]))
P_72h=pstst[pstst[,16]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,17]),]
P_PMA=pstst[pstst[,16]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,17]),]
P_SEB=pstst[pstst[,16]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,17]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:15,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))
ps.log=log2(ps.data)
rownames(ps.log)=genes
rownames(ps.log)
pheatmap(ps.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of psoriasis patients (72h, PMA, SEB)")

#VIT
vitst=rbind(vi, colnames(vi))
vitst[16,]=gsub("V.*_","",vitst[16,])
vitst=vitst[,order(vitst[16,])]
vitst=rbind(vitst, colnames(vitst))
vitst[17,]=gsub("V","",vitst[17,])
vitst[17,]=gsub("_.*","",vitst[17,])
vitst=as.data.frame(t(vitst))
vitst[,17]=as.numeric(as.character(vitst[,17]))
V_72h=vitst[vitst[,16]%in%c("72h"),]
V_72h=V_72h[order(V_72h[,17]),]
V_PMA=vitst[vitst[,16]%in%c("PMA"),]
V_PMA=V_PMA[order(V_PMA[,17]),]
V_SEB=vitst[vitst[,16]%in%c("SEB"),]
V_SEB=V_SEB[order(V_SEB[,17]),]
vitst_ordered=rbind(V_72h,V_PMA,V_SEB)
vi.data=t(vitst_ordered)

vi.data=vi.data[1:15,]
genes=rownames(vi.data)
vi.data=as.matrix(apply(vi.data,2, as.numeric))
vi.log=log2(vi.data)
rownames(vi.log)=genes
rownames(vi.log)
pheatmap(vi.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of vitiligo patients (72h, PMA, SEB)")



#assemble one big ds conatining all groups in controls, psoriasis and vitiligo
data.log.FR=cbind(ctrls.log,ps.log, vi.log)
pheatmap(data.log.FR, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples ofcontrols, psoriasis and vitiligo patients (72h, PMA, SEB)")
save(data.log.FR, file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/All_data_log.RData")

#select only psoriasis and controlls
data.pso.FR=cbind(ctrls.log,ps.log)
save(data.pso.FR, file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/PSO_data_log.RData")

#####analysis of interleukins###########
il=c("IL17A","IL17F", "IL22","IL4", "IL10" )#select interleukins from the data
data.pso=data.pso.FR
#analysis of PMA group
il_pso=data.pso[rownames(data.pso)%in%il,c(25:48,108:142)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression of interleukins in blood samples of PMA/lono stimulated controls and psoriasis ")

##diffexp PSO vs CTRLS
library(limma)
p=il_pso
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
p_c=tt
genes=rownames(p_c)
#Result-one gene
#logFC  AveExpr        t    P.Value adj.P.Val         B
#IL17F 1.991755 1.839656 3.263474 0.00190208 0.0095104 -1.435767

#analysis of SEB group
il_pso=data.pso[rownames(data.pso)%in%il,c(49:72,143:177)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression of interleukins in blood samples of SEB stimulated controls and psoriasis ")

##diffexp PSO vs CTRLS
library(limma)
p=il_pso
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
p_c=tt
genes=rownames(p_c)
#         logFC  AveExpr        t      P.Value   adj.P.Val        B
#IL17F 1.281576 5.793522 3.826515 0.0002655572 0.001327786 0.224565
> 
#####Diff expression analysis of other groups.

il=c("IL17A","IL17F", "IL22","IL4", "IL10")#select interleukins from the data

#analysis of PMA group
g_pso=data.pso[!rownames(data.pso)%in%il,c(25:48,108:142)]
dim(g_pso)
g_pso_t=t(g_pso)
pheatmap(g_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of PMA/lono stimulated controls and psoriasis groups")

##diffexp PSO PMA vs CTRLS
library(limma)
p=g_pso
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
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
p_c=tt
genes=rownames(p_c)
#Result
#logFC   AveExpr        t     P.Value  adj.P.Val         B
#PYCARD 0.8600192 0.7581680 3.288985 0.001772696 0.01751766 -1.365931
#TIGIT2 0.6436355 0.9515226 3.037821 0.003503533 0.01751766 -2.017593

#analysis of SEB group
il_pso=data.pso[!rownames(data.pso)%in%il,c(49:72,143:177)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of SEB stimulated controls and psoriasis groups ")

##diffexp PSO vs CTRLS
library(limma)
p=il_pso
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
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
#0
#analysis of unstimulated group
il_pso=data.pso[!rownames(data.pso)%in%il,c(1:24,73:107)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of unstimulated controls and psoriasis groups ")

##diffexp PSO vs CTRLS
library(limma)
p=il_pso
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
#0


####Diffexp PSO vs PSO
#Pso PMA vs Pso SEB

il_pso=data.pso[!rownames(data.pso)%in%il,c(108:177)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of PMA/lono and SEB stimulated psoriasis groups ")

library(limma)
p=il_pso
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_PMA","PMA")
attribute=str_replace(attribute,".*_SEB","SEB")
coldata=colnames(p)
coldata=str_replace(coldata,".*_PMA","PMA")
coldata=str_replace(coldata,".*_SEB","SEB")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PMA")
cont.vals =c("SEB")
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
#logFC    AveExpr          t      P.Value    adj.P.Val          B
#PYCARD -5.1169147  3.9694781 -27.090240 5.287543e-38 5.287543e-37 76.1725204
#OAS2   -2.2731229  0.5887565 -16.209499 1.737187e-25 8.685933e-25 47.6435957
#FOXP3  -3.3469495  2.6592965  -8.858070 3.613463e-13 1.204488e-12 19.3101110
#IFIH1  -1.8497376  0.2079006  -8.750429 1.635977e-12 4.089944e-12 17.9257149
#CTLA4  -1.9008611  2.4452822  -6.358774 1.602819e-08 3.205638e-08  8.6939931
#PI3     6.1482077 -3.7528815   9.504657 4.324344e-07 7.207240e-07  6.7308700
#TIGIT2 -0.9506678  1.7051872  -5.429030 7.167420e-07 1.023917e-06  4.9571268
#EOMES  -1.3048057  1.9297721  -4.739379 1.039755e-05 1.299694e-05  2.3482101
#INFG    1.1062561 -0.5330277   4.103654 1.117884e-04 1.242093e-04  0.1028141


#Pso PMA vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(73:142)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of PMA/lono stimulated and unstimulated psoriasis groups ")

library(limma)
p=il_pso
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_PMA","PMA")
attribute=str_replace(attribute,".*_72h","72h")
coldata=colnames(p)
coldata=str_replace(coldata,".*_PMA","PMA")
coldata=str_replace(coldata,".*_72h","72h")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PMA")
cont.vals =c("72h")
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

#logFC    AveExpr          t      P.Value    adj.P.Val          B
#INFG   13.2571757 -6.4697991  45.910609 7.812708e-52 7.812708e-51 105.734922
#PYCARD -6.2578653  4.6033395 -31.842601 2.137772e-43 1.068886e-42  88.019560
#CTLA4   4.9687560 -1.0905502  16.752909 4.544025e-27 1.514675e-26  51.243761
#AIM2    3.2421251 -1.2303603  13.072351 4.984670e-21 1.246167e-20  37.403642
#TIGIT2  0.6872963  0.8621174   3.538216 6.951784e-04 1.390357e-03  -1.582502
#PI3     2.7251421 -1.8466957   3.540252 3.251163e-03 4.644519e-03  -1.718287
#FOXP3   1.2531993  0.3006500   3.266381 1.652008e-03 2.753347e-03  -2.387692
#IFIH1  -0.6335594 -0.3742940  -2.508471 1.447848e-02 1.809810e-02  -4.338068
#OAS2   -0.3929808 -0.3455354  -2.364156 2.065961e-02 2.295512e-02  -4.705106

#Pso SEB vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(73:107, 143:177)]
dim(il_pso)
il_pso_t=t(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of SEB stimulated and unstimulated psoriasis groups ")

library(limma)
p=il_pso
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_SEB","SEB")
attribute=str_replace(attribute,".*_72h","72h")
coldata=colnames(p)
coldata=str_replace(coldata,".*_SEB","SEB")
coldata=str_replace(coldata,".*_72h","72h")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("SEB")
cont.vals =c("72h")
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
#logFC     AveExpr         t      P.Value    adj.P.Val          B
#INFG   12.150920 -6.74331099 47.245445 5.047355e-56 5.047355e-55 115.978554
#CTLA4   6.869617 -0.06704967 24.831655 2.153002e-39 1.076501e-38  79.375856
#FOXP3   4.600149  2.01681140 19.419696 7.606103e-32 2.535368e-31  62.143014
#AIM2    2.829977 -1.40926320 12.478931 2.198123e-20 5.495308e-20  35.824080
#OAS2    1.880142  0.75759778 11.355700 3.265775e-18 6.531550e-18  30.843807
#TIGIT2  1.637964  1.34755860  8.618516 4.865705e-13 8.109509e-13  18.946798
#PYCARD -1.140951  6.81413770 -8.114863 4.755903e-12 6.794148e-12  16.680743
#EOMES   1.338170  1.89390169  5.567647 3.329100e-07 4.161376e-07   5.659735
#IFIH1   1.216178  0.51045009  5.223980 1.669105e-06 1.854561e-06   4.174200
#PI3    -3.423066 -5.35995726 -5.096094 8.709805e-05 8.709805e-05   1.599704
