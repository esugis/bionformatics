
#Script finds differentially expressed interleuking in the selected groups
#and other differentially expressed genes in the conditions specified below for each case.

setwd("~/Documents/KAI/cytoskin/BLOOD_BACT")
library(stringr)
#controlls
load(file="~/Documents/KAI/cytoskin/BLOOD_BACT/DataControllsBloodBact_QR.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/BLOOD_BACT/DataPsoriasisBloodBact_QR.RData")     
###Vitiliigo#
#load(file="~/Documents/KAI/cytoskin/BLOOD/datavitiligo_QR.RData")
###PCA plot of the data

ctrls=DataControllsBloodBact_QR
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:15]
ctrls <- apply(ctrls[,2:15],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBloodBact_QR)# psoriasis
rows=rownames(pso[2:15,])
cols=as.character(pso[1,])
ps <- apply(pso[2:15,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
#vit=t(DataVitiligoBlood_QR)# vitiliigo
#rows=rownames(vit[2:16,])
#cols=as.character(vit[1,])
#vi <- apply(vit[2:16,],2,as.numeric)
#colnames(vi)=cols
#rownames(vi)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[15,]=gsub("C.*_","",ctrlstst[15,])
ctrlstst=ctrlstst[,order(ctrlstst[15,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[16,]=gsub("C","",ctrlstst[16,])
ctrlstst[16,]=gsub("_.*","",ctrlstst[16,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,16]=as.numeric(as.character(ctrlstst[,16]))
C_72h=ctrlstst[ctrlstst[,15]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,16]),]
C_PMA=ctrlstst[ctrlstst[,15]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,16]),]
C_SEB=ctrlstst[ctrlstst[,15]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,16]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:14,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
ctrls.log=log2(ctrls.data)
rownames(ctrls.log)=genes
rownames(ctrls.log)
pheatmap(ctrls.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in controls (72h, PMA, SEB)")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[15,]=gsub("P.*_","",pstst[15,])
pstst=pstst[,order(pstst[15,])]
pstst=rbind(pstst, colnames(pstst))
pstst[16,]=gsub("P","",pstst[16,])
pstst[16,]=gsub("_.*","",pstst[16,])
pstst=as.data.frame(t(pstst))
pstst[,16]=as.numeric(as.character(pstst[,16]))
P_72h=pstst[pstst[,15]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,16]),]
P_PMA=pstst[pstst[,15]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,16]),]
P_SEB=pstst[pstst[,15]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,16]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:14,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))
ps.log=log2(ps.data)
rownames(ps.log)=genes
rownames(ps.log)
pheatmap(ps.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of psoriasis patients (72h, PMA, SEB)")

#VIT
#vitst=rbind(vi, colnames(vi))
#vitst[16,]=gsub("V.*_","",vitst[16,])
#vitst=vitst[,order(vitst[16,])]
#vitst=rbind(vitst, colnames(vitst))
#vitst[17,]=gsub("V","",vitst[17,])
#vitst[17,]=gsub("_.*","",vitst[17,])
#vitst=as.data.frame(t(vitst))
#vitst[,17]=as.numeric(as.character(vitst[,17]))
#V_72h=vitst[vitst[,16]%in%c("72h"),]
#V_72h=V_72h[order(V_72h[,17]),]
#V_PMA=vitst[vitst[,16]%in%c("PMA"),]
#V_PMA=V_PMA[order(V_PMA[,17]),]
#V_SEB=vitst[vitst[,16]%in%c("SEB"),]
#V_SEB=V_SEB[order(V_SEB[,17]),]
#vitst_ordered=rbind(V_72h,V_PMA,V_SEB)
#vi.data=t(vitst_ordered)

#vi.data=vi.data[1:15,]
#genes=rownames(vi.data)
#vi.data=as.matrix(apply(vi.data,2, as.numeric))
#vi.log=log2(vi.data)
#rownames(vi.log)=genes
#rownames(vi.log)
#pheatmap(vi.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of vitiligo patients (72h, PMA, SEB)")



#assemble one big ds conatining all groups in controls, psoriasis and vitiligo
data.log=cbind(ctrls.log,ps.log)
pheatmap(data.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of controls and psoriasis patients (72h, PMA, SEB)")
#save(data.log, file="~/Documents/KAI/cytoskin/BLOOD_BACT/data_ctrl_pso_bact_log.RData")

#select only psoriasis and controlls
data.pso=cbind(ctrls.log,ps.log)
save(data.pso, file="~/Documents/KAI/cytoskin/BLOOD_BACT/PSO_data_log.RData")

#####analysis of interleukins###########
il=c("IL17A","IL17F", "IL22","IL4", "IL10" )#select interleukins from the data

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
#Result-0 genes

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
#          logFC  AveExpr        t     P.Value  adj.P.Val         B
#IL17F 0.8908039 2.832813 2.617777 0.009418849 0.04709425 -2.743399


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
#Result 0 genes

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
#Result 0 genes

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
#result 0 genes


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
#            logFC    AveExpr         t      P.Value    adj.P.Val         B
#AIM2    3.3046210 -1.8483556 17.608824 2.318457e-29 2.086611e-28 56.445044
#INFG    3.9620421 -2.1965189 15.433009 4.464210e-25 2.008895e-24 46.678299
#PYCARD -2.0130086  2.3114346 -9.270967 5.358395e-14 1.607519e-13 21.299959
#TIGIT2  1.9516004 -0.3803968  8.368772 1.449677e-12 3.261773e-12 17.939256
#EOMES   1.5621274 -0.1802580  5.437916 5.614608e-07 1.010629e-06  5.239224
#IFIH1   1.0502873 -1.8405355  4.866501 6.412841e-06 9.619261e-06  2.944139
#CTLA4   1.0191574  0.3233478  4.559776 1.804029e-05 2.319466e-05  1.864993
#OAS2    0.6154688 -1.5155683  2.823318 5.985347e-03 6.733516e-03 -3.629992

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
#            logFC     AveExpr          t      P.Value    adj.P.Val          B
#INFG   14.2243472 -7.09824649  46.406259 5.679927e-59 5.111934e-58 120.803183
#CTLA4   6.0122022 -2.18055982  22.748066 2.489877e-38 1.120445e-37  76.687841
#PYCARD -5.0331992  3.96766808 -21.067822 4.148091e-34 1.244427e-33  67.087557
#AIM2    4.2688265 -2.35440492  15.543146 8.316837e-27 1.871288e-26  50.560675
#FOXP3   2.3144430 -0.85543182   6.145711 2.502321e-08 4.504178e-08   8.447203
#TIGIT2  1.7368855 -0.28718142   6.034304 3.825991e-08 5.738987e-08   8.013412
#EOMES   1.0029236  0.09529168   3.084386 2.744601e-03 3.528773e-03  -2.730514
#OAS2    0.6218308 -1.52320917   2.334639 2.187132e-02 2.460524e-02  -4.61783


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
#            logFC    AveExpr          t      P.Value    adj.P.Val           B
#INFG   10.2623051 -9.1431714  36.677706 4.094863e-49 3.685377e-48 100.4808884
#CTLA4   4.9930448 -2.7263200  19.603013 1.923766e-32 8.656947e-32  63.4624158
#PYCARD -3.0201905  4.7306305 -16.609363 1.290137e-27 3.870412e-27  52.4587843
#FOXP3   2.7265943 -0.6493562   9.855854 1.900680e-15 4.276529e-15  24.5961380
#AIM2    0.9642055 -4.0137024   3.730930 3.532875e-04 6.359174e-04  -0.9175114
#IFIH1  -0.7678182 -2.0010899  -3.525794 7.451998e-04 1.117800e-03  -1.5343029
