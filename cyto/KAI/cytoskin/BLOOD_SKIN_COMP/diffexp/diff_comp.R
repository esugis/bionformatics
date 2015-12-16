
#In this script I compare the results of diff expression analysis in Skin and blood samples
#The following tasks are performed:
#1 Find overlapping genes among genes that separate PL and PNL samples and PL and C samples.
library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/BLOOD_SKIN_COMP/diffexp")
#read in the files from previous analysis. theses files contain differentially expressed genes and  some other infor including p-values 

PLvsPNL_up=read.table(file="PSLvsPSNLup_filt.txt", sep="\t", header=T)
pl_pnl_up=as.character(as.vector(PLvsPNL_up[,1]))
PLvsPNL_down=read.table(file="PSLvsPSNLdown_filt.txt", sep="\t", header=T)
pl_pnl_down=as.character(as.vector(PLvsPNL_down[,1]))
pl_pnl=c(pl_pnl_up, pl_pnl_down) #diff genes' names in PLvs PNL
length(pl_pnl)#24

PNLvsC_up=read.table(file="PSNLvsCup_filt.txt", sep="\t", header=T)
pnl_c_up=as.character(as.vector(PNLvsC_up[,1]))
PNLvsC_down=read.table(file="PSNLvsCdown_filt.txt", sep="\t", header=T)   
pnl_c_down=as.character(as.vector(PNLvsC_down[,1]))
pnl_c=c(pnl_c_up, pnl_c_down) #diff genes' names in PNLvsC
length(pnl_c)#8

lm=pnl_c[pnl_c%in%pl_pnl]# "IL1F6"   "CTLA4"   "CXCL10"  "IFNG"    "IL22RA1" "CCL20"  disease markers NL
lm1=pl_c[pl_c%in%pl_pnl]#
# [1] "PI3"     "S100A9"  "S100A8"  "IL1F6"   "LCN2"    "CCL20"   "CXCL10"  "IFNG"    "CXCL1"   "CTLA4"   "CXCL8"  
#[12] "IL1RN"   "CCL2"    "CXCL2"   "FOXP3"   "IL1b"    "AIM2"    "IFIH1"   "IL22RA1" "CCL27" 
lm1[!lm1%in%lm]
#[1] "PI3"    "S100A9" "S100A8" "LCN2"   "CXCL1"  "CXCL8"  "IL1RN"  "CCL2"   "CXCL2"  "FOXP3"  "IL1b"   "AIM2"  
#[13] "IFIH1"  "CCL27" 

PLvsC_up=read.table(file="PSLvsCup_filt.txt", sep="\t", header=T)
pl_c_up=as.character(as.vector(PLvsC_up[,1]))
PLvsC_down=read.table(file="PSLvsCdown_filt.txt", sep="\t", header=T)   
pl_c_down=as.character(as.vector(PLvsC_down[,1]))
pl_c=c(pl_c_up, pl_c_down) #diff genes' names in PNLvsC
length(pl_c)#25

#find overlap between pl_pnl and pnl_c
ol_pl_c_pnl_c=pl_c[pl_c%in%pnl_c]
length(ol_pl_c_pnl_c)#8
#"IL1F6"   "CCL20"   "CXCL10"  "IFNG"  "CTLA4"   "IL22RA1" "EOMES"   "DEFB1" 


#find overlap between pl_pnl and pnl_c
ol_pl_pnl_pnl_c=pl_pnl[pl_pnl%in%pnl_c]
length(ol_pl_pnl_pnl_c)#6
#[1] "IL1F6"   "CCL20"   "CXCL10"  "IFNG"    "CTLA4"   "IL22RA1"
save(ol_pl_pnl_pnl_c, file="ol_pl_pnl_pnl_c.RData")

#2 Find  differentially expressed genes  in PL skin samples  and psoriasis stimulated(SEB, PMA) and unstimulated blood samples.
#load skin and blood data. Select psoriasis part of the data/.
load(file="~/Documents/KAI/cytoskin/BLOOD_SKIN_COMP/data.skin.Rdata")
#data.skin
load(file="~/Documents/KAI/cytoskin/BLOOD/All_data_log.RData")
data.blood=data.log
summary(data.skin)
summary(data.blood)
#select PSO and Controlls part
data_BS_pso=data=cbind(data.skin[rownames(data.skin)%in%rownames(data.blood),], data.blood[rownames(data.blood)%in%rownames(data.skin),])
tdata_BS_pso=t(data_BS_pso)
data_BS_pso=data_BS_pso[,c(1:93,126:302)]
tdata_BS_pso=t(data_BS_pso)
dim(data_BS_pso)#[1]  11 270
#PL_SEB
p=data_BS_pso[,c(25:58,236:270)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PL.*","PL")
attribute=str_replace(attribute,".*_SEB","P_SEB")
coldata=colnames(p)
coldata=str_replace(coldata,"PL.*","PL")
coldata=str_replace(coldata,".*_SEB","P_SEB")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL")
cont.vals =c("P_SEB")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pl_pseb=tt
genes_pl_seb=rownames(pl_pseb)
save(pl_pseb,genes_pl_seb, file="pl_pseb.RData") #save the results of diff expressed analysis and gene names


#PL_PMA
p=data_BS_pso[,c(25:58,201:235)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PL.*","PL")
attribute=str_replace(attribute,".*_PMA","P_PMA")
coldata=colnames(p)
coldata=str_replace(coldata,"PL.*","PL")
coldata=str_replace(coldata,".*_PMA","P_PMA")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL")
cont.vals =c("P_PMA")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pl_ppma=tt
genes_pl_pma=rownames(pl_ppma)
save(pl_ppma,genes_pl_pma, file="pl_ppma.RData") #save the results of diff expressed analysis and gene names


#PL_72h
p=data_BS_pso[,c(25:58,166:200)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PL.*","PL")
attribute=str_replace(attribute,".*_72h","P_72h")
coldata=colnames(p)
coldata=str_replace(coldata,"PL.*","PL")
coldata=str_replace(coldata,".*_72h","P_72h")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL")
cont.vals =c("P_72h")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pl_p72h=tt
genes_pl_72h=rownames(pl_p72h)
save(pl_p72h,genes_pl_72h, file="pl_p72h.RData") #save the results of diff expressed analysis and gene names
#overlap of tose genes
genes_pl_seb_pma=genes_pl_seb[genes_pl_seb%in%genes_pl_pma]
genes_pl_seb_pma_72h= genes_pl_seb_pma[genes_pl_seb_pma%in%genes_pl_72h]
save(genes_pl_seb_pma_72h, file="genes_pl_seb_pma_72h.RData")

#3. Find  differentially expressed genes in psoriasis non lesional (PNL) skin samples  and psoriasis stimulated (SEB, PMA) and unstimulated blood samples.  
#PNL_SEB
p=data_BS_pso[,c(59:93,236:270)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PNL.*","PNL")
attribute=str_replace(attribute,".*_SEB","P_SEB")
coldata=colnames(p)
coldata=str_replace(coldata,"PNL.*","PNL")
coldata=str_replace(coldata,".*_SEB","P_SEB")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PNL")
cont.vals =c("P_SEB")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pnl_pseb=tt
genes_pnl_seb=rownames(pnl_pseb)
save(pnl_pseb,genes_pnl_seb, file="pnl_pseb.RData") #save the results of diff expressed analysis and gene names


#PNL_PMA
p=data_BS_pso[,c(59:93,201:235)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PNL.*","PNL")
attribute=str_replace(attribute,".*_PMA","P_PMA")
coldata=colnames(p)
coldata=str_replace(coldata,"PNL.*","PNL")
coldata=str_replace(coldata,".*_PMA","P_PMA")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PNL")
cont.vals =c("P_PMA")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pnl_ppma=tt
genes_pnl_pma=rownames(pnl_ppma)
save(pnl_ppma,genes_pnl_pma, file="pnl_ppma.RData") #save the results of diff expressed analysis and gene names


#PNL_72h
p=data_BS_pso[,c(59:93,166:200)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,"PNL.*","PNL")
attribute=str_replace(attribute,".*_72h","P_72h")
coldata=colnames(p)
coldata=str_replace(coldata,"PNL.*","PNL")
coldata=str_replace(coldata,".*_72h","P_72h")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PNL")
cont.vals =c("P_72h")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
pnl_p72h=tt
genes_pnl_72h=rownames(pnl_p72h)
save(pnl_p72h,genes_pnl_72h, file="pnl_p72h.RData") #save the results of diff expressed analysis and gene names
#overlap of tose genes
genes_pnl_seb_pma=genes_pnl_seb[genes_pnl_seb%in%genes_pnl_pma]
genes_pnl_seb_pma_72h= genes_pnl_seb_pma[genes_pnl_seb_pma%in%genes_pnl_72h]
save(genes_pnl_seb_pma_72h, file="genes_pnl_seb_pma_72h.RData")

#4. Find differentially expressed genes in controls (C) skin samples and controls stimulated (SEB, PMA) and unstimulated blood samples.
#C_SEB
p=data_BS_pso[,c(1:24,142:165)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_SEB","C_SEB")
attribute=str_replace(attribute,"C0.*","C")
coldata=colnames(p)
coldata=str_replace(coldata,".*_SEB","C_SEB")
coldata=str_replace(coldata,"C0.*","C")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("C")
cont.vals =c("C_SEB")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
c_cseb=tt
genes_c_seb=rownames(c_cseb)
save(c_cseb,genes_c_seb, file="c_cseb.RData") #save the results of diff expressed analysis and gene names

#C_PMA
p=data_BS_pso[,c(1:24,118:141)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_PMA","C_PMA")
attribute=str_replace(attribute,"C0.*","C")
coldata=colnames(p)
coldata=str_replace(coldata,".*_PMA","C_PMA")
coldata=str_replace(coldata,"C0.*","C")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("C")
cont.vals =c("C_PMA")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
c_cpma=tt
genes_c_pma=rownames(c_cpma)
save(c_cpma,genes_c_pma, file="c_cpma.RData") #save the results of diff expressed analysis and gene names

#C_C72h
p=data_BS_pso[,c(1:24,94:117)]
cols=colnames(p)
colnames(p)=cols
attribute=colnames(p)
attribute=str_replace(attribute,".*_72h","C_72h")
attribute=str_replace(attribute,"C0.*","C")
coldata=colnames(p)
coldata=str_replace(coldata,".*_72h","C_72h")
coldata=str_replace(coldata,"C0.*","C")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("C")
cont.vals =c("C_72h")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
#head(p)
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
c_c72h=tt
genes_c_72h=rownames(c_c72h)
save(c_c72h,genes_c_72h, file="c_c72h.RData") #save the results of diff expressed analysis and gene names


#overlap of tose genes
genes_c_seb_pma=genes_c_seb[genes_c_seb%in%genes_c_pma]
genes_c_seb_pma_72h= genes_c_seb_pma[genes_c_seb_pma%in%genes_c_72h]
save(genes_c_seb_pma_72h, file="genes_c_seb_pma_72h.RData")


#Compare overlap genes(4) with the gene lists from (2) and (3).
#genes_c_seb_pma_72h
#genes_pnl_seb_pma_72h
#genes_pl_seb_pma_72h

pl_pnl_seb_pma_72h_ol=genes_pnl_seb_pma_72h[genes_pnl_seb_pma_72h%in%genes_pl_seb_pma_72h]
#"IFIH1" "IL17F" "AIM2"  "IL10" 

pl_pnl_seb_pma_72h_ol_minusC=pl_pnl_seb_pma_72h_ol[!pl_pnl_seb_pma_72h_ol%in%genes_c_seb_pma_72h]

# compare PL vs PNL PL vs C
pl_pnl_c_minus_C=ol_pl_pnl_pnl_c[!ol_pl_pnl_pnl_c%in%genes_c_seb_pma_72h]




