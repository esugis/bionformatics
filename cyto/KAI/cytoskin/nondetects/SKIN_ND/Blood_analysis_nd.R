
library(stringr)
library("HTqPCR")#nondetects uses an objects 
library("nondetects")
#controlls
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/datacontrol_nd_bact.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/datapsoriasis_nd_bact.RData")     
###Vitiliigo#

ctrls=DataControllsBlood_nd_bact
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:15]
ctrls <- apply(ctrls[,2:15],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBlood_nd_bact)# psoriasis
rows=rownames(pso[2:15,])
cols=as.character(pso[1,])
ps <- apply(pso[2:15,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
dim(ctrlstst)#15,216
ctrlstst[15,]=gsub("C.*_","",ctrlstst[15,])
ctrlstst[15,]=gsub(".1","",ctrlstst[15,])
ctrlstst[15,]=gsub(".2","",ctrlstst[15,])
ctrlstst[15,]=gsub("h","72h",ctrlstst[15,])

ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[16,]=gsub("C","",ctrlstst[16,])
ctrlstst[16,]=gsub("_.*","",ctrlstst[16,])
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[17,]=gsub("C.*","C",ctrlstst[17,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,16]=as.numeric(as.character(ctrlstst[,16]))
ctrlstst=ctrlstst[order(ctrlstst$V17),]
colnames(ctrlstst)[16]="ID"
colnames(ctrlstst)[15]="Stimulation"
colnames(ctrlstst)[17]="Group"
ctrlstst_or=t(ctrlstst)


#PSO
pstst=rbind(ps, colnames(ps))
pstst[15,]=gsub("P.*_","",pstst[15,])
pstst[15,]=gsub(".1","",pstst[15,])
pstst[15,]=gsub(".2","",pstst[15,])
pstst[15,]=gsub("h","72h",pstst[15,])

pstst=rbind(pstst, colnames(pstst))
pstst[16,]=gsub("P","",pstst[16,])
pstst[16,]=gsub("_.*","",pstst[16,])
pstst=rbind(pstst, colnames(pstst))
pstst[17,]=gsub("P.*","P",pstst[17,])
pstst=as.data.frame(t(pstst))
pstst[,16]=as.numeric(as.character(pstst[,16]))
pstst=pstst[order(pstst$V17),]
colnames(pstst)[16]="ID"
colnames(pstst)[15]="Stimulation"
colnames(pstst)[17]="Group"
pstst_or=t(pstst)

### Combine C and PSO into one dataframe
cp_or=cbind(ctrlstst_or,pstst_or)
cp_or=data.frame(cp_or)
# Convert data to qPCRset

head(cp_or)

test=(colnames(cp_or))     
Sample=(colnames(cp_or)) 
Sample=gsub("_.*","",Sample)

head(test)
#[1] "C1_72h"   "C1_72h.1" "C1_72h.2" "C1_PMA"   "C1_PMA.1" "C1_PMA.2"

test=gsub("C.*_","",test)
head(test)
#[1] "72h"   "72h.1" "72h.2" "PMA"   "PMA.1" "PMA.2"
test=gsub("P.*_","",test)
head(test)
#[1] "72h"   "72h.1" "72h.2" "PMA"   "PMA.1" "PMA.2"
test=gsub("72h","0",test)
test=gsub("PMA","0",test)
test=gsub("SEB","0",test)
test=gsub("0.2","3",test)
test=gsub("0.1","2",test)
test=gsub("0","1",test)
head(test)
Replicate=test
cp_or=t(cp_or)
cp_or=cbind(cp_or,test)
cp_or=cbind(cp_or,Sample)

cp_or=t(cp_or)
rownames(cp_or)[18]="Replicate"
#rownames(DataControllsBlood)=as.vector(as.character(DataControllsBlood$SAMPLE))
#convert to numeric
cp_or=data.frame(cp_or)

tmp1=cp_or[1:14,]
tmp1=data.frame(tmp1)

tmp1=sapply(tmp1,function(x) as.numeric(as.character(x)))
tmp1[is.na(tmp1)]<-40
#mat=data.matrix(dat)
#dat=sapply(tmp1,function(x) as.numeric(as.character(x)))
dat=tmp1
mat=data.matrix(dat)
str(mat)
rownames(mat)=rownames(cp_or[1:14,])
#mat=t(mat)
dimnames(mat)=NULL
#raw <- new("qPCRset", exprs = mat, featureCategory = as.data.frame(array("OK",  dim=dim(mat))))
raw <- new("qPCRset", exprs = mat)
raw@phenoData@data$sampleNamesShort <-Sample
raw@phenoData@data$sampleName <-as.character(colnames(cp_or))
raw@phenoData@data$Replicate <- Replicate
tmp=t(cp_or)
Group=as.character(tmp[,17])
raw@phenoData@data$Group <- Group
Stimulation=as.character(tmp[,15])
raw@phenoData@data$Stimulation <- Stimulation
featureNames(raw) <- rownames(cp_or)[1:14]
colnames(exprs(raw))=colnames(cp_or)
#set feature category
FC=data.frame(exprs(raw))
FC[!FC==40]<-c("OK")
FC[FC==40]<-c("Undetermined")
featureCategory(raw)=FC#check what features are undetermined or needed to be indicated

#set feature type
FT=replicate(14,"Target")
#FT[15]="Endogenous Control"
featureType(raw)=FT
exprs(raw)#see expression matrix

#Imputation
#data(oncogene2013)
#oncogene2013 <- qpcrImpute(oncogene2013,groupVars=c("sampleType","treatment"))
#colnames(exprs(raw))=as.character(colnames(cp_or))
#tmp1=data.frame(tmp1)
#tmp1[is.na(tmp1)]<-40
#View(tmp1)
#tmp1=sapply(cp_or[1:15,],function(x) as.numeric(as.character(x)))
#tmp1=data.matrix(tmp1)
#raw <- new("qPCRset", exprs = tmp1, featureCategory = as.data.frame(array("OK",  dim=dim(tmp1))))
#exprs(raw)=tmp1


pData(raw)[sapply(pData(raw), is.character)] <- lapply(pData(raw)[sapply(pData(raw), is.character)], 
                                       as.factor)
#cpImp=qpcrImpute(raw,groupVars=c("sampleNamesShort","Group","Stimulation"))
#exprs(cpImp)
#cpImp_SamGrSt_bact=cpImp
#exprs(cpImp_SamGrSt_bact)
#save(cpImp_SamGrSt_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/cpImp_SamGrSt_bact.RData")
cpImp_GrSt=qpcrImpute(raw,groupVars=c("Group","Stimulation"))
exprs(cpImp_GrSt)
cpImp_GrSt_bact=cpImp_GrSt
save(cpImp_GrSt_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/cpImp_GrSt_bact.RData")
#check with other dataset
data(oncogene2013)
#oncogeneImp=qpcrImpute(oncogene2013,groupVars=c("sampleType","treatment"))
exprs(oncogeneImp)
data(sagmb2011)
exprs(sagmb2011)
#### Extract imputed expression part and
#assemble a dataframe with mean values of Ct-s and STD of Ct-s per patient/control

CT=exprs(cpImp_GrSt_bact)
#CT=exprs(cpImp_SamGrSt_bact)
head(CT)

byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}
#calculate means of 3 replicate per patient
CTMeans <- byapply(CT, 3, rowMeans)
dim(CTMeans) #14 177
meanNames=colnames(CT)
meanNames=gsub("h.2","h",meanNames)
meanNames=gsub("B.2","B",meanNames)
meanNames=gsub("A.2","A",meanNames)
meanNames=gsub("h.1","h",meanNames)
meanNames=gsub("B.1","B",meanNames)
meanNames=gsub("A.1","A",meanNames)

meanNames=unique(meanNames)
length(meanNames)
colnames(CTMeans)=meanNames

#calculate SD of 3 replicates per patient
CTStd <- byapply(CT, 3, sd)
mystd<-function(x){apply(x,1,sd)}
CTStd <- byapply(CT, 3, mystd)
dim(CTStd) #14 177
colnames(CTStd)=meanNames
#test if the values are correct
library("matrixStats")

#bind together CTMean and CTStd
CTMean_bact=t(CTMeans)
CTMean_bact=cbind(rownames(CTMean_bact), CTMean_bact)
colnames(CTMean_bact)[1]="SAMPLE"
CTStd_bact=t(CTStd)
CTStd_bact=cbind(rownames(CTStd_bact), CTStd_bact)
colnames(CTStd_bact)[1]="SAMPLE"
#combine with the reference gene expression.
#create a data frame for the reference gene in C/P samples and calibrator.
#Store calibtaror separately.
#dataframes are located in
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefStd_bact.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefMean_bact.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefCal_bact.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefMean_PSO_bact.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefStd_PSO_bact.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefCal_PSO_bact.RData")

#combine CT values from Controls and Psoriasis
CTRefMean_CP_bact=rbind(CTRefMean_bact,CTRefMean_PSO_bact )
CTBlood_bact=merge(CTMean_bact,CTRefMean_CP_bact, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
rownames(CTBlood_bact)=CTBlood_bact[,1]
rowsCTBlood_bact=rownames(CTBlood_bact)
#combine Std from Controls and Psoriasis
CTRefStd_CP_bact=rbind(CTRefStd_bact,CTRefStd_PSO_bact)
CTStdBlood_bact=merge(CTStd_bact,CTRefStd_CP_bact, by.x="SAMPLE", by.y="SAMPLE", all.x=T ) 
rownames(CTStdBlood_bact)=CTStdBlood_bact[,1]
rowsCTStdBlood_bact=rownames(CTStdBlood_bact)
#combine RefCal for C and P
CTRefCal_CP_bact=rbind(CTRefCal_bact,CTRefCal_PSO_bact)
rownames(CTRefCal_CP_bact)=CTRefCal_CP_bact[,1]
rowsCTRefCal_CP_bact=rownames(CTRefCal_CP_bact)

DDCTBlood_bact=CTBlood_bact[,2:15]#only delta-delta CTs for 14 genes
rownames(DDCTBlood_bact)=CTBlood_bact[,1]
rowsDDCTBlood_bact=rownames(DDCTBlood_bact)

SEMCTBlood_bact=CTBlood_bact[,2:15] #only SEM
rownames(SEMCTBlood_bact)=CTBlood_bact[,1]
rowsSEMCTBlood_bact=rownames(SEMCTBlood_bact)

CTBlood_bact=CTBlood_bact[,2:16] #exclude  first column sample names
CTStdBlood_bact=CTStdBlood_bact[,2:16]#exclude  first column sample names

DDCTBlood_bact=sapply(DDCTBlood_bact,function(x) as.numeric(as.character(x)))
rownames(DDCTBlood_bact)=rowsDDCTBlood_bact
SEMCTBlood_bact=sapply(SEMCTBlood_bact,function(x) as.numeric(as.character(x)))
rownames(SEMCTBlood_bact)=rowsSEMCTBlood_bact

CTBlood_bact=sapply(CTBlood_bact,function(x) as.numeric(as.character(x)))
CTStdBlood_bact=sapply(CTStdBlood_bact,function(x) as.numeric(as.character(x)))
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefCal_bact.RData")
CTRefCal_CP_bact=CTRefCal_CP_bact[,2:15]
CTRefCal_CP_bact=sapply(CTRefCal_CP_bact,function(x) as.numeric(as.character(x)))

head(CTRefCal_CP_bact)
for(i in 1:length(rowsCTStdBlood_bact)){
  print(i)
  for(j in 1:14){ #14 genes
    print(j)
    print(CTBlood_bact[i,j])
    print(CTBlood_bact[i,15])
    print(CTRefCal_bact[i,j])
    DDCTBlood_bact[i,j]=CTBlood_bact[i,j]-CTBlood_bact[i,15]-CTRefCal_CP_bact[i,j]
    SEMCTBlood_bact[i,j]=sqrt((CTStdBlood_bact[i,j])^2+(CTStdBlood_bact[i,15])^2)/(sqrt(3+3/2))
  }
}

#save results
save(DDCTBlood_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/DDCTBlood_bact.RData")
save(SEMCTBlood_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/SEMCTBlood_bact.RData")

##Filter Delta-delta CTs based on the SEM<0,2

SEMCTBlood_bact=SEMCTBlood_bact[complete.cases(SEMCTBlood_bact),] # samples were filtered out because the 
DDCTBlood_bact=DDCTBlood_bact[complete.cases(DDCTBlood_bact),]
DDCTBlood_bact_F=DDCTBlood_bact

for(i in 1:length(rownames(DDCTBlood_bact))){
  print(i)
  for(j in 1:14){ #14 genes
    print(j)
    if(SEMCTBlood_bact[i,j]<= 0.2){
    DDCTBlood_bact_F[i,j]=DDCTBlood_bact[i,j]
    } else{
      if(SEMCTBlood_bact[i,j]> 0.2){ 
        DDCTBlood_bact_F[i,j]="NA"
      }
    }
  }
}
#save delta-delta CTs after filtering based on qurtiles of reference gene CTs and SEM<0.2
save(DDCTBlood_bact_F, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/DDCTBlood_bact_F.RData")

###############Differential expression
library("limma")
library("pheatmap")
data=t(DDCTBlood_bact_F)
genenames=colnames(DDCTBlood_bact_F)
data=data.frame(data[complete.cases(data),])
data=sapply(data,function(x) as.numeric(as.character(x)))
data=data*(-1)
head(data)
rownames(data)=genenames
#pheatmap(data, cluster_rows = F, cluster_cols = F,scale ="row", main="Gene expression in blood samples (B2mg) in controls, psoriasis (72h, PMA, SEB)")
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
pheatmap(data_tmp_or, cluster_rows = F, cluster_cols = F,scale ="row", main="Gene expression(delta-delta CT) in blood samples (Bact) in controls, psoriasis (72h, PMA, SEB)")

#####analysis of interleukins###########
il=c("IL17A","IL17F", "IL22","IL4", "IL10" )#select interleukins from the data

#analysis of PMA group
il_pma=data_tmp_or[rownames(data_tmp_or)%in%il,c(25:48,108:142)]
dim(il_pma)
pheatmap(il_pma, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression (DDCT) of interleukins in blood samples (Bact) of PMA/lono stimulated controls and psoriasis ")

##diffexp PSO vs CTRLS
library(limma)
p=il_pma
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
tt
# data frame with 0 columns and 0 rows

#analysis of SEB group
il_seb=data_tmp_or[rownames(data_tmp_or)%in%il,c(49:72,143:177)]
dim(il_seb)#5 59
pheatmap(il_seb, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression (DDCT) of interleukins in blood samples (Bact) of SEB stimulated controls and psoriasis ")

##diffexp PSO vs CTRLS
library(limma)
p=il_seb
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
#data frame with 0 columns and 0 rows

########Diff expression analysis of other groups.

il=c("IL17A","IL17F", "IL22","IL4", "IL10")#select interleukins from the data

#analysis of PMA group
g_pma=data_tmp_or[!rownames(data_tmp_or)%in%il,c(25:48,108:142)]
dim(g_pma)#9 56
pheatmap(g_pma, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of PMA/lono stimulated controls and psoriasis groups")

##diffexp PSO PMA vs CTRLS
library(limma)
p=g_pma
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
tt
# data frame with 0 columns and 0 rows

#analysis of SEB group
g_seb=data_tmp_or[!rownames(data_tmp_or)%in%il,c(49:72,143:177)]
dim(g_seb)# 9,59
pheatmap(g_seb, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of SEB stimulated controls and psoriasis groups ")

##diffexp PSO vs CTRLS
library(limma)
p=g_seb
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
tt
#no diff exp genes

#analysis of unstimulated group
g_72h=data_tmp_or[!(rownames(data_tmp_or))%in%il,c(1:24,73:107)]
dim(g_72h) #9, 59
pheatmap(g_72h, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of unstimulated controls and psoriasis groups ")

##diffexp PSO vs CTRLS
library(limma)
p=g_72h
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
#no diff expressed genes

####Diffexp PSO vs PSO
#Pso PMA vs Pso SEB
data.pso=data_tmp_or
il_pso=data.pso[!rownames(data.pso)%in%il,c(108:177)]
dim(il_pso)
il_pso_t=t(il_pso)#9 70
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of PMA/lono and SEB stimulated psoriasis groups ")

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
tt
#            logFC    AveExpr          t      P.Value    adj.P.Val         B
#AIM2    3.3092586 -1.8272304  17.677108 1.233235e-29 1.109911e-28 57.060326
#INFG    4.0429273 -2.1677510  15.342713 2.485238e-24 1.118357e-23 44.997426
#PYCARD -2.0633264  2.3470229 -10.812216 8.827492e-17 2.648248e-16 27.656335
#TIGIT2  1.9462610 -0.3253182   8.266971 2.502975e-12 5.631693e-12 17.354143
#IFIH1   1.2512741 -1.8133273   6.600262 6.776356e-09 1.219744e-08  9.639860
#EOMES   1.5874306 -0.1879646   5.677696 2.083614e-07 3.125421e-07  6.147322
#CTLA4   1.0185535  0.3249517   4.625147 1.369960e-05 1.761377e-05  2.048925
#OAS2    0.6356072 -1.5198088   2.923586 4.477008e-03 5.036634e-03 -3.438517

#Pso PMA vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(73:142)]
dim(il_pso)# 9 68
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of PMA/lono stimulated and unstimulated psoriasis groups ")

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
tt
#            logFC     AveExpr          t      P.Value    adj.P.Val         B
#INFG   14.3241311 -6.32532429  42.952041 9.264314e-49 8.337883e-48 98.076578
#CTLA4   5.8342517 -2.00448039  25.447136 2.605705e-40 1.172567e-39 81.230027
#PYCARD -5.0849408  4.08829218 -24.180313 3.287923e-36 9.863770e-36 71.853584
#AIM2    4.2682190 -2.35537619  15.770918 1.177722e-26 2.649874e-26 50.244214
#FOXP3   2.6454694 -0.72941537   8.782937 2.653760e-13 4.776769e-13 19.725921
#TIGIT2  1.7726283 -0.28083575   6.187769 2.307944e-08 3.461916e-08  8.454459
#EOMES   1.0750067  0.06428492   3.345145 1.249742e-03 1.606811e-03 -2.072698
#OAS2    0.5911372 -1.49368657   2.174499 3.263276e-02 3.671185e-02 -5.030800

#Pso SEB vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(73:107, 143:177)]
dim(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples (Bact) of SEB stimulated and unstimulated psoriasis groups ")

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
tt
#            logFC    AveExpr          t      P.Value    adj.P.Val          B
#INFG   10.2812038 -8.6242437  32.507593 2.818007e-40 2.536207e-39 80.4682876
#CTLA4   4.8156982 -2.5636169  21.851034 8.965072e-35 4.034282e-34 68.7763155
#PYCARD -3.0216144  4.7321209 -16.867304 9.242041e-28 2.772612e-27 52.8021068
#FOXP3   2.7168955 -0.6738948   9.682517 5.765861e-15 1.297319e-14 23.5019164
#IFIH1  -0.9077791 -1.9933273  -4.095474 1.187749e-04 2.137947e-04  0.2405141
#AIM2    0.9589605 -4.0241078   3.664003 4.504313e-04 6.756469e-04 -1.1544674
