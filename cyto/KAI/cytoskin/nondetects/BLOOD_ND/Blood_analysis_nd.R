
setwd("~/Documents/KAI/cytoskin/nondetects/BLOOD_ND")
library(stringr)
library("HTqPCR")#nondetects uses an objects 
library("nondetects")
#controlls
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/datacontrol_nd.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/datapsoriasis_nd.RData")     
###Vitiliigo#

ctrls=DataControllsBlood_nd
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:16]
ctrls <- apply(ctrls[,2:16],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBlood_nd)# psoriasis
rows=rownames(pso[2:16,])
cols=as.character(pso[1,])
ps <- apply(pso[2:16,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
dim(ctrlstst)#16,216
ctrlstst[16,]=gsub("C.*_","",ctrlstst[16,])
ctrlstst[16,]=gsub(".1","",ctrlstst[16,])
ctrlstst[16,]=gsub(".2","",ctrlstst[16,])
ctrlstst[16,]=gsub("h","72h",ctrlstst[16,])

ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[17,]=gsub("C","",ctrlstst[17,])
ctrlstst[17,]=gsub("_.*","",ctrlstst[17,])
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[18,]=gsub("C.*","C",ctrlstst[18,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,17]=as.numeric(as.character(ctrlstst[,17]))
ctrlstst=ctrlstst[order(ctrlstst$V17),]
colnames(ctrlstst)[17]="ID"
colnames(ctrlstst)[16]="Stimulation"
colnames(ctrlstst)[18]="Group"
ctrlstst_or=t(ctrlstst)

#ctrls.data=ctrls.data[1:15,]
#genes=rownames(ctrls.data)
#ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))

#PSO
pstst=rbind(ps, colnames(ps))
pstst[16,]=gsub("P.*_","",pstst[16,])
pstst[16,]=gsub(".1","",pstst[16,])
pstst[16,]=gsub(".2","",pstst[16,])
pstst[16,]=gsub("h","72h",pstst[16,])

pstst=rbind(pstst, colnames(pstst))
pstst[17,]=gsub("P","",pstst[17,])
pstst[17,]=gsub("_.*","",pstst[17,])
pstst=rbind(pstst, colnames(pstst))
pstst[18,]=gsub("P.*","P",pstst[18,])
pstst=as.data.frame(t(pstst))
pstst[,17]=as.numeric(as.character(pstst[,17]))
pstst=pstst[order(pstst$V17),]
colnames(pstst)[17]="ID"
colnames(pstst)[16]="Stimulation"
colnames(pstst)[18]="Group"
pstst_or=t(pstst)

### Combine C and PSO into one dataframe
cp_or=cbind(ctrlstst_or,pstst_or)
cp_or=data.frame(cp_or)
cp_or=cp_or[c(1:12,14:18),]#Remove row belongs to  PIE3
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
#tmp1=exprs(raw)
#tmp1=cp_or[1:15,]
#tmp1=data.frame(tmp1)
##tmp1[is.na(tmp1)]<-40
#tmp1=sapply(tmp1,function(x) (as.character(x)))

#tmp1[is.na(tmp1)]<-c("Undetermined")
#tmp1=data.frame(tmp1)
#dat=tmp1
#View(tmp1)
#dat=sapply(tmp1,function(x) as.numeric(as.character(x)))
#dat=sapply(tmp1,function(x) (as.character(x)))
#dat=data.matrix(dat)
#raw <- new("qPCRset", exprs = tmp1, featureCategory = as.data.frame(array("OK",  dim=dim(tmp1))))
#exprs(raw)=tmp1
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
cpImp=qpcrImpute(raw,groupVars=c("sampleNamesShort","Group","Stimulation"))
exprs(cpImp)
cpImp_SamGrSt=cpImp
save(cpImp_SamGrSt, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/cpImp_SamGrSt_b2mg.RData")
exprs(cpImp_SamGrSt)
cpImp_GrSt=qpcrImpute(raw,groupVars=c("Group","Stimulation"))
save(cpImp_GrSt, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/cpImp_GrSt_b2mg.RData")
exprs(cpImp_GrSt)
#check with other dataset
data(oncogene2013)
#oncogeneImp=qpcrImpute(oncogene2013,groupVars=c("sampleType","treatment"))
exprs(oncogeneImp)

#### Extract imputed expression part and
#assemble a dataframe with mean values of Ct-s and STD of Ct-s per patient/control

#  pomeniat obratno CT=exprs(cpImp_SamGrSt)
CT=exprs(cpImp_GrSt)
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
#0.16621186
#v=c (25.935,25.89700,26.202)
#sd(v)
#[1] 0.1662117
test=head(CT[,1:6])
byapply(test, 3, rowSds)
#              1           2
#[1,] 0.16621186 0.099846718
#[2,] 0.20545625 0.002645379
#[3,] 0.04025324 0.167011515
#[4,] 0.13497780 0.151187228
#[5,] 0.34513130 0.709478320
#[6,] 0.04324789 0.126645328

#bind together CTMean and CTStd
CTMean=t(CTMeans)
CTMean=cbind(rownames(CTMean), CTMean)
colnames(CTMean)[1]="SAMPLE"
CTStd=t(CTStd)
CTStd=cbind(rownames(CTStd), CTStd)
colnames(CTStd)[1]="SAMPLE"
#combine with the reference gene expression.
#create a data frame for the reference gene in C/P samples and calibrator.
#Store calibtaror separately.
#dataframes are located in
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefStd.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefMean.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefCal.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefMean_PSO.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefStd_PSO.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefCal_PSO.RData")

#combine CT values from Controls and Psoriasis
CTRefMean_CP=rbind(CTRefMean,CTRefMean_PSO )
CTBlood=merge(CTMean,CTRefMean_CP, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
rownames(CTBlood)=CTBlood[,1]
rowsCTBlood=rownames(CTBlood)
#combine Std from Controls and Psoriasis
CTRefStd_CP=rbind(CTRefStd,CTRefStd_PSO)
CTStdBlood=merge(CTStd,CTRefStd_CP, by.x="SAMPLE", by.y="SAMPLE", all.x=T ) 
rownames(CTStdBlood)=CTStdBlood[,1]
rowsCTStdBlood=rownames(CTStdBlood)
#combine RefCal for C and P
CTRefCal_CP=rbind(CTRefCal,CTRefCal_PSO)
rownames(CTRefCal_CP)=CTRefCal_CP[,1]
rowsCTRefCal_CP=rownames(CTRefCal_CP)

DDCTBlood=CTBlood[,2:15]#only delta-delta CTs for 14 genes
rownames(DDCTBlood)=CTBlood[,1]
rowsDDCTBlood=rownames(DDCTBlood)

SEMCTBlood=CTBlood[,2:15] #only SEM
rownames(SEMCTBlood)=CTBlood[,1]
rowsSEMCTBlood=rownames(SEMCTBlood)

CTBlood=CTBlood[,2:16] #exclude  first column sample names
CTStdBlood=CTStdBlood[,2:16]#exclude  first column sample names

DDCTBlood=sapply(DDCTBlood,function(x) as.numeric(as.character(x)))
rownames(DDCTBlood)=rowsDDCTBlood
SEMCTBlood=sapply(SEMCTBlood,function(x) as.numeric(as.character(x)))
rownames(SEMCTBlood)=rowsSEMCTBlood

CTBlood=sapply(CTBlood,function(x) as.numeric(as.character(x)))
CTStdBlood=sapply(CTStdBlood,function(x) as.numeric(as.character(x)))
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/CTRefCal.RData")
CTRefCal_CP=CTRefCal_CP[,2:15]
CTRefCal_CP=sapply(CTRefCal_CP,function(x) as.numeric(as.character(x)))

head(CTRefCal_CP)
for(i in 1:length(rowsCTStdBlood)){
  print(i)
  for(j in 1:14){ #14 genes
    print(j)
    print(CTBlood[i,j])
    print(CTBlood[i,15])
    print(CTRefCal[i,j])
    DDCTBlood[i,j]=CTBlood[i,j]-CTBlood[i,15]-CTRefCal_CP[i,j]
    SEMCTBlood[i,j]=sqrt((CTStdBlood[i,j])^2+(CTStdBlood[i,15])^2)/(sqrt(3+3/2))
  }
}

DDCTBloodB2mg=DDCTBlood
save(DDCTBloodB2mg, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/DDCTBloodB2mg.RData")

SEMCTBloodB2mg=SEMCTBlood
save(SEMCTBloodB2mg, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/SEMCTBloodB2mg.RData")

##Filter Delta-delta CTs based on the SEM<0,2


SEMCTBloodB2mg=SEMCTBloodB2mg[complete.cases(SEMCTBloodB2mg),] #4 samples were filtered out because the 
DDCTBloodB2mg=DDCTBloodB2mg[complete.cases(DDCTBloodB2mg),]
DDCTBloodB2mgF=DDCTBloodB2mg

for(i in 1:length(rownames(DDCTBloodB2mg))){
  print(i)
  for(j in 1:14){ #14 genes
    print(j)
    if(SEMCTBloodB2mg[i,j]<= 0.2){
    DDCTBloodB2mgF[i,j]=DDCTBloodB2mg[i,j]
    } else{
      if(SEMCTBloodB2mg[i,j]> 0.2){ 
        DDCTBloodB2mgF[i,j]="NA"
      }
    }
  }
}
#save delta-delta CTs after filtering based on qurtiles of reference gene CTs and SEM<0.2
save(DDCTBloodB2mgF, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/DDCTBloodB2mgF.RData")

###############Differential expression
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/DDCTBloodB2mgF.RData")
library("limma")
library("pheatmap")
data=t(DDCTBloodB2mgF)
genenames=colnames(DDCTBloodB2mgF)
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
pheatmap(data_tmp_or, cluster_rows = F, cluster_cols = F,scale ="row", main="Gene expression(delta-delta CT) in blood samples (B2mg) in controls, psoriasis (72h, PMA, SEB)")


#####analysis of interleukins###########
il=c("IL17A","IL17F", "IL22","IL4", "IL10" )#select interleukins from the data

#analysis of PMA group
il_pma=data_tmp_or[rownames(data_tmp_or)%in%il,c(24:46,106:138)]
dim(il_pma)
pheatmap(il_pma, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression (DDCT) of interleukins in blood samples of PMA/lono stimulated controls and psoriasis ")

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
#         logFC  AveExpr        t     P.Value adj.P.Val         B
#IL17F 1.808852 2.406207 2.729981 0.009232221 0.0461611 -2.667702
write.table(tt,file="diffexp_il_PpmavsCpma.txt", sep="\t", quote=F)

#analysis of SEB group
il_seb=data_tmp_or[rownames(data_tmp_or)%in%il,c(47:70,139:173)]
dim(il_seb)#5 59
pheatmap(il_seb, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression (DDCT) of interleukins in blood samples of SEB stimulated controls and psoriasis ")
 
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
#         logFC  AveExpr        t      P.Value   adj.P.Val          B
#IL17F 1.420363 5.766965 3.498884 0.0006963734 0.003481867 -0.6302377
write.table(tt,file="diffexp_il_PsebvsCseb.txt", sep="\t", quote=F)

#analysis of 72h group
il_seb=data_tmp_or[rownames(data_tmp_or)%in%il,c(1:23,71:105)]
dim(il_seb)#5 59
pheatmap(il_seb, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression (DDCT) of interleukins in blood samples of unstimulated controls and psoriasis ")

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
#logFC   AveExpr       t     P.Value  adj.P.Val         B
#IL4 1.528434 -6.014534 3.02547 0.002918829 0.01459415 -1.731047
write.table(tt,file="diffexp_il_PusvsCus.txt", sep="\t", quote=F)

#analysis of unstimulated group
g_72h=data_tmp_or[(rownames(data_tmp_or))%in%il,c(1:23,71:105)]
dim(g_72h) #9, 58
pheatmap(g_72h, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of unstimulated controls and psoriasis groups ")
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
#logFC   AveExpr       t     P.Value  adj.P.Val         B
#IL4 1.528434 -6.014534 3.02547 0.002918829 0.01459415 -1.731047
write.table(tt,file="diffexp_шд_PusvsCus.txt", sep="\t", quote=F)


########Diff expression analysis of other groups.

il=c("IL17A","IL17F", "IL22","IL4", "IL10")#select interleukins from the data

#analysis of PMA group
g_pma=data_tmp_or[!rownames(data_tmp_or)%in%il,c(24:46,106:138)]
dim(g_pma)#9 56
pheatmap(g_pma, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of PMA/lono stimulated controls and psoriasis groups")

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
#           logFC    AveExpr        t      P.Value   adj.P.Val          B
#PYCARD 0.9941057  0.8387675 3.796798 0.0004345172 0.003910655 -0.1513306
#INFG   1.1321051 -0.3338262 3.026705 0.0039271045 0.017671970 -2.2082382
#TIGIT2 0.5885267  0.9859568 2.786085 0.0071570451 0.021471135 -2.8362368
#IFIH1  0.6514715 -0.8392808 2.405766 0.0201088863 0.045244994 -3.6172157
rownames(tt)
write.table(tt,file="diffexp_genes_PpmavsCpma.txt", sep="\t", quote=F)

#analysis of SEB group
g_seb=data_tmp_or[!rownames(data_tmp_or)%in%il,c(47:70,139:173)]
dim(g_seb)# 9,59
pheatmap(g_seb, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of SEB stimulated controls and psoriasis groups ")

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
rownames(tt)
write.table(tt,file="diffexp_genes_PsebvsCseb.txt", sep="\t", quote=F)


#analysis of unstimulated group
g_72h=data_tmp_or[!(rownames(data_tmp_or))%in%il,c(1:23,71:105)]
dim(g_72h) #9, 58
pheatmap(g_72h, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of unstimulated controls and psoriasis groups ")
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
write.table(tt,file="diffexp_genes_PusvsCus.txt", sep="\t", quote=F)


####Diffexp PSO vs PSO
#Pso PMA vs Pso SEB
data.pso=data_tmp_or
il_pso=data.pso[!rownames(data.pso)%in%il,c(106:173)]
dim(il_pso)
il_pso_t=t(il_pso)#9 68
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of PMA/lono and SEB stimulated psoriasis groups ")

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
#           logFC    AveExpr          t      P.Value    adj.P.Val          B
#PYCARD -5.003200  4.2044096 -31.421420 2.034668e-40 1.831202e-39 81.6285299
#OAS2   -2.243846  0.5942944 -16.290276 9.008430e-26 4.053794e-25 48.2576147
#FOXP3  -3.056865  2.8184915 -11.641468 4.420063e-18 1.326019e-17 30.5619290
#IFIH1  -1.773391  0.3117506  -8.997323 9.254495e-13 2.082261e-12 18.4694372
#CTLA4  -1.900979  2.4488088  -6.378476 1.477824e-08 2.660083e-08  8.6912405
#TIGIT2 -0.949295  1.6894758  -5.333176 1.107562e-06 1.661343e-06  4.4655602
#EOMES  -1.324098  1.9302238  -4.749922 1.018516e-05 1.309520e-05  2.2913183
#INFG    1.146008 -0.4749510   4.142643 1.100586e-04 1.238159e-04  0.1164017
rownames(tt)#"PYCARD" "OAS2"   "FOXP3"  "IFIH1"  "CTLA4"  "TIGIT2" "EOMES"  "INFG"  
write.table(tt,file="diffexp_genes_PpmavsPseb.txt", sep="\t", quote=F)

#Pso PMA vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(71:138)]
dim(il_pso)# 9 68
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of PMA/lono stimulated and unstimulated psoriasis groups ")

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
#            logFC    AveExpr          t      P.Value    adj.P.Val         B
#PYCARD -6.1408897  4.8793104 -36.581430 3.009439e-45 2.708495e-44 92.124118
#INFG   13.2231708 -5.9313123  44.780269 1.373087e-44 6.178892e-44 89.467489
#CTLA4   4.9050404 -0.9821566  16.360102 5.120426e-26 1.536128e-25 48.844085
#AIM2    3.2422476 -1.2294706  13.191097 3.692520e-21 8.308169e-21 37.657541
#FOXP3   1.5528420  0.4437940   5.814768 1.524817e-07 2.744670e-07  6.490800
#TIGIT2  0.6874188  0.8610098   3.609944 5.528643e-04 8.292965e-04 -1.471818
#IFIH1  -0.7179695 -0.2411006  -3.101904 2.904817e-03 3.734764e-03 -2.907541
#OAS2   -0.3840584 -0.3523445  -2.322687 2.301200e-02 2.588850e-02 -4.893890
rownames(tt)#"PYCARD" "INFG"   "CTLA4"  "AIM2"   "FOXP3"  "TIGIT2" "IFIH1"  "OAS2"
write.table(tt,file="diffexp_genes_PpmavsPus.txt", sep="\t", quote=F)



#Pso SEB vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(71:105, 139:173)]
dim(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of SEB stimulated and unstimulated psoriasis groups ")

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
#           logFC     AveExpr         t      P.Value    adj.P.Val         B
#INFG   12.077163 -6.22645673 41.900846 3.180196e-45 2.862176e-44 91.589408
#CTLA4   6.806020  0.06842157 24.077861 4.128112e-37 1.857651e-36 74.208776
#FOXP3   4.609707  1.99575439 19.304370 6.957990e-31 2.087397e-30 59.971682
#AIM2    2.839973 -1.40350725 12.534398 3.286349e-20 7.394285e-20 35.425140
#OAS2    1.859788  0.78345748 11.217274 1.190968e-17 2.143743e-17 29.563674
#TIGIT2  1.636714  1.32169707  8.453460 1.588350e-12 2.382525e-12 17.781103
#PYCARD -1.137690  6.80845448 -8.325809 2.333218e-12 2.999851e-12 17.374729
#EOMES   1.338324  1.89346667  5.560071 3.704385e-07 4.167434e-07  5.531309
#IFIH1   1.055421  0.64832573  4.802251 9.753400e-06 9.753400e-06  2.473271
rownames(tt)#"INFG"   "CTLA4"  "FOXP3"  "AIM2"   "OAS2"   "TIGIT2" "PYCARD" "EOMES"  "IFIH1" 
write.table(tt,file="diffexp_genes_PsebvsPus.txt", sep="\t", quote=F)


##### Control groups comparison ##########
#Pso SEB vs Pso unstimulated
il_pso=data.pso[!rownames(data.pso)%in%il,c(1:23, 47:70)]
dim(il_pso)
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of SEB stimulated and unstimulated control groups ")

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
write.table(tt,file="diffexp_genes_CsebvsCus.txt", sep="\t", quote=F)
rownames(tt)
#"CTLA4"  "FOXP3"  "INFG"   "PYCARD" "AIM2"   "OAS2"   "TIGIT2" "IFIH1"  "EOMES" 

##CpmavsCus

il_pso=data.pso[!rownames(data.pso)%in%il,c(1:46)]
dim(il_pso)# 9 68
pheatmap(il_pso, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression(DDCT) in blood samples of PMA/lono stimulated and unstimulated control groups ")

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
write.table(tt,file="diffexp_genes_CpmavsCus.txt", sep="\t", quote=F)
rownames(tt)
#"PYCARD" "CTLA4"  "INFG"   "AIM2"   "OAS2"   "IFIH1"  "FOXP3" 