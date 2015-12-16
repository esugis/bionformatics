

library(stringr)
library("HTqPCR")#nondetects uses an objects 
library("nondetects")
#controlls
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataControllsSkin_nd.Rdata")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataPsoriasisSkin_ND.RData")     


ctrls=DataControllsSkin_nd
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:43]
ctrls <- apply(ctrls[,2:43],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisSkin_ND)# psoriasis
rows=rownames(pso[2:43,])
cols=as.character(pso[1,])
ps <- apply(pso[2:43,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows
colnames(ps)=str_replace(colnames(ps),"PH","PL")#sick
colnames(ps)=str_replace(colnames(ps),"PT","PNL")#healthy
rownames(ps)=gsub("IL8b","CXCL8", rownames(ps))


#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
dim(ctrlstst)#43, 72
ctrlstst[43,]=gsub("\\.2","",ctrlstst[43,])
ctrlstst[43,]=gsub("\\.1","",ctrlstst[43,])
ctrlstst[43,]=gsub("C00","",ctrlstst[43,])
ctrlstst[43,]=gsub("C0","",ctrlstst[43,])

ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[44,]=gsub("C.*","C",ctrlstst[44,])
#ctrlstst=as.data.frame(t(ctrlstst))

ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[45,]=gsub("C.*","C",ctrlstst[45,])

ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,43]=as.numeric(as.character(ctrlstst[,43]))
ctrlstst=ctrlstst[order(ctrlstst$V43),]
colnames(ctrlstst)[43]="ID"
colnames(ctrlstst)[44]="Type"
colnames(ctrlstst)[45]="Group"
ctrlstst_or=t(ctrlstst)


#PSO
pstst=rbind(ps, colnames(ps))
pstst[43,]=gsub("\\.2","",pstst[43,])
pstst[43,]=gsub("\\.1","",pstst[43,])
pstst[43,]=gsub("PL00","",pstst[43,])
pstst[43,]=gsub("PNL00","",pstst[43,])
pstst[43,]=gsub("PL0","",pstst[43,])
pstst[43,]=gsub("PNL0","",pstst[43,])

pstst=rbind(pstst, colnames(pstst))
pstst[44,]=gsub("PL.*","PL",pstst[44,])
pstst[44,]=gsub("PNL.*","PNL",pstst[44,])

pstst=rbind(pstst, colnames(pstst))
pstst[45,]=gsub("P.*","P",pstst[45,])

pstst=as.data.frame(t(pstst))
pstst[,43]=as.numeric(as.character(pstst[,43]))
pstst=pstst[order(pstst$V43),]
colnames(pstst)[43]="ID"
colnames(pstst)[44]="Type"
colnames(pstst)[45]="Group"
pstst_or=t(pstst)

### Combine C and PSO into one dataframe
cp_or=cbind(ctrlstst_or,pstst_or)
cp_or=data.frame(cp_or)
# Convert data to qPCRset

head(cp_or)

Sample=(colnames(cp_or)) 
Sample=gsub("\\.1","",Sample)
Sample=gsub("\\.2","",Sample)

cp_or=t(cp_or)
cp_or=cbind(cp_or,Sample)

cp_or=t(cp_or)
cp_or=data.frame(cp_or)
cp_or=t(cp_or)
tmp1=cp_or[1:42,]
tmp1=data.frame(tmp1)

tmp1=sapply(tmp1,function(x) as.numeric(as.character(x)))
tmp1[is.na(tmp1)]<-40
#mat=data.matrix(dat)
#dat=sapply(tmp1,function(x) as.numeric(as.character(x)))
dat=tmp1
mat=data.matrix(dat)
str(mat)
rownames(mat)=rownames(cp_or[1:42,])
#mat=t(mat)
dimnames(mat)=NULL
#raw <- new("qPCRset", exprs = mat, featureCategory = as.data.frame(array("OK",  dim=dim(mat))))
raw <- new("qPCRset", exprs = mat)
raw@phenoData@data$sampleNamesShort <-Sample
raw@phenoData@data$sampleName <-as.character(colnames(cp_or))
tmp=t(cp_or)
Group=as.character(tmp[,45])
raw@phenoData@data$Group <- Group
Type=as.character(tmp[,44])
raw@phenoData@data$Type <- Type
featureNames(raw) <- rownames(cp_or)[1:42]
colnames(exprs(raw))=colnames(cp_or)
#set feature category
FC=data.frame(exprs(raw))
FC[!FC==40]<-c("OK")
FC[FC==40]<-c("Undetermined")
featureCategory(raw)=FC#check what features are undetermined or needed to be indicated

#set feature type
FT=replicate(42,"Target")
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

cpImp_GrSt=qpcrImpute(raw,groupVars=c("Group","Type"))
exprs(cpImp_GrSt)
cpImp_GrSt_skin=cpImp_GrSt
save(cpImp_GrSt_skin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/cpImp_GrSt_skin.RData")
#check with other dataset
data(oncogene2013)
#oncogeneImp=qpcrImpute(oncogene2013,groupVars=c("sampleType","treatment"))
exprs(oncogeneImp)
data(sagmb2011)
exprs(sagmb2011)
#### Extract imputed expression part and
#assemble a dataframe with mean values of Ct-s and STD of Ct-s per patient/control

CT=exprs(cpImp_GrSt_skin)
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
dim(CTMeans) #42 93
meanNames=colnames(CT)
meanNames=gsub("\\.2","",meanNames)
meanNames=gsub("\\.1","",meanNames)
meanNames=unique(meanNames)
length(meanNames) #93
colnames(CTMeans)=meanNames

#calculate SD of 3 replicates per patient
CTStd <- byapply(CT, 3, sd)
mystd<-function(x){apply(x,1,sd)}
CTStd <- byapply(CT, 3, mystd)
dim(CTStd) #42 93
colnames(CTStd)=meanNames
#test if the values are correct
library("matrixStats")

#bind together CTMean and CTStd
CTMean_skin=t(CTMeans)
CTMean_skin=cbind(rownames(CTMean_skin), CTMean_skin)
colnames(CTMean_skin)[1]="SAMPLE"
CTStd_skin=t(CTStd)
CTStd_skin=cbind(rownames(CTStd_skin), CTStd_skin)
colnames(CTStd_skin)[1]="SAMPLE"
#combine with the reference gene expression.
#create a data frame for the reference gene in C/P samples and calibrator.
#Store calibtaror separately.
#dataframes are located in
!!!!!
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefStdNDRef.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMeanSkinNDRef.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCalNDRef.RData")

load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMean_skin.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefStd_skin.RData")
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCal_skin.RData")

#combine CT values from Controls and Psoriasis
colnames(CTRefMeanSkinNDRef)[2]="CTRef"
CTRefMean_CP_skin=rbind(CTRefMean_skin,CTRefMeanSkinNDRef)
CTSkin=merge(CTMean_skin,CTRefMean_CP_skin, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
rownames(CTSkin)=CTSkin[,1]
rowsCTSkin=rownames(CTSkin)
#combine Std from Controls and Psoriasis
colnames(CTRefStdNDRef)[2]="STDRef"
CTRefStd_CP_skin=rbind(CTRefStd_skin,CTRefStdNDRef)
CTStdSkin=merge(CTStd_skin,CTRefStd_CP_skin, by.x="SAMPLE", by.y="SAMPLE", all.x=T ) 
rownames(CTStdSkin)=CTStdSkin[,1]
rowsCTStdSkin=rownames(CTStdSkin)




#!!!!!!
#combine RefCal for C and P
CTRefCal_CP_skin=rbind(CTRefCal_skin,CTRefCalNDRef)
rownames(CTRefCal_CP_skin)=CTRefCal_CP_skin[,1]
rowsCTRefCal_CP_skin=rownames(CTRefCal_CP_skin)

DDCTSkin=CTSkin[,2:43]#only delta-delta CTs for 42 genes
rownames(DDCTSkin)=CTSkin[,1]
rowsDDCTSkin=rownames(DDCTSkin)

SEMCTSkin=CTSkin[,2:43] #only SEM
rownames(SEMCTSkin)=CTSkin[,1]
rowsSEMCTSkin=rownames(SEMCTSkin)

CTSkin=CTSkin[,2:44] #exclude  first column sample names
CTStdSkin=CTStdSkin[,2:44]#exclude  first column sample names

DDCTSkin=sapply(DDCTSkin,function(x) as.numeric(as.character(x)))
rownames(DDCTSkin)=rowsDDCTSkin
SEMCTSkin=sapply(SEMCTSkin,function(x) as.numeric(as.character(x)))
rownames(SEMCTSkin)=rowsSEMCTSkin

CTSkin=sapply(CTSkin,function(x) as.numeric(as.character(x)))
CTStdSkin=sapply(CTStdSkin,function(x) as.numeric(as.character(x)))
#!!!
#load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCal_skin.RData")
CTRefCal_CP_skin=CTRefCal_CP_skin[,2:43]
CTRefCal_CP_skin=sapply(CTRefCal_CP_skin,function(x) as.numeric(as.character(x)))

head(CTRefCal_CP_skin)
for(i in 1:length(rowsCTStdSkin)){
  print(i)
  for(j in 1:42){ #42  genes
    print(j)
    print(CTSkin[i,j])
    print(CTSkin[i,15])
   # print(CTRefCal_skin[i,j])
    DDCTSkin[i,j]=CTSkin[i,j]-CTSkin[i,15]-CTRefCal_CP_skin[i,j]
    SEMCTSkin[i,j]=sqrt((CTStdSkin[i,j])^2+(CTStdSkin[i,15])^2)/(sqrt(3+3/2))
  }
}

#save results
save(DDCTSkin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DDCTSkin.RData")
save(SEMCTSkin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/SEMCTSkin.RData")

##Filter Delta-delta CTs based on the SEM<0,2

SEMCTSkin=SEMCTSkin[complete.cases(SEMCTSkin),] # samples were filtered out because the 
DDCTSkin=DDCTSkin[complete.cases(DDCTSkin),]
DDCTSkin_F=DDCTSkin

for(i in 1:length(rownames(DDCTSkin))){
  print(i)
  for(j in 1:42){ #42 genes
    print(j)
    if(SEMCTSkin[i,j]<= 0.2){
    DDCTSkin_F[i,j]=DDCTSkin[i,j]
    } else{
      if(SEMCTSkin[i,j]> 0.2){ 
        DDCTSkin_F[i,j]="NA"
      }
    }
  }
}
#save delta-delta CTs after filtering based on qurtiles of reference gene CTs and SEM<0.2
save(DDCTSkin_F, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DDCTSkin_F.RData")

###############Differential expression
library("limma")
library("pheatmap")
data=t(DDCTSkin_F)
genenames=colnames(DDCTSkin_F)
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
colnames(data_tmp)[c(43,44,45)]=c("Group","ID","Type")
library(stringr)
data_tmp=data.frame(data_tmp)

#!!!!!!!!
#  !!!!!!!!!
#extract sample group
data_tmp$Group=str_replace(data_tmp$Group,"C.*","C")
data_tmp$Group=str_replace(data_tmp$Group,"P.*","P")
#extract sample number ID
data_tmp$ID=str_replace(data_tmp$ID,"\\.2","")
data_tmp$ID=str_replace(data_tmp$ID,"\\.1","")
data_tmp$ID=str_replace(data_tmp$ID,"PL00","")#some of the psoriasis numeric IDs have 0 in the beginning
data_tmp$ID=str_replace(data_tmp$ID,"PNL00","")
data_tmp$ID=str_replace(data_tmp$ID,"PL0","")
data_tmp$ID=str_replace(data_tmp$ID,"PNL0","")
data_tmp$ID=str_replace(data_tmp$ID,"C00","")
data_tmp$ID=str_replace(data_tmp$ID,"C0","")
data_tmp$ID=as.numeric(as.character(data_tmp$ID))
#extract sample type
data_tmp$Type=str_replace(data_tmp$Type,"PL.*","PL")
data_tmp$Type=str_replace(data_tmp$Type,"PNL.*","PNL")
data_tmp$Type=str_replace(data_tmp$Type,"C.*","C")
#order based on Group, Treatment, ID
data_tmp=data_tmp[order(data_tmp$Group,data_tmp$Type, data_tmp$ID),]
data_tmp_or=sapply(data_tmp,function(x) as.numeric(as.character(x)))
data_tmp_or=data_tmp_or[,1:42]
rownames(data_tmp_or)=rownames(data_tmp)
data_tmp_or=t(data_tmp_or)
pheatmap(data_tmp_or, cluster_rows = F, cluster_cols = F,scale ="row", main="Gene expression(delta-delta CT) in skin samples in C, PL, PNL")
data_skin_cp_nd=data_tmp_or
save(data_skin_cp_nd, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
!!!!!!take gene expression from skin
#####analysis of interleukins###########
#Differential expression analysis of skin samples after imputation is described in file
# ~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DiffExp_S_ND_S