#This script reads csv files that contain info about expression of controlls, psoriasis and vitiligo
#for each of the genes.
#Files are in the format gene_names_1.csv, gene_names_2.csv, gene_names_3.csv
#where: 


#Csv are obtained from xlsx provided by the collaborators.
#Script
#-Reads the files
library(functional)
#Csv are obtained from xlsx provided by the collaborators.
#Script
#-Reads the files
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_skin", full.names = TRUE)
#take corresponding _1,_2,_3 for further analysis
controlls=csvFiles[grepl("*_1.csv", csvFiles)]
psoriasis=csvFiles[grepl("*_2.csv", csvFiles)]
vitiliigo=csvFiles[grepl("*_3.csv", csvFiles)]  

###Psofiasis#

data.cols=read.csv(file=psoriasis[2])
data.cols=data.cols[7:nrow(data.cols),]
rowsP=unique(as.character(data.cols[,2]))
DataPsoriasisSkin=rowsP
DataPsoriasisSkin=as.data.frame(DataPsoriasisSkin)
colnames(DataPsoriasisSkin)=c("SAMPLE")
DataPsoriasisSkin$SAMPLE=make.names(as.vector(as.character(DataPsoriasisSkin$SAMPLE)), unique=TRUE)
CTRefMean=rowsP
CTRefMean=as.data.frame(CTRefMean)
colnames(CTRefMean)=c("SAMPLE")
CTRefMean$SAMPLE=as.character(CTRefMean$SAMPLE)
CTRefStd=CTRefMean
CTRefCal=CTRefMean

library(functional)
#for each in psoriasis
for(i in 1:length(psoriasis)){  
  #for(i in 1:3){
  f=psoriasis[i] 
  #f=psoriasis[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_2.csv","",gene)
  CTdata=read.csv(file=f)
  Calref=as.numeric(as.character(CTdata[5,11]))
  CTdata$CalRef=replicate(length(row.names(CTdata)), Calref)
  
  CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
  CTdata$SAMPLE=as.character(CTdata$SAMPLE)
  
  df=CTdata
  dft=t(df)
  dfr=as.data.frame(dft)
  dfr=t(dfr[c(F,F,T)])
  rows=as.character(dfr[,1])
  dim(dfr)
  
  ## filter on quartiles
  d=as.data.frame(dfr)
  d[,8]=as.numeric(as.character(d[, 8]))
  dfr.f=d[d[,8] > quantile(d[,8], .25) - 1.5*IQR(d[,8]) & d[,8] < quantile(d[,8], .75) + 1.5*IQR(d[,8]),]
  
  CTRef=as.data.frame(dfr.f[,c(1,8,7,17)])
  CTRef$SAMPLE=as.character(CTRef$SAMPLE)
  CTRef$X.CT..1=as.numeric(as.character(CTRef$X.CT..1))
  CTRef$STD.1=as.numeric(as.character(CTRef$STD.1))
  CTRef$CalRef=as.numeric(as.character(CTRef$CalRef))
  
  colnames(CTRef)[2]=gene 
  colnames(CTRef)[3]=gene 
  colnames(CTRef)[4]=gene 
  
  CTRefStdTmp=CTRef[,c(1,3)]
  str(CTRefStdTmp)
  CTRefMeanTmp=CTRef[,c(1,2)]
  CTRefCalTmp=CTRef[,c(1,4)]
  # print(CTRef)
  # print(CTRefMean)
  # print(CTRefMeanTmp)
  CTRefMean=merge(CTRefMean,CTRefMeanTmp, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  CTRefStd=merge(CTRefStd,CTRefStdTmp, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  CTRefCal=merge(CTRefCal,CTRefCalTmp, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
}

CTRefMeanSkinNDRef=CTRefMean[,1:2]
write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMeanSkinNDRef.txt",CTRefMeanSkinNDRef, sep="\t", quote=F )      
save(CTRefMeanSkinNDRef, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMeanSkinNDRef.RData")

CTRefStdNDRef=CTRefStd[,1:2]
write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefStdNDRef.txt",CTRefStdNDRef, sep="\t", quote=F )      
save(CTRefStdNDRef, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefStdNDRef.RData")

CTRefCalNDRef=CTRefCal
write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCalNDRef.txt",CTRefCalNDRef, sep="\t", quote=F )      
save(CTRefCalNDRef, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCalNDRef.RData")

