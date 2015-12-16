#This script reads csv files that contain info about expression of controlls, psoriasis and vitiligo
#for each of the genes.
#Files are in the format gene_names_1.csv, gene_names_2.csv, gene_names_3.csv
#where: 

library(functional)
#Csv are obtained from xlsx provided by the collaborators.
#Script
#-Reads the files
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_skin", full.names = TRUE)
#take corresponding _1,_2,_3 for further analysis
controlls=csvFiles[grepl("*_1.csv", csvFiles)]
psoriasis=csvFiles[grepl("*_2.csv", csvFiles)]
vitiliigo=csvFiles[grepl("*_3.csv", csvFiles)]  

CTdata=read.csv(file=controlls[1], na.strings = "")
CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
rowsC=unique(as.character(as.vector(CTdata[,1])))
DataControllsSkin=rowsC
DataControllsSkin=as.data.frame(DataControllsSkin)
colnames(DataControllsSkin)=c("SAMPLE")
DataControllsSkin$SAMPLE=make.names(as.vector(as.character(DataControllsSkin$SAMPLE)), unique=TRUE)
CTRefMean=rowsC
CTRefMean=as.data.frame(CTRefMean)
colnames(CTRefMean)=c("SAMPLE")
CTRefMean$SAMPLE=as.character(CTRefMean$SAMPLE)
CTRefStd=CTRefMean
CTRefCal=CTRefMean
for(i in 1:length(controlls)){  
  #for(i in 1:3){
  f=controlls[i] 
  #print(f)
  #f=controlls[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_1.csv","",gene)
  print(gene)
  CTdata=read.csv(file=f, na.strings = "")
  
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
  dfr.f=d[d[,8] > quantile(d[,8], .25, na.rm = T) - 1.5*IQR(d[,8],na.rm = T) & d[,8] < quantile(d[,8], .75, na.rm = T) + 1.5*IQR(d[,8],na.rm = T),]
  
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
CTRefMean_skin=CTRefMean[c(1:23,25),c(1,2)]
dim(CTRefMean_skin) # 
colnames(CTRefMean_skin)[2]="CTRef"
colnames(CTRefMean_skin)[1]="SAMPLE"

save(CTRefMean_skin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMean_skin.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefMean_skin.txt",CTRefMean_skin, sep="\t", quote=F )      

CTRefStd_skin=CTRefStd[c(1:23,25),c(1,2)]#
dim(CTRefStd_skin)#
colnames(CTRefStd_skin)[2]="STDRef"
colnames(CTRefStd_skin)[1]="SAMPLE"
save(CTRefStd_skin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefStd_skin.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTReStd_skin.txt",CTRefStd_skin, sep="\t", quote=F )      
CTRefCal_skin=CTRefCal[c(1:23,25),]#
dim(CTRefCal_skin) # 24, 43
save(CTRefCal_skin, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCal_skin.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/CTRefCal_skin.txt",CTRefCal_skin, sep="\t", quote=F )      
