#This script reads csv files that contain info about expression of controlls, psoriasis and vitiligo
#for each of the genes.
#Files are in the format gene_names_1.csv, gene_names_2.csv, gene_names_3.csv
#where: 
#*_1.csv corresponds to controls1
#*_2.csv controls2
#*_3.csv psoriasis1
#*_4.csv psoriasis2
#*_5.csv vitiliigo1
#*_6.csv vitiliigo2

#Csv are obtained from xlsx provided by the collaborators.
#Script
#-Reads the files

#list all files with csv from the folder
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_blood_bact", full.names = TRUE)
#take corresponding _1,_2,_3, _4, _5, _6 for further analysis
controlls1=csvFiles[grepl("*_1.csv", csvFiles)]
controlls2=csvFiles[grepl("*_2.csv", csvFiles)]
psoriasis1=csvFiles[grepl("*_3.csv", csvFiles)]
psoriasis2=csvFiles[grepl("*_4.csv", csvFiles)]
vitiliigo1=csvFiles[grepl("*_5.csv", csvFiles)]  
vitiliigo2=csvFiles[grepl("*_6.csv", csvFiles)]  
# Assemble the matrix for each of the groups. Combine controlls1,2, psoriasis 1,2, vitiligo 1,2
#read in row names
CTdata1=read.csv(file=psoriasis1[1], na.strings = "")
Calref1=as.numeric(as.character(CTdata1[5,11]))
CTdata1=CTdata1[7:159,2:16]

CTdata2=read.csv(file=psoriasis2[1], na.strings = "")
Calref1=as.numeric(as.character(CTdata1[5,11]))
CTdata2=CTdata2[7:168,2:16]

rn1=as.character(as.vector(CTdata1[,1]))#sample names from the first csv
rn2=as.character(as.vector(CTdata2[,1]))#sample names from the second csv
rowsC=unique(c(rn1, rn2))
length(rowsC)
#[1] 105
CTRefMean=rowsC
CTRefMean=as.data.frame(CTRefMean)
colnames(CTRefMean)=c("SAMPLE")
CTRefMean$SAMPLE=as.character(CTRefMean$SAMPLE)
CTRefStd=CTRefMean
CTRefCal=CTRefMean
library(functional)
#for each of controlls
#selected.patient=c()
#selected.expression=c()
for(i in 1:length(psoriasis1)){  
  #for(i in 3:5){
  f1=psoriasis1[i] 
  f2=psoriasis2[i]
  #f1=psoriasis1[8]
  #f2=psoriasis2[8]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_blood_bact/CtkalkulaatorPBMC_","",f1)
  gene=gsub("_Bact_3.csv","",gene)
  print(gene)
  
  CTdata1=read.csv(file=f1, na.strings = "")
  Calref1=as.numeric(as.character(CTdata1[5,11]))
  CTdata1=CTdata1[7:159,2:16]
  head(CTdata1)
  CTdata1$CalRef=replicate(length(row.names(CTdata1)), Calref1)
 
  CTdata2=read.csv(file=f2, na.strings = "")
  Calref2=as.numeric(as.character(CTdata2[5,11])) 
  CTdata2=CTdata2[7:168,2:16]
  CTdata2$CalRef=replicate(length(row.names(CTdata2)), Calref2)
  CTdata=rbind(CTdata1, CTdata2)
 
  df=CTdata
  dft=t(df)
  dfr=as.data.frame(dft)
  dfr=t(dfr[c(F,F,T)])#### selects only the rows containing mean values and values related to STD etc from the data.
  rows=as.character(dfr[,1])
  dim(dfr)
    
## filter based on quartiles
  d=as.data.frame(dfr)
  d[,8]=as.numeric(as.character(d[, 8]))
  dfr.f=d[d[,8] > quantile(d[,8], .25, na.rm=T) - 1.5*IQR(d[,8],na.rm=T) & d[,8] < quantile(d[,8], .75, na.rm=T) + 1.5*IQR(d[,8],na.rm=T),]
  print("dfr.f dim")
  print(dim(dfr.f))
 
CTRef=as.data.frame(dfr.f[,c(1,8,7,16)])
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


CTRefMean_PSO_bact=CTRefMean[,c(1,2)]
dim(CTRefMean_PSO_bact)
colnames(CTRefMean_PSO_bact)[2]="CTRef"
colnames(CTRefMean_PSO_bact)[1]="SAMPLE"
#CTRefMean$Order=gsub("_.*","",CTRefMean$)
save(CTRefMean_PSO_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefMean_PSO_bact.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefMean_PSO_bact.txt",CTRefMean_PSO_bact, sep="\t", quote=F )      
CTRefStd_PSO_bact=CTRefStd[,c(1,2)]
dim(CTRefStd_PSO_bact)
colnames(CTRefStd_PSO_bact)[2]="STDRef"
colnames(CTRefStd_PSO_bact)[1]="SAMPLE"
save(CTRefStd_PSO_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefStd_PSO_bact.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_N_bact/CTReStd_PSO_bact.txt",CTRefStd_PSO_bact, sep="\t", quote=F )      
CTRefCal_PSO_bact=CTRefCal
dim(CTRefCal_PSO_bact)
save(CTRefCal_PSO_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefCal_PSO_bact.RData")
#write.table(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/CTRefCal_PSO_bact.txt",CTRefCal_PSO_bact, sep="\t", quote=F )   

