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
CTdata1=CTdata1[7:159,2:16]
CTdata2=read.csv(file=psoriasis2[1], na.strings = "")
CTdata2=CTdata2[7:168,2:16]
rn1=as.character(as.vector(CTdata1[,1]))#sample names from the first csv
rn2=as.character(as.vector(CTdata2[,1]))#sample names from the second csv
rowsC=c(rn1, rn2)
length(rowsC)
#[1] 315
DataPsoriasisBlood=rowsC
DataPsoriasisBlood=as.data.frame(DataPsoriasisBlood)
colnames(DataPsoriasisBlood)=c("SAMPLE")
DataPsoriasisBlood$SAMPLE=make.names(as.vector(as.character(DataPsoriasisBlood$SAMPLE)), unique=TRUE)
library(functional)
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
  CTdata1=CTdata1[7:159,2:3]
  CTdata1$SAMPLE=as.character(CTdata1$SAMPLE)
  CTdata1[,2]=as.numeric(as.character(CTdata1$CT))
  CTdata2=read.csv(file=f2, na.strings = "")
  CTdata2=CTdata2[7:168,2:3]
  CTdata2$SAMPLE=as.character(CTdata2$SAMPLE)
  CTdata2[,2]=as.numeric(as.character(CTdata2$CT))
  CTdata=rbind(CTdata1, CTdata2)
  df=CTdata
  colnames(df)[2]=gene 
  df$SAMPLE=make.names(as.vector(as.character(df$SAMPLE)), unique=TRUE)
  DataPsoriasisBlood=merge(DataPsoriasisBlood,df, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  #print(DataPsoriasisBlood)
}
DataPsoriasisBlood_nd_bact=DataPsoriasisBlood
save(DataPsoriasisBlood_nd_bact, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/datapsoriasis_nd_bact.RData")
write.table(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/datapsoriasis_nd_bact.txt",DataPsoriasisBlood_nd_bact, sep="\t", quote=F )      

