#This script reads csv files that contain info about expression of controlls
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
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_blood2", full.names = TRUE)
#take corresponding _1,_2,_3, _4, _5, _6 for further analysis
controlls1=csvFiles[grepl("*_1.csv", csvFiles)]
controlls2=csvFiles[grepl("*_2.csv", csvFiles)]
psoriasis1=csvFiles[grepl("*_3.csv", csvFiles)]
psoriasis2=csvFiles[grepl("*_4.csv", csvFiles)]
vitiliigo1=csvFiles[grepl("*_5.csv", csvFiles)]  
vitiliigo2=csvFiles[grepl("*_6.csv", csvFiles)]  
# Assemble the matrix for each of the groups. 
#read in row names
CTdata1=read.csv(file=controlls1[1], na.strings = "")
CTdata1=CTdata1[7:114,2:ncol(CTdata1)]
CTdata2=read.csv(file=controlls2[1], na.strings = "")
CTdata2=CTdata2[7:123,2:ncol(CTdata2)]
rn1=as.character(as.vector(CTdata1[,1]))#sample names from the first csv
rn2=as.character(as.vector(CTdata2[,1]))#sample names from the second csv
rowsC=c(rn1, rn2)
length(rowsC)
#[1] 225
DataControllsBlood=rowsC
DataControllsBlood=as.data.frame(DataControllsBlood)
colnames(DataControllsBlood)=c("SAMPLE")
DataControllsBlood$SAMPLE=make.names(as.vector(as.character(DataControllsBlood$SAMPLE)), unique=TRUE)
library(functional)
for(i in 1:length(controlls1)){  
  #for(i in 1:2){
  #i=1
  f1=controlls1[i] 
  f2=controlls2[i]
  #f1=controlls1[8]
  #f2=controlls2[8]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_blood2/CtkalkulaatorPBMC_","",f1)
  gene=gsub("_B2mg_1.csv","",gene)
  print(gene)
  CTdata1=read.csv(file=f1, na.strings = "")
  CTdata1=CTdata1[7:114,2:3]
 # dim(CTdata1)
  CTdata1$SAMPLE=as.character(CTdata1$SAMPLE)
  CTdata1[,2]=as.numeric(as.character(CTdata1$CT))
  CTdata2=read.csv(file=f2, na.strings = "")
  CTdata2=CTdata2[7:123,2:3]
 # dim(CTdata2)
  CTdata2$SAMPLE=as.character(CTdata2$SAMPLE)
  CTdata2[,2]=as.numeric(as.character(CTdata2$CT))
  CTdata=rbind(CTdata1, CTdata2)
  df=CTdata
 # dim(CTdata)
  colnames(df)[2]=gene 
 # print(head(df))
 # str(df)
  df$SAMPLE=make.names(as.vector(as.character(df$SAMPLE)), unique=TRUE)
 # print(head(DataControllsBlood))
 # str(DataControllsBlood)
  DataControllsBlood=merge(DataControllsBlood,df, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
 # print(DataControllsBlood)
}
#rownames(DataControllsBlood)=make.names(as.vector(as.character(DataControllsBlood$SAMPLE)), unique=TRUE)
DataControllsBlood_nd=DataControllsBlood[c(1:144,154:225),]#remove C024 from the data.
dim(DataControllsBlood_nd)
save(DataControllsBlood_nd, file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/datacontrol_nd.RData")
write.table(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/datacontrol_nd.txt",DataControllsBlood_nd, sep="\t", quote=F )      



