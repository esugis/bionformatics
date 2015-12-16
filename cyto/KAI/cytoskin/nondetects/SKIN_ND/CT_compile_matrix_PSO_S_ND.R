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
rowsP=as.character(data.cols[,2])
DataPsoriasisSkin=rowsP
DataPsoriasisSkin=as.data.frame(DataPsoriasisSkin)
colnames(DataPsoriasisSkin)=c("SAMPLE")
DataPsoriasisSkin$SAMPLE=make.names(as.vector(as.character(DataPsoriasisSkin$SAMPLE)), unique=TRUE)
library(functional)
#for each in psoriasis
for(i in 1:length(psoriasis)){  
  #for(i in 1:3){
  f=psoriasis[i] 
  #f=psoriasis[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_2.csv","",gene)
  CTdata=read.csv(file=f)
  CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
  CTdata$SAMPLE=as.character(CTdata$SAMPLE)
  CTdata[,2]=as.numeric(as.character(CTdata$CT))
  df=CTdata[,1:2]
  df$SAMPLE=make.names(as.vector(as.character(df$SAMPLE)), unique=TRUE)
  colnames(df)[2]=gene 
  DataPsoriasisSkin=merge(DataPsoriasisSkin,df, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  ##print(DataControllsSkin)
}
DataPsoriasisSkin_ND=DataPsoriasisSkin
write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataPsoriasisSkin_ND.txt",DataPsoriasisSkin_ND, sep="\t", quote=F )      
save(DataPsoriasisSkin_ND, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataPsoriasisSkin_ND.RData")
