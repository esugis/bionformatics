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
rowsC=as.character(as.vector(CTdata[,1]))
DataControllsSkin=rowsC
DataControllsSkin=as.data.frame(DataControllsSkin)
colnames(DataControllsSkin)=c("SAMPLE")
DataControllsSkin$SAMPLE=make.names(as.vector(as.character(DataControllsSkin$SAMPLE)), unique=TRUE)
for(i in 1:length(controlls)){  
  #for(i in 1:3){
  f=controlls[i] 
  #print(f)
  #f=controlls[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_1.csv","",gene)
  print(gene)
  CTdata=read.csv(file=f, na.strings = "")
  CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
  CTdata$SAMPLE=as.character(CTdata$SAMPLE)
  CTdata[,2]=as.numeric(as.character(CTdata$CT))
  df=CTdata[,1:2]
  df$SAMPLE=make.names(as.vector(as.character(df$SAMPLE)), unique=TRUE)
  colnames(df)[2]=gene 
  DataControllsSkin=merge(DataControllsSkin,df, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  #print(DataControllsSkin)
}
DataControllsSkin_nd=DataControllsSkin[c(1:69,73:75),]
write.table(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataControllsSkin_nd.txt",DataControllsSkin_nd, sep="\t", quote=F )      
save(DataControllsSkin_nd, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/DataControllsSkin_nd.Rdata")

