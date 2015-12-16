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

rowsC=c("C001","C002","C003","C004","C005","C006","C007","C008","C009","C010","C011","C012","C013","C014","C015","C016","C017","C018","C019","C020","C021","C022","C023","C024","C025")
DataControllsSkin=rowsC
DataControllsSkin=as.data.frame(DataControllsSkin)
colnames(DataControllsSkin)=c("SAMPLE")

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
  df=CTdata
  df[,1]=as.character(df[,1])
  df[,2]=as.character(df[,2])
  df[,2]=gsub("Undetermined","NA", df[,2],)
  selected.patient=c()
  selected.expression=c()
  #check if there are more than one NAs. selected.patients has patients ID that have at least 2 normal samples.
  for (i in seq(1, length(df[,2]), 3)){  
    patient=c(df[i,1],df[i+1,1],df[i+2,1])
    expression=c(df[i,2],df[i+1,2],df[i+2,2])
    expression=as.numeric(expression)
    print(patient)
    print(expression)
    s=sum(is.na(expression))
    print(s)
    if (sum(is.na(expression)) < 2 ){
      selected.patient=c(selected.patient, patient)
      print("selected.patient")
      print(selected.patient)
      #selected.expression=c(selected.expression, expression)
    }   
  } 
  #filter out the patients that have 2 out of 3 NAs in the samples.
  df=df[df[,1]%in%selected.patient,]
  df[,6]=as.numeric(as.character(df[,6]))
  df[,2]=as.numeric(as.character(df[,2]))
  #check if there are rows in the df
  if(nrow(df)>0){
    print("nrows > 0")
    print(selected.patient)
    #print(head(df))
    dft=t(df)
    dfr=as.data.frame(dft)
    dfr=t(dfr[c(F,F,T)])
    rows=as.character(dfr[,1])
    #filter out patients based on the high varience between the samples of one patient.
    
    #filter the data based on mean +2 std
    if(nrow(dfr)>1){
      d=as.data.frame(dfr)
      d[,6]=as.numeric(as.character(d[, 6]))
      dfr.f=d[d[,6] > quantile(d[,6], .25) - 1.5*IQR(d[,6]) & d[,6] < quantile(d[,6], .75) + 1.5*IQR(d[,6]),]
      print("dfr.f dim")
      print(dim(dfr.f))
      print("outliers")  
      } else{
      if(nrow(dfr)==1){
        dfr.f=dfr
        dfr.f=(as.data.frame(dfr.f)) 
        dfr.f[,15]=as.numeric(as.character(dfr.f[,15]))
      }
    }
    dfr.f[,15]=as.numeric(as.character(dfr.f[,15]))
    dfr.f=dfr.f[dfr.f[,15]< 0.2,]# keep samples with SEM<20%
    print("dfr.f")
    print(dfr.f)
    geneExpr=as.data.frame(dfr.f[,c(1,14)])
  } else{
    if(nrow(df)==0){
      geneExpr=as.data.frame(df[,c(1,14)])
    }
  }
  colnames(geneExpr)[2]=gene 
  print("new geneEXP")
  print(geneExpr)
  DataControllsSkin=merge(DataControllsSkin,geneExpr, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  #print(DataControllsSkin)
}
DataControllsSkin=DataControllsSkin[c(1:23,25),]
write.table(file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.control_FR.txt",DataControllsSkin, sep="\t", quote=F )      
save(DataControllsSkin, file="~/Documents/KAI/cytoskin/SKIN_filter_relaxed/data.control_FR.Rdata")

