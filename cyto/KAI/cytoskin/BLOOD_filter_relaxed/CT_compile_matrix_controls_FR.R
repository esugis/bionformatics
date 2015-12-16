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
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_blood2", full.names = TRUE)
#take corresponding _1,_2,_3, _4, _5, _6 for further analysis
controlls1=csvFiles[grepl("*_1.csv", csvFiles)]
controlls2=csvFiles[grepl("*_2.csv", csvFiles)]
psoriasis1=csvFiles[grepl("*_3.csv", csvFiles)]
psoriasis2=csvFiles[grepl("*_4.csv", csvFiles)]
vitiliigo1=csvFiles[grepl("*_5.csv", csvFiles)]  
vitiliigo2=csvFiles[grepl("*_6.csv", csvFiles)]  
# Assemble the matrix for each of the groups. Combine controlls1,2, psoriasis 1,2, vitiligo 1,2
#read in row names
CTdata1=read.csv(file=controlls1[1], na.strings = "")
CTdata1=CTdata1[7:114,2:ncol(CTdata1)]
CTdata2=read.csv(file=controlls2[1], na.strings = "")
CTdata2=CTdata2[7:123,2:ncol(CTdata2)]
rn1=as.character(as.vector(CTdata1[,1]))#sample names from the first csv
rn2=as.character(as.vector(CTdata2[,1]))#sample names from the second csv
rowsC=unique(c(rn1, rn2))
length(rowsC)
#[1] 75
DataControllsBlood=rowsC
DataControllsBlood=as.data.frame(DataControllsBlood)
colnames(DataControllsBlood)=c("SAMPLE")
library(functional)
#for each of controlls
#selected.patient=c()
#selected.expression=c()
for(i in 1:length(controlls1)){  
  #for(i in 3:5){
  f1=controlls1[i] 
  f2=controlls2[i]
  #f1=controlls1[8]
  #f2=controlls2[8]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_blood2/CtkalkulaatorPBMC_","",f1)
  gene=gsub("_B2mg_1.csv","",gene)
 
  print(gene)
  CTdata1=read.csv(file=f1, na.strings = "")
  #CTdata1$SAMPLE=as.character(CTdata1$SAMPLE)
 # CTdata1[,2:17]=apply(CTdata1[,2:17],2, as.numeric)
  CTdata1=CTdata1[7:114,2:16]
  CTdata2=read.csv(file=f2, na.strings = "")
  CTdata2=CTdata2[7:123,2:16]
  CTdata=rbind(CTdata1, CTdata2)
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
    if (sum(is.na(expression)) < 3 ){# here is the change from 2 to 3 NA samples
      selected.patient=c(selected.patient, patient)
      print("selected.patient")
      print(selected.patient)
      #selected.expression=c(selected.expression, expression)
    }   
  } 
  #filter out the patients that have 2 out of 3 NAs in the samples.
  df=df[df[,1]%in%selected.patient,]
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
      ##usingMAD
      #filt=median(as.numeric(dfr[,6]),na.rm=T)+2*mad(as.numeric(dfr[,6]),na.rm=T)#filter based on the reference gene cycles
      #dfr.f=dfr[dfr[,6]<filt,]#keep samples with reference genes number of cycles< meadian+2mad
      
      ##alternatively filter on quartiles
      d=as.data.frame(dfr)
      d[,8]=as.numeric(as.character(d[, 8]))
      dfr.f=d[d[,8] > quantile(d[,8], .25) - 1.5*IQR(d[,8]) & d[,8] < quantile(d[,8], .75) + 1.5*IQR(d[,8]),]
      print("dfr.f dim")
      print(dim(dfr.f))
      print("outliers")
      #print(dfr[dfr[,6]>filt,])
    } else{
      if(nrow(dfr)==1){
        dfr.f=dfr
        dfr.f=(as.data.frame(dfr.f)) 
        dfr.f[,15]=as.character(dfr.f[,15])
      }
    }
    dfr.f[,15]=as.numeric(as.character(dfr.f[, 15]))
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
  #print("new geneEXP")
  #print(geneExpr)
  DataControllsBlood=merge(DataControllsBlood,geneExpr, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  #print(DataControllsBlood)
}
##using MAD
#DataControllsBlood=DataControllsBlood[c(1:48,52:75),]
#save(DataControllsBlood, file="~/Documents/KAI/cytoskin/BLOOD/datacontrol.RData")
#write.table(file="~/Documents/KAI/cytoskin/BLOOD/datacontrol.txt",DataControllsBlood, sep="\t", quote=F )      
##using Quqrtiles
DataControllsBlood_QR=DataControllsBlood[c(1:48,52:75),]
dim(DataControllsBlood_QR)
save(DataControllsBlood_QR, file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/datacontrol_QR.RData")
write.table(file="~/Documents/KAI/cytoskin/BLOOD_filter_relaxed/datacontrol_QR.txt",DataControllsBlood_QR, sep="\t", quote=F )      


