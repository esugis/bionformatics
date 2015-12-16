#This script reads csv files that contain info about expression of controlls, psoriasis and vitiligo
#for each of the genes.
#Files are in the format gene_names_1.csv, gene_names_2.csv, gene_names_3.csv
#where: 
#*_1.csv corresponds to controls
#*_2.csv psoriasis
#*_3.csv vitiliigo
#Csv are obtained from xlsx provided by the collaborators.
#Script
#-Reads the files,
#-Proceeses the data(filtering based on 2 or more of the parralell experiments are undefined, 
#refenece gene is not stable or had very low expression,SEM should be < 0,2)
#-Assembles data set from the filtered data.
#-Scales the data(convert to Z scores)
#-FIlteres out rows and cols with >50%missing
#-Imputes missing values
#-Does K-means clustering and finds GO annotations for each cluster.
#-Does herarchical clustering and GO annotations for selected groups.
# Additionally using pvclust() finds the clusters that are strongly supported by the data.

library(xlsx)
require(gdata)
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
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
#rms= function(x){
#  rms = rowMeans(x, na.rm = TRUE)
#  rowSums(is.na(mydf)) <= 3
#}
#1. Rename the files in the directory
#oldNames = list.files("~/Documents/KAI/cytoskin/CT_skin", full.names = TRUE)
#newNames=gsub("\\(","_",oldNames)
#newNames=gsub("\\)","",newNames)
#file.rename(from = file.path(oldNames), to = file.path( newNames))

#for(i in 1:length(newNames)){ # for all the files in the folder
##for(i in 1:1){ 
#file=newNames[i]
#f=gsub(".xlsx","",file)
#for(j in 1:3){
#system(sprintf("xlsx2csv -s %s %s %s_%s.csv",j,file,f,j));#split xlsx to 3 csv
#where *_1.csv corresponds to controls
#*_2.csv psoriasis
#*_3.csv vitiliigo
#}
#}
#list all files with csv from the folder
csvFiles = list.files("~/Documents/KAI/cytoskin/CT_skin", full.names = TRUE)
#take corresponding _1,_2,_3 for further analysis
controlls=csvFiles[grepl("*_1.csv", csvFiles)]
psoriasis=csvFiles[grepl("*_2.csv", csvFiles)]
vitiliigo=csvFiles[grepl("*_3.csv", csvFiles)]  

rowsC=c("C001","C002","C003","C004","C005","C006","C007","C008","C009","C010","C011","C012","C013","C014","C015","C016","C017","C018","C019","C020","C021","C022","C023","C024","C025")
DataControllsSkin=rowsC
DataControllsSkin=as.data.frame(DataControllsSkin)
colnames(DataControllsSkin)=c("SAMPLE")
library(functional)
#for each of controlls
#selected.patient=c()
#selected.expression=c()
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
       filt=median(as.numeric(dfr[,6]),na.rm=T)+2*sd(as.numeric(dfr[,6]),na.rm=T)#filter based on the reference gene cycles
       dfr.f=dfr[dfr[,6]<filt,]#keep samples with reference genes number of cycles< meadian+2std
      } else{
             if(nrow(dfr)==1){
             dfr.f=dfr
             dfr.f=(as.data.frame(dfr.f)) 
             dfr.f[,15]=as.character(dfr.f[,15])
             }
        }
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
write.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.control_newfilt.txt",DataControllsSkin, sep="\t", quote=F )      

###Psofiasis#

data.cols=read.csv(file=psoriasis[2])
data.cols=data.cols[7:nrow(data.cols),]
rowsP=as.character(data.cols[,2])
rowsP=rowsP[c(F,F,T)]
DataPsoriasisSkin=rowsP
DataPsoriasisSkin=as.data.frame(DataPsoriasisSkin)
colnames(DataPsoriasisSkin)=c("SAMPLE")
library(functional)
#for each in psoriasis
for(i in 1:length(psoriasis)){  
  #for(i in 1:3){
  f=psoriasis[i] 
  #print(f)
  #f=psoriasis[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_2.csv","",gene)
  CTdata=read.csv(file=f)
  CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
  df=CTdata
  df[,1]=as.character(df[,1])
  df[,2]=as.character(df[,2])
  df[,2]=gsub("Undetermined","NA", df[,2])
  x=df[,2]
  z <- gsub("\\s+", "", x)
  x[z==""] <- NA 
  df[,2]=x
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
      filt=median(as.numeric(dfr[,6]),na.rm=T)+2*sd(as.numeric(dfr[,6]),na.rm=T)#filter based on the reference gene cycles
      dfr.f=dfr[dfr[,6]<filt,]#keep samples with reference genes number of cycles< meadian+2std
    } else{
      if(nrow(dfr)==1){
        dfr.f=dfr
        dfr.f=(as.data.frame(dfr.f)) 
        dfr.f[,15]=as.character(dfr.f[,15])
      }
    }
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
  DataPsoriasisSkin=merge(DataPsoriasisSkin,geneExpr, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
  ##print(DataControllsSkin)
}
write.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.psoriasis_newfilt.txt",DataPsoriasisSkin, sep="\t", quote=F )      


###Vitiliigo#

data.cols=read.csv(file=vitiliigo[1])
data.cols=data.cols[7:nrow(data.cols),]
rowsP=as.character(data.cols[,2])
rowsP=rowsP[c(F,F,T)]
DataVitiliigoSkin=rowsP
DataVitiliigoSkin=as.data.frame(DataVitiliigoSkin)
colnames(DataVitiliigoSkin)=c("SAMPLE")
library(functional)
#for each in vitiliigo
for(i in 1:length(vitiliigo)){  
  #for(i in 1:3){
  f=vitiliigo[i] 
  #print(f)
  #f=vitiliigo[1]
  gene=gsub("/Users/nikolaeva/Documents/KAI/cytoskin/CT_skin/Ctkalkulaatornahk_","",f)
  gene=gsub("_3.csv","",gene)
  CTdata=read.csv(file=f)
  CTdata=CTdata[7:nrow(CTdata),2:ncol(CTdata)]
  df=CTdata
  df[,1]=as.character(df[,1])
  df[,2]=as.character(df[,2])
  df[,2]=gsub("Undetermined","NA", df[,2])
  x=df[,2]
  z <- gsub("\\s+", "", x)
  x[z==""] <- NA 
  df[,2]=x
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
      filt=median(as.numeric(dfr[,6]),na.rm=T)+2*sd(as.numeric(dfr[,6]),na.rm=T)#filter based on the reference gene cycles
      dfr.f=dfr[dfr[,6]<filt,]#keep samples with reference genes number of cycles< meadian+2std
    } else{
      if(nrow(dfr)==1){
        dfr.f=dfr
        dfr.f=(as.data.frame(dfr.f)) 
        dfr.f[,15]=as.character(dfr.f[,15])
      }
    }
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
  
  DataVitiliigoSkin=merge(DataVitiliigoSkin,geneExpr, by.x="SAMPLE", by.y="SAMPLE", all.x=T )
 # print(DataVitiliigoSkin)
}
write.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.vitiliigo_newfilt.txt",DataVitiliigoSkin, sep="\t", quote=F )     
###PCA plot of the dat
ctrls=DataControllsSkin
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:43]
ctrls <- apply(ctrls[,2:43],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisSkin)# psoriasis
rows=rownames(pso[2:43,])
cols=as.character(pso[1,])
ps <- apply(pso[2:43,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiliigoSkin)# vitiliigo
rows=rownames(vit[2:43,])
cols=as.character(vit[1,])
vi <- apply(vit[2:43,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

library(pheatmap)
data.all=cbind(ctrls, ps,vi)
summary(data.all)

#data.pca=scale(data.all)
data.log=log2(data.all)
data.pca=data.log
#data.pca=scale(data.log)#col scaling
#data.pca=scale_rows(data.log)

pheatmap(data.pca, cluster_rows = F, cluster_cols = F,scale = "row", main="All data scaled, contains NA",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log_NAs.png")
#pheatmap(data.pca, cluster_rows = T, cluster_cols = F,scale = "row",
# main="All data scaled, contains NA")

write.table(data.pca,file="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_all_log2_with_NA.txt",sep="\t",row.names=T,quote=FALSE)

data.pca1=data.pca[,! apply( data.pca , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data.pca1, cluster_rows = F, cluster_cols = F,scale = "row", main="data")#data.pca=data.pca[!rowSums(!is.finite(data.pca)),]#remove nas

# remove rows with more than 50 NAs
numNAs <- apply(data.pca1, 1, function(z) sum(is.na(z)))
data.pca.new=data.pca1[!(numNAs > 0.5*length(colnames(data.pca1))),]
pheatmap(data.pca.new, cluster_rows = F, cluster_cols = F, scale = "row", main="data")
dim(data.pca.new)
#remove cols with more than 50% NAs
data.pca=t(data.pca.new)
numNAs <- apply(data.pca, 1, function(z) sum(is.na(z)))
data.pca.new=data.pca[!(numNAs >0.5*length(colnames(data.pca))),]#previously it was 9
data.pca=t(data.pca.new)
pheatmap(data.pca, cluster_rows = F, cluster_cols = F,scale = "row",main="All data scaled, contains NA",
         filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log_filt_NAs.png")

dim(data.pca)#35 118
data_log2_filt=data.pca
save(data_log2_filt, file = "~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")

#Impute the missing data using k nearest neighbours all data together
library(impute)
library(RColorBrewer)
data=data.pca
data.pca.imp=impute.knn(data ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)

##plot the data distribution
M=as.matrix(data.pca.imp)
corm = cor(t(M))
plot(density(corm))
##

data.pca.imp=as.data.frame(data.pca.imp$data)
pheatmap(data.pca.imp, cluster_rows = T, cluster_cols = F,clustering_distance_rows = "correlation",
clustering_distance_cols = "correlation",clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="row",
main="Filtered data after imputation. KNN, k=4", filename="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/filtered_imputed_knn4.png")
##check the difference between euclidean with NA and without NA
#pheatmap(data.pca.imp, cluster_rows = T, cluster_cols = F,clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="row", main="All data after imputation. KNN, k=4")

library(amap)

hc.data=hcluster(data.pca.imp, method = "correlation", diag = FALSE, upper = FALSE,
                 link = "complete", members = NULL, nbproc = 2,
                 doubleprecision = TRUE)

plot(hc.data)
ct=cutree(hc.data, k=5)
sort(ct)
hc.cl.1=as.character(names(sort(ct[ct==1])))
hc.cl.2=as.character(names(sort(ct[ct==2])))
hc.cl.3=as.character(names(sort(ct[ct==3])))
hc.cl.4=as.character(names(sort(ct[ct==4])))
hc.cl.5=as.character(names(sort(ct[ct==5])))
#hc.cl.6=as.character(names(sort(ct[ct==6])))
#GO summaries for HC

library(GOsummaries)
#gs_kmeans = gosummaries(km, components = 1:2, exp = data_clust, annotation = annotation)
#plot(gs_kmeans, fontsize = 8, classes = "Groups", filename = "figure3.pdf")

g1=hc.cl.1
g2=hc.cl.2
g3=hc.cl.3
g4=hc.cl.4
g5=hc.cl.5
#g6=hc.cl.6

gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5)#, Cluster6 = g6)
gs = gosummaries(gl)
plot(gs, fontsize = 8)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)
setwd("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/")
plot(gs_exp, fontsize = 8, classes = "Groups", filename = "Gosummaries_data_log2_knn4_hc.pdf")

save(ct,hc, file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/noscale/hc.RData")

#find significant clusters
library(pvclust)
result=pvclust(t(data.pca.imp), method.hclust="complete",
               method.dist="correlation",
               nboot=10000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE)
plot(result)
pvrect(result, alpha=0.95)

###Impute each data part separately ctrls, vt, ps
#library(impute)

###impute missing Was performing worse(PCA and pvclust) in finding patterms in the data
#data=data.pca
#data.ctrls=data[,1:24]
#data.ph=data[,25:57]
#data.pt=data[,58:89]
#data.vh=data[,90:105]
#data.vt=data[,106:120]  
  
#ctrls.imp=impute.knn(data.ctrls ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
#ctrls.imp=as.data.frame(ctrls.imp$data)
#summary(ctrls.imp)

#ph.imp=impute.knn(data.ph ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
#ph.imp=as.data.frame(ph.imp$data)

#pt.imp=impute.knn(data.pt ,k = 4, rowmax = 0.6, colmax = 0.8, maxp = 1500, rng.seed=362436069)
#pt.imp=as.data.frame(pt.imp$data)

#vh.imp=impute.knn(data.vh ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
#vh.imp=as.data.frame(vh.imp$data)

#vt.imp=impute.knn(data.vt ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
#vt.imp=as.data.frame(vt.imp$data)
##Combine the imputed parts together

#data.pca.imp=cbind(ctrls.imp,ph.imp,pt.imp,vh.imp,vt.imp)
#pheatmap(data.pca.imp, cluster_rows = T, cluster_cols = F,sclustering_distance_rows = "correlation",clustering_distance_cols = "correlation",clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="row", main="All data after imputation. KNN, k=4")

##Quantile normalization(erases most of the variation in the data,PCA makes no sense)
#library(preprocessCore)
#data.pca.imp.m=data.matrix(data.pca.imp)
#data.qnorm=normalize.quantiles(data.pca.imp.m)
#data.pca.imp=data.qnorm
#k-means on imputed data

####http://www.statmethods.net/advstats/cluster.html
#K-means clustering is the most popular partitioning method. 
#It requires the analyst to specify the number of clusters to extract. 
#A plot of the within groups sum of squares by number of clusters extracted 
#can help determine the appropriate number of clusters.

mydata= data.pca.imp
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


library(amap)
s=set.seed(1234)
km=Kmeans(data.pca.imp,5,iter.max = 100500, nstart =50,
          method = "correlation")
head(km)
data_clust <- cbind(data.pca.imp,km$cluster)#add column cluster number
o<- order(data_clust[,119])###sort by cluster number
data_clust <- data_clust[o, ]###sort by cluster number
head(data_clust)
data_km=data_clust
library(pheatmap)
library(plyr)
colnames(data_clust)[119]="Cluster"###add column name
#Go summaries
library(stringr)
Groups=colnames(data_clust)
Groups=str_replace(Groups,"C.*","C")
Groups=str_replace(Groups,"PH.*","PH")
Groups=str_replace(Groups,"PT.*","PT")
Groups=str_replace(Groups,"VH.*","VH")
Groups=str_replace(Groups,"VT.*","VT")
annotation=cbind(Groups,Groups )
rownames(annotation)=colnames(data_clust)
annot=as.data.frame(annotation[1:118,1])
colnames(annot)="Groups"

Groups=Groups[1:118]
data_clust=data_clust[,1:118]
#colnames(data_clust)=Groups
head(data_clust)

write.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024//data_log2_kmeans4_knn4.txt",data_km, sep="\t", quote=F )      
save(km, file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024//km5.RData")

library(GOsummaries)
#gs_kmeans = gosummaries(km, components = 1:2, exp = data_clust, annotation = annotation)
#plot(gs_kmeans, fontsize = 8, classes = "Groups", filename = "figure3.pdf")
list=sort(km$cluster)
g1=names(list[list==1])
g2=names(list[list==2])
g3=names(list[list==3])
g4=names(list[list==4])
g5=names(list[list==5])
#g6=names(list[list==6])
gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5)#, Cluster6 = g6)
gs = gosummaries(gl)
plot(gs, fontsize = 8)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)

plot(gs_exp, fontsize = 8, classes = "Groups")
##end go summaries


## Calculate the within groups sum of squared error (SSE) for the number of cluster solutions selected by the use
#wss=c()
#for (i in 2:10){
#  wss.single=c()
#  for (j in 1:100){
#    v=sum(Kmeans(data.pca.imp,i,iter.max = 1500, nstart = 100, method = "correlation")$withinss)
#  wss.single=c(wss.single,v)
#}
#w=sum(wss.single)/100
#wss=c(wss,w)
#}
#y=c(2,3,4,5,6,7,8,9,10)
#x=wss
#xy=cbind(y,x)
#g_range=range(2,10)
#plot(xy,type="o", xlab="Number of clusters", ylab="Within Groups Sum of Squares", main="Cluster Solutions", xlim=g_range)
#sum(Kmeans(data.pca.imp,7,iter.max = 1500, nstart = 100, method = "correlation")$withinss)
##2 0.7821383
##3 0.4033993
##4 0.403408
##5  0.2413656
##6 0.1171066
##7 0.303912


#hclust
hc=hcluster(data.pca.imp, method = "correlation", diag = FALSE, upper = FALSE,
         link = "complete", members = NULL, nbproc = 2,
         doubleprecision = TRUE)
#hc=hcluster(data.pca.imp, method = "correlation", diag = FALSE, upper = FALSE,
#            link = "single", members = NULL, nbproc = 2,
#            doubleprecision = TRUE)
plot(hc)
hc1=hcluster(data.all, method = "correlation", diag = FALSE, upper = FALSE,
            link = "complete", members = NULL, nbproc = 2,
            doubleprecision = TRUE)
#hc=hcluster(data.pca.imp, method = "correlation", diag = FALSE, upper = FALSE,
#            link = "single", members = NULL, nbproc = 2,
#            doubleprecision = TRUE)
plot(hc1)



#parPvclust(cl=NULL, data.pca.imp, method.hclust="complete",
#           method.dist="correlation", use.cor="pairwise.complete.obs",
#           nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE,
#           init.rand=TRUE, seed=NULL, iseed=NULL)
##Dynamic Tree Cut

#library(dynamicTreeCut)
#cutreeDynamic(
#  hc, cutHeight = NULL, minClusterSize =3,
  # Basic tree cut options
#  method = "hybrid")

#hc=hclust(d = dissim, method = "average")
#library(moduleColor )

#DetectedColors = NULL;
#DetectedColors = cbind(DetectedColors,
#                       labels2colors(cutreeDynamic(hc, cutHeight = NULL,
#                                                   minClusterSize = 3,
#                                                   method = "hybrid", deepSplit = 3,
#                                                   pamStage = TRUE, maxDistToLabel = 0,
#                                                   verbose = 0)));
#hclustplotn(hc,DetectedColors, RowLabels = Methods, main="");
 

##Mclust on imputed data?
#library(pheatmap)
#library(plyr)
#library(mclust)
#DataMclust=Mclust(data.pca.imp)

#Matrix_mclust<-cbind(data.pca.imp,DataMclust$classification)
#head(Matrix_mclust)
#o<- order(Matrix_mclust[,116 ])###sort by cluster number
#Matrix_mclust <- Matrix_mclust[o, ]
#colnames(Matrix_mclust)[116]="Cluster"
#Mclust_plot=ddply(Matrix_mclust, .(Cluster), function(x) apply(x, 2, mean))
#rownames(Mclust_plot)=Mclust_plot[,116]
#pheatmap(Mclust_plot, scale ="row",cluster_cols = TRUE)
###

##
data.pca.knn=data.pca.imp
#data.pca=data.pca.new
data.pca.knn=t(data.pca.knn)
pca.all=prcomp(data.pca.knn, scale=T, center=T)#, na.action=na.omit)
#pca.all=prcomp(~., data=data.pca, scale=T, center=T, na.action=na.omit)
summary(pca.all)
#knn4 all data scale 


#Importance of components:
#  PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8    PC9    PC10    PC11    PC12    PC13   PC14    PC15    PC16    PC17
#Standard deviation     3.6332 2.3631 1.5951 1.39040 1.18670 1.12397 1.02496 0.91810 0.8815 0.80443 0.74860 0.72263 0.70550 0.6482 0.64711 0.60137 0.58146
#Proportion of Variance 0.3772 0.1596 0.0727 0.05523 0.04024 0.03609 0.03002 0.02408 0.0222 0.01849 0.01601 0.01492 0.01422 0.0120 0.01196 0.01033 0.00966
#Cumulative Proportion  0.3772 0.5367 0.6094 0.66464 0.70488 0.74097 0.77099 0.79507 0.8173 0.83576 0.85178 0.86670 0.88092 0.8929 0.90488 0.91522 0.92488
#PC18    PC19   PC20    PC21   PC22    PC23    PC24    PC25    PC26    PC27    PC28    PC29    PC30   PC31    PC32    PC33    PC34
#Standard deviation     0.54979 0.53643 0.5225 0.49602 0.4621 0.43105 0.41192 0.40781 0.38479 0.37187 0.33191 0.31039 0.28799 0.2647 0.25734 0.22072 0.15514
#Proportion of Variance 0.00864 0.00822 0.0078 0.00703 0.0061 0.00531 0.00485 0.00475 0.00423 0.00395 0.00315 0.00275 0.00237 0.0020 0.00189 0.00139 0.00069
#Cumulative Proportion  0.93351 0.94173 0.9495 0.95656 0.9627 0.96797 0.97282 0.97757 0.98181 0.98576 0.98890 0.99166 0.99403 0.9960 0.99792 0.99931 1.00000
#PC35
#Standard deviation     1.613e-15
#Proportion of Variance 0.000e+00
#Cumulative Proportion  1.000e+00

#plot(pca.all$x[,1],pca.all$x[,2])
#plot(pca.all$x[,1],pca.all$x[,3])
####Plot 3first PCs and legend

#libraries used
library(ggplot2)
library(stringr)
library(gridBase)
library(gridExtra)
Groups=rownames(na.omit(data.pca.knn))
Groups=str_replace(Groups,"C.*","C")
Groups=str_replace(Groups,"PH.*","PH")
Groups=str_replace(Groups,"PT.*","PT")
Groups=str_replace(Groups,"VH.*","VH")
Groups=str_replace(Groups,"VT.*","VT")
Samples=rownames(na.omit(data.pca.knn))  
a=qplot(pca.all$x[,1],pca.all$x[,2],xlab="Principal Component 1",ylab= "Principal Component 2",  colour = Groups, label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
a1=a+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") + theme(legend.position="none")

b=qplot(pca.all$x[,1],pca.all$x[,3],xlab="Principal Component 1",ylab= "Principal Component 3",   colour = Groups, label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
b1=b+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") + theme(legend.position="none")
c=qplot(pca.all$x[,2],pca.all$x[,3],xlab="Principal Component 2",ylab= "Principal Component 3",  colour = Groups, label=Samples) +scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
c1=c+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") + theme(legend.position="none")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(a)
grid.arrange(a1,b1,c1, legend, 
             ncol=2, nrow=2, widths=c(1/2,1/2))

save(pca.all,Groups,Samples,file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/pca.RData")
#save(data.log,Groups,Samples,file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/pca.RData")

 
