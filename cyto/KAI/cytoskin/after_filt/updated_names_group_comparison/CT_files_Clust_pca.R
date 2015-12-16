#This script reads in the files computed in DiffExpPVC_removedC024.R 
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
setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/updated_names_PL_PLN/")
#controlls
DataControllsSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.control_newfilt.txt", sep="\t", header=T)
###Psofiasis#
DataPsoriasisSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.psoriasis_newfilt.txt", sep="\t", header=T )      
###Vitiliigo#
DataVitiliigoSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.vitiliigo_newfilt.txt", sep="\t", header=T )      

###PCA plot of the data

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
dim(data.all)
summary(data.all)

#data.pca=scale(data.all)
data.log=log2(data.all)
data.pca=data.log
library(stringr)
colnames(data.pca)=str_replace(colnames(data.pca),"PH","PL")#sick
colnames(data.pca)=str_replace(colnames(data.pca),"PT","PNL")#healthy
colnames(data.pca)=str_replace(colnames(data.pca),"VH","VL")#sick
colnames(data.pca)=str_replace(colnames(data.pca),"VT","VNL")#healthy
rownames(data.pca)=gsub("IL8b","CXCL8", rownames(data.pca))
save(data.pca, file="data.pca.Rdata")



pheatmap(data.pca, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL, PNL, VL and VNL groups",
         filename="data_log_NAs.png")
#pheatmap(data.pca, cluster_rows = T, cluster_cols = F,scale = "row",
# main="All data scaled, contains NA")

write.table(data.pca,file="data_all_log2_with_NA.txt",sep="\t",row.names=T,quote=FALSE)

data.pca1=data.pca[,! apply( data.pca , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data.pca1, cluster_rows = F, cluster_cols = F,scale = "none", main="data")#data.pca=data.pca[!rowSums(!is.finite(data.pca)),]#remove nas

# remove rows with more than 50 NAs
numNAs <- apply(data.pca1, 1, function(z) sum(is.na(z)))
data.pca.new=data.pca1[!(numNAs > 0.5*length(colnames(data.pca1))),]
pheatmap(data.pca.new, cluster_rows = F, cluster_cols = F, scale = "row", main="data")
dim(data.pca.new)#35 117
#remove cols with more than 50% NAs
data.pca=t(data.pca.new)
numNAs <- apply(data.pca, 1, function(z) sum(is.na(z)))
data.pca.new=data.pca[!(numNAs >0.5*length(colnames(data.pca))),]#previously it was 9
data.pca=t(data.pca.new)
pheatmap(data.pca, cluster_rows = F, cluster_cols = F,scale = "none",main="Gene expression. Rows and columns containing >50% missing values were removed",
         filename="data_log_filt_NAs.png")

dim(data.pca)#35 117
data_log2_filt=data.pca
save(data_log2_filt,file = "data_log2_filt.RData")

#setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison")
#rownames(data_log2_filt)=gsub("IL8b","CXCL8", rownames(data_log2_filt))
#Impute the missing data using k nearest neighbours all data together
library(impute)
library(RColorBrewer)
data=data_log2_filt
data.pca.imp=impute.knn(data ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=36243606)
save(data.pca.imp,file="imputed_data_knn4.RData")
##plot the data distribution
#M=as.matrix(data.pca.imp)
#corm = cor(t(M))
#plot(density(corm))
##

data.pca.imp=as.data.frame(data.pca.imp$data)
pheatmap(data.pca.imp, cluster_rows = T, cluster_cols = F,clustering_distance_rows = "correlation",
clustering_distance_cols = "correlation",clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="none",
main="Filtered log-transformed data after imputation. KNN, k=4", filename="filtered_imputed_knn4.png")
##check the difference between euclidean with NA and without NA
#pheatmap(data.pca.imp, cluster_rows = T, cluster_cols = F,clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="row", main="All data after imputation. KNN, k=4")

library(amap)

hc.data=hcluster(data.pca.imp, method = "correlation", diag = FALSE, upper = FALSE,
                 link = "complete", members = NULL, nbproc = 2,
                 doubleprecision = TRUE)

plot(hc.data)
ct=cutree(hc.data, k=6)

save(hc.data,ct,file="hclust.RData")
sort(ct)
hc.cl.1=as.character(names(sort(ct[ct==1])))
hc.cl.2=as.character(names(sort(ct[ct==2])))
hc.cl.3=as.character(names(sort(ct[ct==3])))
hc.cl.4=as.character(names(sort(ct[ct==4])))
hc.cl.5=as.character(names(sort(ct[ct==5])))
hc.cl.6=as.character(names(sort(ct[ct==6])))
#GO summaries for HC

library(GOsummaries)

g1=hc.cl.1
g2=hc.cl.2
g3=hc.cl.3
g4=hc.cl.4
g5=hc.cl.5
g6=hc.cl.6

gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5, Cluster6 = g6)
gs = gosummaries(gl,  go_branches = c("BP","CC", "ke", "re"), ordered_query = T,  min_set_size = 2, max_set_size = 10000)
plot(gs, fontsize = 8,  term_length=150)
###add expression. change the formulas from km###
data.pca.imp=cbind(rownames(data.pca.imp), data.pca.imp)
colnames(data.pca.imp)[1]="rn"
data_clust1 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.1,],1)#add column cluster number
data_clust1=data_clust1[, c(-1)]
colnames(data_clust1)[118]="Cluster"
data_clust2 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.2,],2)
data_clust2=data_clust2[, c(-1)]
colnames(data_clust2)[118]="Cluster"
data_clust3 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.3,],3)
data_clust3=data_clust3[, c(-1)]
colnames(data_clust3)[118]="Cluster"
data_clust4 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.4,],4)
data_clust4=data_clust4[, c(-1)]
colnames(data_clust4)[118]="Cluster"
data_clust5 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.5,],5)
data_clust5=data_clust5[, c(-1)]
colnames(data_clust5)[118]="Cluster"
data_clust6 <- cbind(data.pca.imp[data.pca.imp$rn%in%hc.cl.6,],6)
data_clust6=data_clust6[, c(-1)]
colnames(data_clust6)[118]="Cluster"
data_clust=rbind(data_clust1,data_clust2,data_clust3,data_clust4,data_clust5,data_clust6)
o<- order(data_clust[,118])###sort by cluster number
data_clust <- data_clust[o, ]###sort by cluster number
head(data_clust)
data_hc=data_clust
library(pheatmap)
library(plyr)
colnames(data_clust)[118]="Cluster"###add column name
#Go summaries
library(stringr)
Groups=colnames(data_clust)
Groups=str_replace(Groups,"C.*","Controls")
Groups=str_replace(Groups,"PL.*","Psoriasis lesional")
Groups=str_replace(Groups,"PNL.*","Psoriasis nonlesional")
Groups=str_replace(Groups,"VL.*","Vitiligo lesional")
Groups=str_replace(Groups,"VNL.*","Vitiligo nonlesional")
annotation=cbind(Groups,Groups )
rownames(annotation)=colnames(data_clust)
annot=as.data.frame(annotation[1:117,1])
colnames(annot)="Groups"

Groups=Groups[1:117]
data_clust=data_clust[,1:117]
#colnames(data_clust)=Groups
head(data_clust)
#################################################

gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)
#setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/")
plot(gs_exp, fontsize = 8, classes = "Groups", term_length=150, filename = "Gosummaries_data_log2_knn4_hc.pdf")

#save(ct,hc.data,data_hc, file = "~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/hc.RData")
save(ct,hc.data,data_hc, file = "hc.RData")

#find significant clusters
library(pvclust)
result=pvclust(t(data.pca.imp[,2:118]), method.hclust="complete",
               method.dist="correlation",
               nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE)
plot(result)
pvrect(result, alpha=0.95)
save(result,file="pvclust.RData")

mydata= data.pca.imp[,2:118]
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

data.pca.imp=data.pca.imp[,c(-1)]
library(amap)
s=set.seed(1234)
km=Kmeans(data.pca.imp,6,iter.max = 100500, nstart =100,
          method = "correlation")

# if needed to compare load(file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/km6.RData")
head(km)
data_clust <- cbind(data.pca.imp,km$cluster)#add column cluster number
o<- order(data_clust[,118])###sort by cluster number
data_clust <- data_clust[o, ]###sort by cluster number
head(data_clust)
data_km=data_clust
library(pheatmap)
library(plyr)
colnames(data_clust)[118]="Cluster"###add column name
#Go summaries
library(stringr)
Groups=colnames(data_clust)
Groups=str_replace(Groups,"C.*","Controls")
Groups=str_replace(Groups,"PL.*","Psoriasis lesional")
Groups=str_replace(Groups,"PNL.*","Psoriasis nonlesional")
Groups=str_replace(Groups,"VL.*","Vitiligo lesional")
Groups=str_replace(Groups,"VNL.*","Vitiligo nonlesional")
annotation=cbind(Groups,Groups )
rownames(annotation)=colnames(data_clust)
annot=as.data.frame(annotation[1:117,1])
colnames(annot)="Groups"

Groups=Groups[1:117]

data_clust=data_clust[,1:117]
#colnames(data_clust)=Groups
head(data_clust)

library(GOsummaries)
#gs_kmeans = gosummaries(km, components = 1:2, exp = data_clust, annotation = annotation)
#plot(gs_kmeans, fontsize = 8, classes = "Groups", filename = "figure3.pdf")
list=sort(km$cluster)
g1=names(list[list==1])
g2=names(list[list==2])
g3=names(list[list==3])
g4=names(list[list==4])
g5=names(list[list==5])
g6=names(list[list==6])
gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5, Cluster6 = g6)
gs = gosummaries(gl,go_branches = c("BP","CC", "ke", "re"), ordered_query = T,  min_set_size = 3, max_set_size = 10000)
plot(gs, fontsize = 8,term_length=150)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)
#setwd("~/Documents/KAI/cytoskin/after_filt/figure_clustering_log2_data/")
#plot(gs_exp, fontsize = 10, term_length=50, classes = "Groups", filename = "~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/Gosummaries_kmean6_knn4.pdf")
plot(gs_exp, fontsize = 10, term_length=50, classes = "Groups", filename = "Gosummaries_kmean6_knn4.pdf")
#
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
#Importance of components:
#  PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10    PC11    PC12    PC13    PC14
#Standard deviation     3.637 2.3725 1.58903 1.39216 1.19033 1.12655 1.02180 0.91926 0.88332 0.80445 0.74400 0.71852 0.70203 0.64975
#Proportion of Variance 0.378 0.1608 0.07214 0.05537 0.04048 0.03626 0.02983 0.02414 0.02229 0.01849 0.01582 0.01475 0.01408 0.01206
#Cumulative Proportion  0.378 0.5388 0.61093 0.66631 0.70679 0.74305 0.77288 0.79703 0.81932 0.83781 0.85362 0.86837 0.88246 0.89452
#PC15   PC16    PC17    PC18    PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27    PC28
#Standard deviation     0.63787 0.6005 0.57894 0.54318 0.52421 0.51367 0.49475 0.45871 0.43199 0.41478 0.40824 0.38189 0.37045 0.33238
#Proportion of Variance 0.01163 0.0103 0.00958 0.00843 0.00785 0.00754 0.00699 0.00601 0.00533 0.00492 0.00476 0.00417 0.00392 0.00316
#Cumulative Proportion  0.90614 0.9164 0.92602 0.93445 0.94230 0.94984 0.95683 0.96285 0.96818 0.97309 0.97786 0.98202 0.98594 0.98910
#PC29    PC30    PC31    PC32    PC33    PC34      PC35
#Standard deviation     0.31044 0.28215 0.26260 0.25282 0.22023 0.15541 1.615e-15
#Proportion of Variance 0.00275 0.00227 0.00197 0.00183 0.00139 0.00069 0.000e+00
#Cumulative Proportion  0.99185 0.99413 0.99610 0.99792 0.99931 1.00000 1.000e+00

#plot(pca.all$x[,1],pca.all$x[,2])
#plot(pca.all$x[,1],pca.all$x[,3])
####Plot 3first PCs and legend

#libraries used
library(ggplot2)
library(stringr)
library(gridBase)
library(gridExtra)
#changed the group names to 

#PS could be PSL (psoriasis lesional), 
#PH could be PSNL (psoriasis nonlesional)
#instead of VS then VL (vitiligo lesional), 
#VH accordingly VNL (vitiligo nonlesional)

Groups=rownames(na.omit(data.pca.knn))
Groups=str_replace(Groups,"C.*","Controls")
Groups=str_replace(Groups,"PL.*","Psoriasis lesional")
Groups=str_replace(Groups,"PNL.*","Psoriasis nonlesional")
Groups=str_replace(Groups,"VL.*","Vitiligo lesional")
Groups=str_replace(Groups,"VNL.*","Vitiligo nonlesional")

Samples=rownames(na.omit(data.pca.knn))  
a=qplot(pca.all$x[,1],pca.all$x[,2],xlab="Principal Component 1",ylab= "Principal Component 2",  colour = Groups, label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))+theme_bw()
a1=a+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")

b=qplot(pca.all$x[,1],pca.all$x[,3],xlab="Principal Component 1",ylab= "Principal Component 3",   colour = Groups, label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
b1=b+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")
c=qplot(pca.all$x[,2],pca.all$x[,3],xlab="Principal Component 2",ylab= "Principal Component 3",  colour = Groups, label=Samples) +scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
c1=c+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")

##install grid extra!!!
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(a)
grid.arrange(a1,b1,c1, legend, 
             ncol=2, nrow=2, widths=c(1/2,1/2))

save(pca.all,Groups,Samples,file="pca.RData")
#save(data.log,Groups,Samples,file="~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/pca.RData")

 ##Plot without lables

a=qplot(pca.all$x[,1],pca.all$x[,2],xlab="Principal Component 1",ylab= "Principal Component 2",  colour = Groups)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))+theme_bw()
a1=a+theme_bw()+theme(legend.position="none")+ stat_ellipse()

b=qplot(pca.all$x[,1],pca.all$x[,3],xlab="Principal Component 1",ylab= "Principal Component 3",   colour = Groups)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
b1=b+theme_bw()+theme(legend.position="none")+ stat_ellipse()
c=qplot(pca.all$x[,2],pca.all$x[,3],xlab="Principal Component 2",ylab= "Principal Component 3",  colour = Groups) +scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
c1=c +theme_bw()+ theme(legend.position="none")+ stat_ellipse()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(a)
grid.arrange(a1,b1,c1, legend, 
             ncol=2, nrow=2, widths=c(1/2,1/2))
####elipse for pso only 

a=qplot(pca.all$x[,1],pca.all$x[,2],xlab="Principal Component 1",ylab= "Principal Component 2",  colour = Groups)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))+theme_bw()
a1=a+theme_bw()+theme(legend.position="none")+ stat_ellipse()

b=qplot(pca.all$x[,1],pca.all$x[,3],xlab="Principal Component 1",ylab= "Principal Component 3",   colour = Groups)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
b1=b+theme_bw()+theme(legend.position="none")+ stat_ellipse()
c=qplot(pca.all$x[,2],pca.all$x[,3],xlab="Principal Component 2",ylab= "Principal Component 3",  colour = Groups) +scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
c1=c +theme_bw()+ theme(legend.position="none")+ stat_ellipse()


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(a)
grid.arrange(a1,b1,c1, legend, 
             ncol=2, nrow=2, widths=c(1/2,1/2))

#CHECK!!!!###Genes importance for each of the components
####Genes with larger positive and negative weights have more impact on the component
genes=rownames(pca.all$rotation)

plot_rotation1 <- cbind(pca.all$rotation[,1],genes)
genes_pc1=sort(abs(pca.all$rotation[,1]), decreasing=T)#take the names of first top 10% of genes=genes that separate samples the most.
sep_genes_pc1=names(genes_pc1[1:ceiling(0.2*length(genes))])
head(plot_rotation1)
dim(plot_rotation1)
colnames(plot_rotation1)[1]="Weights"
tmp1 <- plot_rotation1[,1]#for the 1st component
dd1 <- data.frame(y=genes, x=tmp1)
dd1$x = as.numeric(as.character(dd1$x))
w1=qplot(y, x,data=dd1, geom="bar",fill = I("grey50"),width=.5, stat="identity") + coord_flip()+labs(y="Rotation matrix weights for the first principal component", x="Gene name")+theme_bw()+theme (axis.text.y = element_text(colour="grey20",size=8,face="plain"),axis.title.x = element_text(size=10,face="plain"))


plot_rotation2 <- cbind(pca.all$rotation[,2],genes)
genes_pc2=sort(abs(pca.all$rotation[,2]), decreasing=T)#take the names of first top 20% of genes=genes that separate samples the most.
sep_genes_pc2=names(genes_pc2[1:ceiling(0.2*length(genes))])
colnames(plot_rotation2)[1]="Weights"
tmp2=plot_rotation2[,1]##for the 2nd component
dd2 <- data.frame(y=genes, x=tmp2)
dd2$x = as.numeric(as.character(dd2$x))
w2=qplot(y, x,data=dd2, geom="bar",fill = I("grey50"),width=.5, stat="identity") + coord_flip()+labs(y="Rotation matrix weights for the second principal component", x="Gene name")+theme_bw()+theme (axis.text.y = element_text(colour="grey20",size=8,face="plain"),axis.title.x = element_text(size=10,face="plain"))



plot_rotation3 <- cbind(pca.all$rotation[,3],genes)
genes_pc3=sort(abs(pca.all$rotation[,3]), decreasing=T)#take the names of first top 20% of genes=genes that separate samples the most.
sep_genes_pc3=names(genes_pc3[1:ceiling(0.2*length(genes))])
colnames(plot_rotation3)[1]="Weights"
tmp3=plot_rotation3[,1]##for the 3rd component
dd3<- data.frame(y=genes, x=tmp3)
dd3$x = as.numeric(as.character(dd3$x))
w3=qplot(y, x,data=dd3, geom="bar",fill = I("grey50"),width=.5, stat="identity") + coord_flip()+labs(y="Rotation matrix weights for the third principal component", x="Gene name")+theme_bw()+theme(axis.text.y = element_text(colour="grey20",size=8,face="plain"),axis.title.x = element_text(size=10,face="plain"))

l=textGrob("Genes with the largest positive\nor negative weights have stronger\nimpact on the corresponding component ")
#width=stringWidth(string), height=stringHeight(string)

#grid.arrange(w1,w2,w3,l, ncol=2, nrow=2, widths=c(1/2,1/2))

grid.arrange(w1,w2,w3, ncol=2, nrow=2, widths=c(1/2,1/2))




