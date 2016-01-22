#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
library(plyr)
library(amap)
library(impute)
library(RColorBrewer)
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")

#rownames(data_skin_cp_nd)=gsub("IL8b","CXCL8", rownames(data_skin_cp_nd))
#rownames(data_skin_cp_nd)[16]="INFG"
#save(data_skin_cp_nd, file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

#impute missing data using KNN 
#remove cols with more than 50% NAs
data.pca=t(data.skin)
numNAs <- apply(data.pca, 1, function(z) sum(is.na(z)))
data.pca.new=data.pca[!(numNAs >0.5*length(colnames(data.pca))),]
data.pca=t(data.pca.new)
pheatmap(data.pca, cluster_rows = F, cluster_cols = F,scale = "none",main="Gene expression.Columns containing >50% missing values were removed")

dim(data.pca)#42 77
data_skin_na_rm=data.pca
save(data_skin_na_rm,file = "data_skin_na_rm.RData")

#setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison")
#rownames(data_log2_filt)=gsub("IL8b","CXCL8", rownames(data_log2_filt))
#Impute the missing data using k nearest neighbours all data together

data=data_skin_na_rm
data.pca.imp=impute.knn(data ,k = 4, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=36243606)
data_imp_knn=data.pca.imp
save(data_imp_knn,file="imputed_data_knn4.RData")

data.imp=as.data.frame(data_imp_knn$data)
pheatmap(data.imp, cluster_rows = F, cluster_cols = F,clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="none",
         main="Filtered log-transformed data after imputation. KNN, k=4", filename="filtered_imputed_knn4.png")
##check the difference between euclidean with NA and without NA
#pheatmap(data.imp, cluster_rows = T, cluster_cols = F,clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="row", main="All data after imputation. KNN, k=4")

#####
genes=data.imp
#hc.data=hcluster(genes, method = "spearman", diag = FALSE, upper = FALSE,#correlation=centered pearson 1-corr(xy)
#                 link = "average", members = NULL, nbproc = 2, #Ward's minimum variance method aims at finding compact, spherical clusters.
#                 doubleprecision = TRUE)
#plot(hc.data)
hc.data=hcluster(genes, method = "euclidean", diag = FALSE, upper = FALSE,
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
#gs_kmeans = gosummaries(km, components = 1:2, exp = data_clust, annotation = annotation)
#plot(gs_kmeans, fontsize = 8, classes = "Groups", filename = "figure3.pdf")

g1=hc.cl.1
g2=hc.cl.2
g3=hc.cl.3
g4=hc.cl.4
g5=hc.cl.5
g6=hc.cl.6

gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5, Cluster6 = g6)
gs = gosummaries(gl)
plot(gs, fontsize = 8)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)
plot(gs_exp, fontsize = 8, classes = "Groups", filename = "Gosummaries_data_knn4_hc.pdf")


#strong clusters
hc.cl.1=as.character(c("CXCL1","CXCL2","IL17F","LCN2"))
hc.cl.2=as.character(c("AIM2","CXCL10","IFIH1"))
hc.cl.3=as.character(c("IL1F6","S100A8","S100A9"))
hc.cl.4=as.character(c("CCL2","TNF"))
hc.cl.5=as.character(c("CTLA4","FOXP3"))
hc.cl.6=as.character(c("CASP1", "KLRK1","WIPI1"))
hc.cl.7=as.character(c("CCL20","EOMES","IFNAR1","IFNGR","IL10","IL1b","MICB" ))
#GO summaries for HC

library(GOsummaries)
#gs_kmeans = gosummaries(km, components = 1:2, exp = data_clust, annotation = annotation)
#plot(gs_kmeans, fontsize = 8, classes = "Groups", filename = "figure3.pdf")

g1=hc.cl.1
g2=hc.cl.2
g3=hc.cl.3
g4=hc.cl.4
g5=hc.cl.5
g6=hc.cl.6
g7=hc.cl.7
gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5, Cluster6 = g6,Cluster7 = g7)
gs = gosummaries(gl)
plot(gs, fontsize = 8)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)
plot(gs_exp, fontsize = 8, classes = "Groups", filename = "Gosummaries_data_knn4_hc_strong_clusters.pdf")



#find significant clusters
library(pvclust)
tdata=t(data.imp)
summary(tdata)
result=pvclust(tdata, method.hclust="complete",
               method.dist="euclidean",
               #use.cor="pairwise.complete.obs",
               nboot=10000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE)
plot(result)
pvrect(result, alpha=0.94)
save(result,file="pvclust.RData")


mydata= data.imp
dim(mydata)
wss <- nrow(mydata-1)*sum(apply(mydata,2,var))
for (i in 1:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
ylab="Within groups sum of squares")

s=set.seed(1234)
km=Kmeans(data.imp,6,iter.max = 1000500, nstart =77,
          method = "euclidean")
head(km)
data_clust <- cbind(data.imp,km$cluster)#add column cluster number
o<- order(data_clust[,78])###sort by cluster number
data_clust <- data_clust[o, ]###sort by cluster number
head(data_clust)
data_km=data_clust
library(pheatmap)
library(plyr)
colnames(data_clust)[78]="Cluster"###add column name
#Go summaries
library(stringr)
Groups=colnames(data_clust)
Groups=str_replace(Groups,"C.*","C")
Groups=str_replace(Groups,"PNL.*","PNL")
Groups=str_replace(Groups,"PL.*","PL")

annotation=cbind(Groups,Groups )
rownames(annotation)=colnames(data_clust)
annot=as.data.frame(annotation[1:77,1])
colnames(annot)="Groups"

Groups=Groups[1:77]
data_clust=data_clust[,1:77]
#colnames(data_clust)=Groups
head(data_clust)
genes_clust=cbind(gene=rownames(data_clust),cluster=data_km[,78])
save(genes_clust, file = "genes_clust.RData")
write.table(file="genes_clust.txt",genes_clust, sep="\t", quote=F )
write.table(file="KM6_knn4.txt",data_km, sep="\t", quote=F )      
save(km, file = "km6_knn4.RData")

pheatmap(data.imp, cluster_rows = F, cluster_cols = F,clustering_distance_rows="eucledian",kmeans_k=4,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250), scale="none",
         main="Kmeans 6 ")
list=sort(km$cluster)
g1=names(list[list==1])
g2=names(list[list==2])
g3=names(list[list==3])
g4=names(list[list==4])
g5=names(list[list==5])
g6=names(list[list==6])
gl=list(Cluster1 = g1, Cluster2 = g2, Cluster3 = g3, Cluster4 = g4, Cluster5 = g5, Cluster6 = g6)
gs = gosummaries(gl)
plot(gs, fontsize = 8)
gs_exp = add_expression.gosummaries(gs, exp = data_clust,
                                    annotation = annot)

plot(gs_exp, fontsize = 8, classes = "Groups", filename = "Gosummaries_data_knn4_km.pdf")

