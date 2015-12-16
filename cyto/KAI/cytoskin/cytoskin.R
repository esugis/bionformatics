setwd(setwd("~/Documents/KAI/cytoskin"))
library(pheatmap)

#read in the data c-for controlls, p-psoriasis, v-vitiliigo
c=read.table(file="controlls_all.txt",sep = "\t", header=T)
head(c)
str(c)

p=read.delim2(file="psoriasis_all.txt", sep="\t", header=F, dec = ",")
head(p)
str(p)

v=read.table(file="vitiliigo_all.txt", sep="\t", header=F)
head(v)
str(v)

#split data from c,p,v to:
#        metadata
#        skin
#        blood
#### CONTROLLS #########
rownames(c)=c[,1]
###########################metadata controlls
c.meta=c[,1:23]
head(c.meta)
###########################skin controlls
c.skin=c[,24:64]
head(c.skin)
c.skin=t(c.skin)
genes=rownames(c.skin)
c.skin=sub(",",".",(c.skin))
c.skin=as.matrix(c.skin)
c.skin <- apply(c.skin,2,as.numeric)
rownames(c.skin)=genes
pheatmap(c.skin, cluster_rows = F, cluster_cols = F, scale = "row", main="Controlls.skin")

############################blood controlls
c.blood=c[,65:106]
head(c.blood)
####separate SEB, PMA stimulated if needed
c.blood.pma1=t(c.blood[c(T,F,F)]) #PMA1part
cpma1=colnames(c.blood.pma1)
cpma1=sub("C0","PMA1.",cpma1)
colnames(c.blood.pma1)=cpma1

genes=rownames(c.blood.pma1)
genes=sub(".PMA_I.","",genes)
rownames(c.blood.pma1)=genes

c.blood.seb=t(c.blood[c(F,T,F)])#SEB part
cseb=colnames(c.blood.seb)
cseb=sub("C0","SEB.",cseb)
colnames(c.blood.seb)=cseb
rownames(c.blood.seb)=genes

c.blood.q72h=t(c.blood[c(F,F,T)])#Q72h part
cq72h=colnames(c.blood.q72h)
cq72h=sub("C0","Q72h.",cq72h)
colnames(c.blood.q72h)=cq72h
rownames(c.blood.q72h)=genes

c.blood.sort=cbind(c.blood.pma1, c.blood.seb,c.blood.q72h)#merged, transformed, genes as rows
c.blood.sort=sub(",",".",(c.blood.sort))
c.blood.sort=as.matrix(c.blood.sort)
c.blood.sort <- apply(c.blood.sort,2,as.numeric)
rownames(c.blood.sort)=genes
pheatmap(c.blood.sort, cluster_rows = F, cluster_cols = F, scale = "row", main="Controlls.Blood.")

####

##### PSORIASIS ##########
##########################metadata pso
p=p[1:37,]
p.meta=p[,1:37]
head(p.meta)
colnames(p.meta)=as.character(unname(unlist(p.meta[2,])))
p.meta=p.meta[3:37,]
rownames(p.meta)=as.character(unname(unlist(p.meta[,1])))
##########################blod pso
p.blood=p[,120:161]
head(p.blood)
colnames(p.blood)=as.character(unname(unlist(p.blood[1,])))
p.blood=p.blood[3:37,]
summary(p.blood)
str(p.blood)
head(p.blood)
rownames(p.blood)=as.character(unname(unlist(p.meta[,1])))

p.blood.pma1=t(p.blood[c(T,F,F)]) #PMA1part
pma1=colnames(p.blood.pma1)
pma1=sub("P0","PMA1.",pma1)
colnames(p.blood.pma1)=pma1
genes=rownames(p.blood.pma1)
genes=sub("?.PMA_I)","",genes)
rownames(p.blood.pma1)=genes

p.blood.seb=t(p.blood[c(F,T,F)])#SEB part
seb=colnames(p.blood.seb)
seb=sub("P0","SEB.",seb)
colnames(p.blood.seb)=seb
rownames(p.blood.seb)=genes

p.blood.q72h=t(p.blood[c(F,F,T)])#Q72h part
q72h=colnames(p.blood.q72h)
q72h=sub("P0","Q72h.",q72h)
colnames(p.blood.q72h)=q72h
rownames(p.blood.q72h)=genes

p.blood.sort=cbind(p.blood.pma1, p.blood.seb,p.blood.q72h)#merged, transformed, genes as rows
p.blood.sort=sub(",",".",(p.blood.sort))
p.blood.sort=as.matrix(p.blood.sort)
p.blood.sort <- apply(p.blood.sort,2,as.numeric)
rownames(p.blood.sort)=genes
pheatmap(p.blood.sort, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriaas.Blood")


##########################end blood pso


#########################skin pso
p.skin=p.skin[3:37,]
p.skin=p[,38:119]
head(p.skin)
#take every other column as
p.skin.h=p.skin[c(T,F)]#PH haige
p.skin.t=p.skin[c(F,T)]#PT terve
#name columns as genes
colnames(p.skin.h)=as.character(unname(unlist(p.skin.h[1,])))
colnames(p.skin.t)=as.character(unname(unlist(p.skin.h[1,])))
p.skin.h=p.skin.h[3:37,]
p.skin.t=p.skin.t[3:37,]
#name rows as samples(patients)
rownames(p.skin.h)=as.character(unname(unlist(p.meta[,1])))
rownames(p.skin.t)=as.character(unname(unlist(p.meta[,1])))
#merge PH and PT
p.skin.h=t(p.skin.h)
haige=colnames(p.skin.h)
haige=sub("","H",haige)
colnames(p.skin.h)=haige

p.skin.t=t(p.skin.t)
terve=colnames(p.skin.t)
terve=sub("P","T",terve)
colnames(p.skin.t)=terve

p.skin.ht=cbind(p.skin.h, p.skin.t)
str(p.skin.ht)
p.skin.ht=sub(",",".",(p.skin.ht))
genes=rownames(p.skin.ht)
p.skin.ht <- apply(p.skin.ht,2,as.numeric)
rownames(p.skin.ht)=genes
library(functional)
p.skin.ht.na=p.skin.ht[apply(p.skin.ht, 1, Compose(is.finite, all)),]##remove all rows with NAs in matrix
#p.skin.ht.na=as.data.frame(p.skin.ht)
#p.skin.ht.na=p.skin.ht.na[complete.cases(p.skin.ht.na),]#remove all rows with NAs in DF

#Visualise haige&terve nahk

pheatmap(p.skin.ht, cluster_rows = F, cluster_cols = F, scale = "row", main="Psoriaas.Skin")# REG3A behaves weird in the raw data
df=as.data.frame(p.skin.ht)

df=df[ , ! apply( df , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
df <- apply(df,2,as.numeric)
pheatmap(df, cluster_rows = F, cluster_cols = F, scale = "row")

#df=df[ , ! apply( df , 1 , function(x) all(is.na(x)) ) ] removes rows with full NAs
#df <- apply(df,2,as.numeric)
#p.skin.ht=t(p.skin.ht)
#pheatmap(p.skin.ht, cluster_rows = F, cluster_cols = F, scale = "row")
summary(df)
#Add normalization? Quantile? or just standardization as for pheatmap (x-mean)/std
#p.skin.ht_qntl=normalize.quantiles(p.skin.ht)
#pheatmap(p.skin.ht_qntl, cluster_rows = F, cluster_cols = F,scale = "row")
############################ end skin pso







