library("limma")
library(stringr)
library(pheatmap)
library(reshape2)
setwd("~/Documents/KAI/cytoskin/BLOOD_SKIN_COMP/")
###########################Put Controlls, Pso, Vit in one matrix##########
#load(file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data_log2_filt.RData")
#data=data_log2_filt
#pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="All data after standardization")
#controlls
DataControllsSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.control_newfilt.txt", sep="\t", header=T)
###Psofiasis#
DataPsoriasisSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.psoriasis_newfilt.txt", sep="\t", header=T )      
###Vitiliigo#
DataVitiliigoSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.vitiliigo_newfilt.txt", sep="\t", header=T )      
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

data.all=cbind(ctrls, ps,vi)
summary(data.all)

data.log=log2(data.all)
data.skin=data.log

colnames(data.skin)=str_replace(colnames(data.skin),"PH","PL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"PT","PNL")#healthy
colnames(data.skin)=str_replace(colnames(data.skin),"VH","VL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"VT","VNL")#healthy
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))
save(data.skin, file="data.skin.Rdata")

pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL, PNL, VL and VNL groups",
         filename="data_log_NAs.png")

##take only PSL samples

psl=data.skin[,25:58]
pheatmap(psl, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in psoriatic lesional skin",
         filename="PSL_log_NAs.png")




#For interleukins Il17A IL17F skin blood  draw and calculate the values for correlation profiles in PL and PMA, SEB, unstimulated

psl_il=psl[rownames(psl)%in%c("IL17A", "IL17F"),]#data for correlation comparison skin

####Load blood sample' data
#controlls
load(file="~/Documents/KAI/cytoskin/BLOOD/datacontrol_QR.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/BLOOD/datapsoriasis_QR.RData")     
###Vitiliigo#
load(file="~/Documents/KAI/cytoskin/BLOOD/datavitiligo_QR.RData")
###PCA plot of the data

ctrls=DataControllsBlood_QR
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:16]
ctrls <- apply(ctrls[,2:16],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBlood_QR)# psoriasis
rows=rownames(pso[2:16,])
cols=as.character(pso[1,])
ps <- apply(pso[2:16,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiligoBlood_QR)# vitiliigo
rows=rownames(vit[2:16,])
cols=as.character(vit[1,])
vi <- apply(vit[2:16,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[16,]=gsub("C.*_","",ctrlstst[16,])
ctrlstst=ctrlstst[,order(ctrlstst[16,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[17,]=gsub("C","",ctrlstst[17,])
ctrlstst[17,]=gsub("_.*","",ctrlstst[17,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,17]=as.numeric(as.character(ctrlstst[,17]))
C_72h=ctrlstst[ctrlstst[,16]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,17]),]
C_PMA=ctrlstst[ctrlstst[,16]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,17]),]
C_SEB=ctrlstst[ctrlstst[,16]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,17]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:15,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
ctrls.log=log2(ctrls.data)
rownames(ctrls.log)=genes
rownames(ctrls.log)
pheatmap(ctrls.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in controls (72h, PMA, SEB)")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[16,]=gsub("P.*_","",pstst[16,])
pstst=pstst[,order(pstst[16,])]
pstst=rbind(pstst, colnames(pstst))
pstst[17,]=gsub("P","",pstst[17,])
pstst[17,]=gsub("_.*","",pstst[17,])
pstst=as.data.frame(t(pstst))
pstst[,17]=as.numeric(as.character(pstst[,17]))
P_72h=pstst[pstst[,16]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,17]),]
P_PMA=pstst[pstst[,16]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,17]),]
P_SEB=pstst[pstst[,16]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,17]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:15,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))
ps.log=log2(ps.data)
rownames(ps.log)=genes
rownames(ps.log)
pheatmap(ps.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of psoriasis patients (72h, PMA, SEB)")
ps.blood.il=ps.log[rownames(ps.log)%in%c("IL17A", "IL17F"),]#data cor correlation comparson R


#plot the data in PL and blood  72h, seb, pma
skin_blood=cbind(psl[rownames(psl)%in%rownames(ps.log),], ps.log[rownames(ps.log)%in%rownames(psl),])
pheatmap(skin_blood, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in blood samples of psoriasis patients (PL,72h, PMA, SEB)")


#######Correlation of the interleukins' profiles IL17A, IL17F in PL and SEB, PMA, 72h


#about pairwise complete observations http://www.r-bloggers.com/pairwise-complete-correlation-considered-dangerous/
to.upper<-function(X) t(X)[lower.tri(X,diag=FALSE)]
genes=rownames(skin_blood)
correlations=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=t(as.data.frame(correlations))
colnames(correlations)=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=correlations[-1,]

head(psl_il)
head(ps.blood.il)
ps.blood.il.72h=ps.blood.il[,1:35]
ps.blood.il.pma=ps.blood.il[,36:70]
ps.blood.il.seb=ps.blood.il[,71:105]

colnames(ps.blood.il.seb)=gsub("_SEB","",colnames(ps.blood.il.seb))
colnames(ps.blood.il.seb)=gsub("P0","P",colnames(ps.blood.il.seb))
colnames(ps.blood.il.72h)=gsub("_72h","",colnames(ps.blood.il.72h))
colnames(ps.blood.il.72h)=gsub("P0","P",colnames(ps.blood.il.72h))
colnames(ps.blood.il.pma)=gsub("_PMA","",colnames(ps.blood.il.pma))
colnames(ps.blood.il.pma)=gsub("P0","P",colnames(ps.blood.il.pma))
colnames(psl_il)=gsub("PL00","P",colnames(psl_il))
colnames(psl_il)=gsub("PL0","P",colnames(psl_il))
length(ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17A"), ])#35
length(psl_il[rownames(psl_il)%in%c("IL17A"),])#34

#chech which columns are not present in one of the data
all.colnames=unique(c(colnames(ps.blood.il.72h), colnames(ps.blood.il.pma), colnames(ps.blood.il.seb), colnames(psl_il)))#34
length(all.colnames)#35
length(colnames(ps.blood.il.72h))#35
length(colnames(ps.blood.il.pma))#35
length(colnames(ps.blood.il.seb))#35
length(colnames(psl_il))#34
#in the dataframe where the colname is missing add colname and set it to NA. This will be used later for correlations
pls.extra=all.colnames[!all.colnames%in%colnames(psl_il)]# "P9" 
psl_il=cbind(psl_il,replicate(length(rownames(psl_il)),"NA"))
colnames(psl_il)[35]=c("P9")
g="IL17A"
DF=rbind(psl_il[rownames(psl_il)%in%g,],
         ps.blood.il.seb[rownames(ps.blood.il.seb)%in%g, ],
         ps.blood.il.pma[rownames(ps.blood.il.pma)%in%g, ],
         ps.blood.il.72h[rownames(ps.blood.il.72h)%in%g, ])
DF=data.frame(DF)
DF=sapply(DF,function(x) as.numeric(as.character(x)))
rownames(DF)=c("PL", "SEB","PMA", "72h")
#cor(t(DF),method="spearman",use="pairwise.complete.obs")
#print(DF)
CR=corr.test(t(DF), method="spearman",adjust ="fdr") 
#CR=corr.test(t(DF), method="spearman",adjust ="holm")
adj.pval=to.upper(CR$p)
cor.gene=cbind(CR$ci, adj.pval)
cor.gene=cbind(cor.gene,gene=replicate(length(rownames(cor.gene)),g))
cor.gene=cbind(groups=rownames(cor.gene), cor.gene)
correlations=rbind(correlations, cor.gene)

correlations=correlations[complete.cases(correlations),]
correlations.adj.pval.sign= correlations[correlations$adj.pval<=0.05, ] #no significant results for IL17A

g="IL17F"
DF=rbind(psl_il[rownames(psl_il)%in%g,],
         ps.blood.il.seb[rownames(ps.blood.il.seb)%in%g, ],
         ps.blood.il.pma[rownames(ps.blood.il.pma)%in%g, ],
         ps.blood.il.72h[rownames(ps.blood.il.72h)%in%g, ])
DF=data.frame(DF)
DF=sapply(DF,function(x) as.numeric(as.character(x)))
rownames(DF)=c("PL", "SEB","PMA", "72h")
#cor(t(DF),method="spearman",use="pairwise.complete.obs")
#print(DF)
CR=corr.test(t(DF), method="spearman",adjust ="fdr") 
#CR=corr.test(t(DF), method="spearman",adjust ="holm")
adj.pval=to.upper(CR$p)
cor.gene=cbind(CR$ci, adj.pval)
cor.gene=cbind(cor.gene,gene=replicate(length(rownames(cor.gene)),g))
cor.gene=cbind(groups=rownames(cor.gene), cor.gene)
correlations=rbind(correlations, cor.gene)

correlations=correlations[complete.cases(correlations),]
correlations.adj.pval.sign= correlations[correlations$adj.pval<=0.05, ] #no significant results for IL17A

