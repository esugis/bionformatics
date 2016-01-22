library("limma")
library(stringr)
library(pheatmap)
library(reshape2)
library("psych")
library("ggplot2")
library("GGally")
setwd("~/Documents/KAI/cytoskin/nondetects/BLOOD_B2mg_SKIN_COMP/")
###########################Put Controlls, Pso, Vit in one matrix##########
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

##take only PSL samples
psl=data.skin[,24:56]
pheatmap(psl, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in psoriatic lesional skin",
         filename="PSL_ND.png")

#For interleukins Il17A IL17F skin blood  draw and calculate the values for correlation profiles in PL and PMA, SEB, unstimulated
#il=c("IL17A","IL17F", "IL22","IL4", "IL10" )
#psl_il=psl[rownames(psl)%in%il,]#data for correlation comparison skin

####Load blood sample' data
load("~/Documents/KAI/cytoskin/nondetects/BLOOD_ND/DDCTBloodB2mgF.RData")

ctrls=DDCTBloodB2mgF[1:70,]
rows=rownames(ctrls)
cols=colnames(ctrls)
ctrls <- apply(ctrls,2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)*(-1)#because gene expr formula is 2^(-DDCT) and we need log(gexp)

#PSO
pso=t(DDCTBloodB2mgF[71:173,])# psoriasis. because gene expr formula is 2^(-DDCT) and we need log(gexp)
rows=rownames(pso)
cols=colnames(pso)
ps <- apply(pso,2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows
ps=ps*(-1)

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[15,]=gsub("C.*_","",ctrlstst[15,])
ctrlstst=ctrlstst[,order(ctrlstst[15,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[16,]=gsub("C","",ctrlstst[16,])
ctrlstst[16,]=gsub("_.*","",ctrlstst[16,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,16]=as.numeric(as.character(ctrlstst[,16]))
C_72h=ctrlstst[ctrlstst[,15]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,16]),]
C_PMA=ctrlstst[ctrlstst[,15]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,16]),]
C_SEB=ctrlstst[ctrlstst[,15]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,16]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:14,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
rownames(ctrls.data)=genes
rownames(ctrls.data)
pheatmap(ctrls.data, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in controls (72h, PMA, SEB). Blood B2mg after imputation")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[15,]=gsub("P.*_","",pstst[15,])
pstst=pstst[,order(pstst[15,])]
pstst=rbind(pstst, colnames(pstst))
pstst[16,]=gsub("P","",pstst[16,])
pstst[16,]=gsub("_.*","",pstst[16,])
pstst=as.data.frame(t(pstst))
pstst[,16]=as.numeric(as.character(pstst[,16]))
P_72h=pstst[pstst[,15]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,16]),]
P_PMA=pstst[pstst[,15]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,16]),]
P_SEB=pstst[pstst[,15]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,16]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:14,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))

rownames(ps.data)=genes
rownames(ps.data)
pheatmap(ps.data, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in blood samples (B2mg) of psoriasis patients (72h, PMA, SEB) after imputation")
#ps.blood.il=ps.data[rownames(ps.data)%in%il,]#data cor correlation comparson R


#plot the data in PL and blood  72h, seb, pma
skin_blood=cbind(psl[rownames(psl)%in%rownames(ps.data),], ps.data[rownames(ps.data)%in%rownames(psl),])
pheatmap(skin_blood, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in skin and blood samples (B2mg) of psoriasis patients (PL,72h, PMA, SEB) after imputation")

tskin_blood=t(skin_blood)

ps.blood.72h=skin_blood[,34:68]
ps.blood.pma=skin_blood[,69:101]
ps.blood.seb=skin_blood[,102:136]
psl=skin_blood[,1:33]
colnames(ps.blood.seb)=gsub("_SEB","",colnames(ps.blood.seb))
colnames(ps.blood.seb)=gsub("P0","P",colnames(ps.blood.seb))
colnames(ps.blood.72h)=gsub("_72h","",colnames(ps.blood.72h))
colnames(ps.blood.72h)=gsub("P0","P",colnames(ps.blood.72h))
colnames(ps.blood.pma)=gsub("_PMA","",colnames(ps.blood.pma))
colnames(ps.blood.pma)=gsub("P0","P",colnames(ps.blood.pma))
colnames(psl)=gsub("PL00","P",colnames(psl))
colnames(psl)=gsub("PL0","P",colnames(psl))

#chech which columns are not present in one of the data
all.colnames=unique(c(colnames(ps.blood.72h), colnames(ps.blood.pma), colnames(ps.blood.seb), colnames(psl)))#35
length(colnames(ps.blood.72h))#35
length(colnames(ps.blood.pma))#33
length(colnames(ps.blood.seb))#35
length(colnames(psl))#33
#in the dataframe where the colname is missing add colname and set it to NA. This will be used later for correlations
pls.extra=all.colnames[!all.colnames%in%colnames(psl)]# "P9"  "P15"
psl=cbind(psl,replicate(length(rownames(psl)),"NA"))
psl=cbind(psl,replicate(length(rownames(psl)),"NA"))
colnames(psl)[34:35]=c("P9","P15")

ps.blood.pma.extra=all.colnames[!all.colnames%in%colnames(ps.blood.pma)]# "P3"  "P14"
ps.blood.pma=cbind(ps.blood.pma,replicate(length(rownames(ps.blood.pma)),"NA"))
ps.blood.pma=cbind(ps.blood.pma,replicate(length(rownames(ps.blood.pma)),"NA"))
colnames(ps.blood.pma)[34:35]=c("P3","P14")

psl=cbind(psl, group=replicate(length(rownames(psl)),"PSL"))
ps.blood.72h=cbind(ps.blood.72h, group=replicate(length(rownames(ps.blood.72h)),"72h"))
ps.blood.pma=cbind(ps.blood.pma, group=replicate(length(rownames(ps.blood.pma)),"PMA"))
ps.blood.seb=cbind(ps.blood.seb, group=replicate(length(rownames(ps.blood.seb)),"SEB"))


skin_blood_group=rbind(psl,ps.blood.72h,ps.blood.pma, ps.blood.seb)

#about pairwise complete observations http://www.r-bloggers.com/pairwise-complete-correlation-considered-dangerous/
to.upper<-function(X) t(X)[lower.tri(X,diag=FALSE)]
genes=rownames(skin_blood)
correlations=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=t(as.data.frame(correlations))
colnames(correlations)=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=correlations[-1,]
for (i in 1:length(genes)){
  #g="IL17A"
  g=genes[i]
  print(g)
  DF=rbind(psl[rownames(psl)%in%g,],
            ps.blood.seb[rownames(ps.blood.seb)%in%g, ],
            ps.blood.pma[rownames(ps.blood.pma)%in%g, ],
            ps.blood.72h[rownames(ps.blood.72h)%in%g, ])
  DF=data.frame(DF)
  DF=sapply(DF,function(x) as.numeric(as.character(x)))
  rownames(DF)=c("PL", "SEB","PMA", "72h")
  #cor(t(DF),method="spearman",use="pairwise.complete.obs")
  #print(DF)
 # CR=corr.test(t(DF), method="spearman",adjust ="fdr") 
  CR=corr.test(t(DF), method="pearson",adjust ="fdr")
  adj.pval=to.upper(CR$p)
  cor.gene=cbind(CR$ci, adj.pval)
  cor.gene=cbind(cor.gene,gene=replicate(length(rownames(cor.gene)),g))
  cor.gene=cbind(groups=rownames(cor.gene), cor.gene)
  correlations=rbind(correlations, cor.gene)
  #chart.Correlation((t(DF)),method="spearman", main="Spearman correlation between expression profiles in PL, SEB(B2mg), PMA(B2mg), 72h(B2mg) ")#check this one
}
#scale_rows = function(x){
#m = apply(x, 1, mean, na.rm = T)
#s = apply(x, 1, sd, na.rm = T)
#return((x - m) / s)
#}
corpl=correlations
groups=as.character(correlations$groups)
cc <- strsplit(groups,'-')
part1 <- unlist(cc)[2*(1:length(groups))-1]
corpl=cbind(part1, correlations)
correlations_pl=corpl[corpl$part1%in%"PL",c(-1)]
save(correlations_pl, file="correlations_skin_blood_pl.RData")
write.table(file="correlations_pl.txt",correlations_pl, sep="\t", quote=F ) 
correlations_pl=correlations_pl[complete.cases(correlations_pl),]
correlations.adj.pval.sign_pl= correlations_pl[correlations_pl$adj.pval<=0.05, ]#11,7
dim(correlations.adj.pval.sign_pl)
save(correlations.adj.pval.sign_pl, file="correlations.adj.pval.sign_pl.RData")
write.table(file="correlations_adj_pval_sign_pl.txt",correlations.adj.pval.sign_pl, sep="\t", quote=F ) 
#correlations.adj.pval.sign
#correlations.p.sign=correlations[correlations$p<=0.05, ]#14,7

#plot results only for the significant correlation results
#1 select th genes and corresponding groups from the signif cor

genes.sign=unique(as.character(correlations.adj.pval.sign_pl$gene))
#"CTLA4" 
skin_blood_group=cbind(skin_blood_group,gene=rownames(skin_blood_group))
skin_blood_group=data.frame(skin_blood_group)
skin_blood_group=skin_blood_group[skin_blood_group$gene%in%genes.sign,]

for (i in 1:length(genes.sign)){
  g=genes.sign[i]
  #g="AIM2"
  groups.sign=as.character(correlations.adj.pval.sign[correlations.adj.pval.sign$gene%in%g,1])
  groups.sign
  groups=unique(unlist(strsplit(groups.sign,"-")))
  dat=skin_blood_group[skin_blood_group$gene%in%g,]
  dat.t=(t((dat))
  class(dat.t)="numeric"    
  dat.t=data.frame(dat.t)
  colnames(dat.t)=dat$group
  colnames(dat.t)[2]="Unstimulated"
  
  ggplot(dat.t[1:35,], aes(x=PMA, y=SEB)) + geom_point(shape=1) +geom_smooth(method=lm)+
     theme_bw()+
     geom_text(data=data.frame(correlations.adj.pval.sign[1,]), aes(label=paste("r=", round(r,5), sep="")), x=1, y=-0.25)+
     geom_text(data=data.frame(correlations.adj.pval.sign[1,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=1, y=-0.5)+
     ggtitle("Correlation between expression profiles of AIM2 gene \n in groups SEB and PMA/Iono") 
     
  ggplot(dat.t[1:35,], aes(x=PMA, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
    theme_bw()+
    geom_text(data=data.frame(correlations.adj.pval.sign[2,]), aes(label=paste("r=", round(r,5), sep="")), x=1, y=2)+
    geom_text(data=data.frame(correlations.adj.pval.sign[2,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=1, y=1.75)+
    ggtitle("Correlation between expression profiles of AIM2 gene \n in groups SEB and Unstimulated") 

#"CTLA4"  
  g="CTLA4"
dat=skin_blood_group[skin_blood_group$gene%in%g,]
dat.t=t((dat))
       class(dat.t)="numeric"    
       dat.t=data.frame(dat.t)
       colnames(dat.t)=dat$group
       colnames(dat.t)[2]="Unstimulated"
       colnames(dat.t)[1]="PL"
ggplot(dat.t[1:35,], aes(x=PL, y=PMA)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[3,]), aes(label=paste("r=", round(r,5), sep="")), x=14, y=2.25)+
  geom_text(data=data.frame(correlations.adj.pval.sign[3,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=14, y=2)+
  ggtitle("Correlation between expression profiles of CTLA4 gene \n in groups PL and PMA/Iono") 

ggplot(dat.t[1:35,], aes(x=SEB, y=PMA)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[4,]), aes(label=paste("r=", round(r,5), sep="")), x=-1, y=-3)+
  geom_text(data=data.frame(correlations.adj.pval.sign[4,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-1, y=-3.25)+
  ggtitle("Correlation between expression profiles of CTLA4 gene \n in groups SEB and PMA") 

ggplot(dat.t[1:35,], aes(x=SEB, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[5,]), aes(label=paste("r=", round(r,5), sep="")), x=-1, y=2.25)+
  geom_text(data=data.frame(correlations.adj.pval.sign[5,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-1, y=2)+
  ggtitle("Correlation between expression profiles of CTLA4 gene \n in groups SEB and Unstimulated") 


#EOMES
g="EOMES"
dat=skin_blood_group[skin_blood_group$gene%in%g,]
dat.t=t((dat))
class(dat.t)="numeric"    
dat.t=data.frame(dat.t)
colnames(dat.t)=dat$group
colnames(dat.t)[2]="Unstimulated"
colnames(dat.t)[1]="PL"
ggplot(dat.t[1:35,], aes(x=SEB, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[6,]), aes(label=paste("r=", round(r,5), sep="")), x=-1.5, y=-2)+
  geom_text(data=data.frame(correlations.adj.pval.sign[6,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-1.5, y=-2.25)+
  ggtitle("Correlation between expression profiles of EOMES gene \n in groups SEB and Unstimulated") 


#FOXP3

g="FOXP3"
dat=skin_blood_group[skin_blood_group$gene%in%g,]
dat.t=t((dat))
class(dat.t)="numeric"    
dat.t=data.frame(dat.t)
colnames(dat.t)=dat$group
colnames(dat.t)[2]="Unstimulated"
colnames(dat.t)[1]="PL"
ggplot(dat.t[1:35,], aes(x=SEB, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[7,]), aes(label=paste("r=", round(r,5), sep="")), x=-2.5, y=-0.5)+
  geom_text(data=data.frame(correlations.adj.pval.sign[7,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-2.5, y=-0.75)+
  ggtitle("Correlation between expression profiles of FOXP3 gene \n in groups SEB and Unstimulated") 

#IFIH1

g="IFIH1"
dat=skin_blood_group[skin_blood_group$gene%in%g,]
dat.t=t((dat))
class(dat.t)="numeric"    
dat.t=data.frame(dat.t)
colnames(dat.t)=dat$group
colnames(dat.t)[2]="Unstimulated"
colnames(dat.t)[1]="PL"
ggplot(dat.t[1:35,], aes(x=SEB, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[9,]), aes(label=paste("r=", round(r,5), sep="")), x=-0.5, y=-0.75)+
  geom_text(data=data.frame(correlations.adj.pval.sign[9,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-0.5, y=-1)+
  ggtitle("Correlation between expression profiles of IFIH1 gene \n in groups SEB and Unstimulated") 

#SEB-PMA 8
ggplot(dat.t[1:35,], aes(x=SEB, y=PMA)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[8,]), aes(label=paste("r=", round(r,5), sep="")), x=-0.5, y=0)+
  geom_text(data=data.frame(correlations.adj.pval.sign[8,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=-0.5, y=-0.25)+
  ggtitle("Correlation between expression profiles of IFIH1 gene \n in groups SEB and PMA") 
#PMA-72h 10
ggplot(dat.t[1:35,], aes(x=PMA, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[10,]), aes(label=paste("r=", round(r,5), sep="")), x=1.75, y=-1)+
  geom_text(data=data.frame(correlations.adj.pval.sign[10,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=1.75, y=-1.25)+
  ggtitle("Correlation between expression profiles of IFIH1 gene \n in groups PMA and Unstimulated") 

#IL22
#PMA-72h
g="IL22"
dat=skin_blood_group[skin_blood_group$gene%in%g,]
dat.t=t((dat))
class(dat.t)="numeric"    
dat.t=data.frame(dat.t)
colnames(dat.t)=dat$group
colnames(dat.t)[2]="Unstimulated"
colnames(dat.t)[1]="PL"

#PMA-72h 10
ggplot(dat.t[1:35,], aes(x=PMA, y=Unstimulated)) + geom_point(shape=1) +geom_smooth(method=lm)+
  theme_bw()+
  geom_text(data=data.frame(correlations.adj.pval.sign[11,]), aes(label=paste("r=", round(r,5), sep="")), x=0, y=7)+
  geom_text(data=data.frame(correlations.adj.pval.sign[11,]), aes(label=paste("pval=", round(adj.pval,5), sep="")), x=0, y=6.75)+
  ggtitle("Correlation between expression profiles of IL22 gene \n in groups PMA and Unstimulated") 
