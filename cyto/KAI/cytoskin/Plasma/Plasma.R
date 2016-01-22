setwd("~/Documents/KAI/cytoskin/Plasma")
library("stringr")
library("limma")
data_pl=read.csv2(file="~/Documents/KAI/cytoskin/Plasma_data/LUMINEX_PLASMA.csv")


cnames_c=colnames(data_pl)
cnames_c=str_replace(cnames_c,"..pg.ml.","")
cnames_c[19:22]=str_replace(cnames_c[19:22],"\\..?{2}\\..*","")
cnames_c[1:19]=str_replace(cnames_c[1:19],"\\.","")
cnames_c=toupper(cnames_c)
colnames(data_pl)=cnames_c
colnames(data_pl)[22]="LCN2"
data_c=data_pl[1:24,]
data_p=data_pl[27:67,]
data_pl=rbind(data_c,data_p)
row.names(data_pl)=as.character(data_pl$SAMPLE)
data_pl_t=t(data_pl)
Group=colnames(data_pl_t)
Group=str_replace(Group,"C.*","C")
Group=str_replace(Group,"P.*","P")
data_pl_t=rbind(data_pl_t,Group)
data_pl=t(data_pl_t)
data_pl=data.frame(data_pl)
data_pl[,2]=as.numeric(as.character((data_pl[,2])))
rows=rownames(data_pl)
data_pl[,2:22]=sapply(data_pl[,2:22],function(x) as.numeric(as.character(x)))
rownames(data_pl)=rows

data_pl_tmp=data_pl

#data_pl_tmp[,2:22]=sapply(data_pl_tmp[,2:22],function(x) as.numeric(as.character(x)))
#data_pl_tmp[,1]=data_pl[,1]

data_pl_m=melt(data_pl_tmp)
colnames(data_pl_m)[3]="Protein"
colnames(data_pl_m)[4]="Level"
data_pl_m[,4]=log2(data_pl_m[,4])
ggplot(data_pl_m, aes(x=Level))+
  facet_wrap(~Protein)+
  theme_bw()+
  geom_density(aes(group=Group, colour=Group,fill=Group), alpha=0.5)+
  ggtitle("Concentrations of cytokins in psoriasis (P) and control (C) group in plasma")

data_pl_tmp_q=normalizeQuantiles(data_pl_tmp[,2:22])

data_pl_tmp_q=cbind(data_pl_tmp_q,data_pl_tmp[,1],data_pl_tmp[,23])

data_pl_m=melt(data_pl_tmp_q)
colnames(data_pl_m)[3]="Protein"
colnames(data_pl_m)[4]="Level"
colnames(data_pl_m)[2]="Group"
data_pl_m[,4]=log2(data_pl_m[,4])
ggplot(data_pl_m, aes(x=Level))+
  facet_wrap(~Protein)+
  theme_bw()+
  geom_density(aes(group=Group, colour=Group,fill=Group), alpha=0.5)+
  ggtitle("Concentrations of cytokins in psoriasis (P) and control (C) group in plasma after quantile normalization")

##diff expressed P vs C
library("limma")
#### PL(lesional) vs PNL(non-lesional)
data_pl_tmp_q_log=log2(data_pl_tmp_q[,1:21])
p.test=t(data_pl_tmp_q_log) #
p.test=p.test[1:21,]
p=p.test
cols=colnames(p)
rows=rownames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"C0.*","C")
attribute=str_replace(attribute,"P0.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C0.*","C")
coldata=str_replace(coldata,"P0.*","P")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
p_c=tt
genes=rownames(p_c)
p_c_heatm=p[rownames(p)%in%genes,]
colnames(p_c_heatm)=cols
p_c_heatm=p_c_heatm[ , ! apply( p_c_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(p_c_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Differentially expressed cytokins in psoriasis (P) vs control (C) group in plasma.")
Gene=rownames(tt)
tt=cbind(Gene,tt)
tt$Gene=as.character(tt$Gene)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PvsC_up_QN.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PvsC_down_QN.txt",sep="\t",row.names=F,quote=FALSE)
selected_genes=Gene
save(selected_genes,p_c,file="diffexp_genes_PvsC_QN.RData")

#diffexp without quantile normalisation
p.test=t(log2(data_pl_tmp[,2:22])) #
p.test[p.test=="-Inf"]="NA"
class(p.test)="numeric"
p=p.test
cols=colnames(p)
rows=rownames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"C0.*","C")
attribute=str_replace(attribute,"P0.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C0.*","C")
coldata=str_replace(coldata,"P0.*","P")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
p_c_log2=tt
genes=rownames(p_c_log2)
p_c_heatm=p[rownames(p)%in%genes,]
colnames(p_c_heatm)=cols
p_c_heatm=p_c_heatm[ , ! apply( p_c_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(p_c_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Differentially expressed cytokins in psoriasis (P) vs control (C) group in plasma.")
Gene=rownames(tt)
tt=cbind(Gene,tt)
tt$Gene=as.character(tt$Gene)
down_log2=tt[tt[,2]<0,]
up_log2=tt[tt[,2]>0,] 
write.table(up_log2,file="PvsC_up_log2.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down_log2,file="PvsC_down_log2.txt",sep="\t",row.names=F,quote=FALSE)
selected_genes_log2=Gene
save(selected_genes_log2,p_c,file="diffexp_genes_PvsC_log2.RData")
#with and without Quantile normalization the differentially expressed genes are the same.

###wilcoxon test
p.val=c()
data=t(data_pl_tmp_q[,1:21])
for(i in 1: nrow(data)){
  u=wilcox.test(data[i,25:65],data[i,1:24])$p.value  
  p.val=c(p.val,u)
}
p.val
p.val.adj=p.adjust(p.val, method="fdr")  
genes=cbind(rownames(data), p.val.adj)
genes=as.data.frame(genes)
genes[,1]=as.character(genes[,1])
genes[,2]=as.numeric(as.character(genes[,2]))
genes=genes[genes$p.val.adj< 0.05,]
wlx_genes=genes
#######################



##########outliers detection
data_pl_ol=data_pl_tmp_q_log
data_pl_ol[data_pl_ol=="-Inf"]="NA"
#filter out controls
data_pl_ol=cbind(Sample=rownames(data_pl_ol),data_pl_ol)
data_pl_ol$Sample=as.character(data_pl_ol$Sample)


#for controls
data_pl_ol_c=data_pl_tmp_q_log[1:24,]
dim(data_pl_tmp_q_log)
data_pl_ol_c=cbind(Sample=rownames(data_pl_ol_c),data_pl_ol_c)
data_f_c=data.frame(rownames(data_pl_ol_c))
rownames(data_f_c)=rownames(data_pl_ol_c)
colnames(data_f_c)="Sample"
for(i in 2:length(colnames(data_pl_ol_c))){
  #head(data_pl_ol_c[,c(1,i)])
  
  data_pl_ol_f=data_pl_ol_c[data_pl_ol_c[,i] > quantile(data_pl_ol_c[,i], .25, na.rm = T) - 1.5*IQR(data_pl_ol_c[,i],na.rm = T) & data_pl_ol_c[,i] < quantile(data_pl_ol_c[,i], .75, na.rm = T) + 1.5*IQR(data_pl_ol_c[,i],na.rm = T),c(1,i)]
  print(head(data_pl_ol_f))
  data_f_c=merge(data_f_c,data_pl_ol_f,by.x="Sample",by.y="Sample",all.x=TRUE)
  #print(head(data_f_c))
}
dim(data_f_c)

##for psoriasis
data_pl_ol_p=data_pl_tmp_q_log[25:65,]
data_pl_ol_p=cbind(Sample=rownames(data_pl_ol_p),data_pl_ol_p)
data_f_p=data.frame(rownames(data_pl_ol_p))
rownames(data_f_p)=rownames(data_pl_ol_p)
colnames(data_f_p)="Sample"
for(i in 2:length(colnames(data_pl_ol_p))){
  
  data_pl_ol_f=data_pl_ol_p[data_pl_ol_p[,i] > quantile(data_pl_ol_p[,i], .25, na.rm = T) - 1.5*IQR(data_pl_ol_p[,i],na.rm = T) & data_pl_ol_p[,i] < quantile(data_pl_ol_p[,i], .75, na.rm = T) + 1.5*IQR(data_pl_ol_p[,i],na.rm = T),c(1,i)]
  print(head(data_pl_ol_f))
  data_f_p=merge(data_f_p,data_pl_ol_f,by.x="Sample",by.y="Sample",all.x=TRUE)
  #print(head(data_f_c))
}
dim(data_f_p)
###plot distributions after filtering out outliers
data_f=rbind(data_f_c,data_f_p)
data_f=cbind(Group=data_f$Sample,data_f)
data_f$Group=str_replace(data_f$Group,"C.*","C")
data_f$Group=str_replace(data_f$Group,"P.*","P")

data_f_m=melt(data_f)
colnames(data_f_m)[3]="Protein"
colnames(data_f_m)[4]="Level"

ggplot(data_f_m, aes(x=Level))+
  facet_wrap(~Protein)+
  theme_bw()+
  geom_density(aes(group=Group, colour=Group,fill=Group), alpha=0.5)+
  ggtitle("Concentrations of cytokins in psoriasis (P) and control (C) group in plasma after quantile normalization")


#########
data_f=rbind(data_f_c,data_f_p)
#colnames(data_f)[23]="LCN2"
rownames(data_f)=data_f$Sample
p.test=t(data_f[,2:22]) #
p=p.test
cols=colnames(p)
rows=rownames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"C0.*","C")
attribute=str_replace(attribute,"P0.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C0.*","C")
coldata=str_replace(coldata,"P0.*","P")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("P")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(p)=coldata
head(p)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
fit = lmFit(p, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
#tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt = topTable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
p_c=tt
genes=rownames(p_c)
p_c_heatm=p[rownames(p)%in%genes,]
colnames(p_c_heatm)=cols
p_c_heatm=p_c_heatm[ , ! apply( p_c_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(p_c_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Differentially expressed cytokins in psoriasis (P) vs control (C) group in plasma.")
Gene=rownames(tt)
tt=cbind(Gene,tt)
tt$Gene=as.character(tt$Gene)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PvsC_up_QF.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PvsC_down_QF.txt",sep="\t",row.names=F,quote=FALSE)
selected_genes_QF=Gene
save(selected_genes_QF,p_c,file="diffexp_genes_PvsC_QF.RData")
###save lasma data
data_plasma=data_f
save(data_plasma, file="data_plasma.RData")
####Plasma correlation with skin PL PNL
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
rownames(data.skin)[16]="INFG"
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

##take only PL samples
pl=data.skin[,24:56]
colnames(pl)=gsub("PL00","P",colnames(pl))
colnames(pl)=gsub("PL0","P",colnames(pl))

###take only PNL samples
pnl=data.skin[,57:87]
colnames(pnl)=gsub("PNL00","P",colnames(pnl))
colnames(pnl)=gsub("PNL0","P",colnames(pnl))

###take only psoriasis plasma samples
pp=data_plasma[25:65,]
pp=t(pp)
colnames(pp)=gsub("P00","P",colnames(pp))
colnames(pp)=gsub("P0","P",colnames(pp))

#select mutual cols and rows
#since there are no measurements for P36-P41 in pl and pnl, they are discarded.
pp=pp[c(-1),1:35]
all.colnames=unique(c(colnames(pp), colnames(pl), colnames(pnl)))#35

ps.extra=all.colnames[!all.colnames%in%colnames(pl)]#  "P9"  "P15" 
pl=cbind(pl,replicate(length(rownames(pl)),"NA"))
pl=cbind(pl,replicate(length(rownames(pl)),"NA"))
colnames(pl)[34:35]=c("P9","P15")

ps.extra=all.colnames[!all.colnames%in%colnames(pnl)]#  "P9"  "P14" "P18" "P26" 
pnl=cbind(pnl,replicate(length(rownames(pnl)),"NA"))
pnl=cbind(pnl,replicate(length(rownames(pnl)),"NA"))
pnl=cbind(pnl,replicate(length(rownames(pnl)),"NA"))
pnl=cbind(pnl,replicate(length(rownames(pnl)),"NA"))
colnames(pnl)[32:35]=c("P9","P14","P18","P26")

length(colnames(pp))#35
#about pairwise complete observations http://www.r-bloggers.com/pairwise-complete-correlation-considered-dangerous/
to.upper<-function(X) t(X)[lower.tri(X,diag=FALSE)]
genes_pp=rownames(pp)
rownames(pl)=gsub("IL8b","CXCL8", rownames(pl))
rownames(pnl)=gsub("IL8b","CXCL8", rownames(pnl))
genes_pl_pnl=rownames(pl)
genes=genes_pl_pnl[genes_pl_pnl%in%genes_pp]#"CXCL10" "INFG"   "IL10"   "IL17A"  "IL17F"  "IL22"   "LCN2" 
correlations=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=t(as.data.frame(correlations))
colnames(correlations)=c("lower", "r", "upper", "p", "adj.pval", "gene")
correlations=correlations[-1,]
for (i in 1:length(genes)){
  #g="IL17A"
  g=genes[i]
  print(g)
  DF=rbind(pl[rownames(pl)%in%g,],
           pnl[rownames(pnl)%in%g,],
           pp[rownames(pp)%in%g, ])
  DF=data.frame(DF)
  DF=sapply(DF,function(x) as.numeric(as.character(x)))
  rownames(DF)=c("PL","PNL", "PP")
  #cor(t(DF),method="spearman",use="pairwise.complete.obs")
  #print(DF)
  CR=corr.test(t(DF), method="spearman",,adjust ="fdr") 
  #CR=corr.test(t(DF), method="pearson",adjust ="fdr")
  adj.pval=to.upper(CR$p)
  cor.gene=cbind(CR$ci, adj.pval)
  cor.gene=cbind(cor.gene,gene=replicate(length(rownames(cor.gene)),g))
  cor.gene=cbind(groups=rownames(cor.gene), cor.gene)
  correlations=rbind(correlations, cor.gene)
 }

corpl=correlations
groups=as.character(correlations$groups)
cc <- strsplit(groups,'-')
part2 <- unlist(cc)[2*(1:length(groups))]
corpl=cbind(part2, correlations)
correlations_pl_pp=corpl[corpl$part2%in%"PP",c(-1)]
correlations_pl_pp=correlations_pl_pp[,c(1,3,5,6,7)]
save(correlations_pl_pp, file="correlations_pl_pp.RData")
write.table(file="correlations_pl_pp.txt",correlations_pl_pp, sep="\t", quote=F ) 
#correlations_pl_pp=correlations_pl_pp[complete.cases(correlations_pl_pp),]
correlations.adj.pval.sign_pl_pp= correlations_pl_pp[correlations_pl_pp$adj.pval<=0.05, ]
dim(correlations.adj.pval.sign_pl_pp)#0
save(correlations.adj.pval.sign_pl_pp, file="correlations.adj.pval.sign_pl_pp.RData")
write.table(file="correlations_adj_pval_sign_pl_pp.txt",correlations.adj.pval.sign_pl_pp, sep="\t", quote=F ) 

head(correlations_pl_pp)
conc=rbind(pl[rownames(pl)%in%genes,],pnl[rownames(pnl)%in%genes,],pp[rownames(pp)%in%genes,])
conc=cbind(rownames(conc), conc)
colnames(conc)[1]="Group"
conc[1:9,1]="PL"
conc[10:17,1]="PNL"
conc[18:14,1]="PP"
