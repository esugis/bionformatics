setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND/anova")
#Load annotations for controls and patients
load(file="~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/anova/annotation.pso.h.t.RData")
head(annotation.pso.h.t)
annotation.pso.h.t$Sample=gsub("PSL","PL",annotation.pso.h.t$Sample)
annotation.pso.h.t$Sample=gsub("PSNL","PNL",annotation.pso.h.t$Sample)
rownames(annotation.pso.h.t)=gsub("PSL","PL",rownames(annotation.pso.h.t))
rownames(annotation.pso.h.t)=gsub("PSNL","PNL",rownames(annotation.pso.h.t))
head(annotation.pso.h.t)

load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

##diff expressed PL+PNL vs C
library("limma")
data=data.skin
datat=t(data)
#### PL(lesional) vs PNL(non-lesional)
p.test=data.skin #
p=p.test
cols=colnames(p)
rows=rownames(p)
colnames(p)=cols
rownames(p)=rows
attribute=colnames(p)
attribute=str_replace(attribute,"C0.*","C")
attribute=str_replace(attribute,"PL0.*","P")
attribute=str_replace(attribute,"PNL0.*","P")
coldata=colnames(p)
coldata=str_replace(coldata,"C0.*","C")
coldata=str_replace(coldata,"PL0.*","P")
coldata=str_replace(coldata,"PNL0.*","P")
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
pheatmap(p_c_heatm, cluster_rows = F, cluster_cols = F, scale = "none", main="Differentially expressed genes in psoriasis lesional skin (PL) and nonlesional skin (PNL) vs controls (C) after imputation")
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="PL_PNLvsC_up.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="PL_PNLvsC_down_filt.txt",sep="\t",row.names=F,quote=FALSE)
selected_genes=Gene
save(selected_genes,file="diffexp_genes_PL_PLN_vs_C.RData")

################################# Density plots of gene expression in selected genes in PL and PNL  ##############
annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t=cbind(annotation.pso.h.t, interactions=interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset))
ann=cbind(annotation.pso.h.t, sample=rownames(annotation.pso.h.t))
data=as.data.frame(data.skin)
dt=t(data)
dt=dt[24:87,]
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[43]="sample"
dt[,43]=as.character(dt[,43])
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[10]="expression"

tmp <- subset(dt.merged.melt,variable%in%selected_genes)
save(dt,file="dt.RData")
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed data after imputation of non-detects.")
