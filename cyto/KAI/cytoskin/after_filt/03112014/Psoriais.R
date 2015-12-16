#This script contains analysis of  Psoriasis vs Controlls
#Read filtered scaled full data set with NAs
data_all=read.table(file="~/Documents/KAI/cytoskin/after_filt/31102014/data_all_scaled_with_NA.txt", header=T)
library(pheatmap)
rows=rownames(data_all)
cols=colnames(data_all)

data_all=apply(data_all,2,as.numeric)
colnames(data_all)=cols
rownames(data_all)=rows
pheatmap(data_all, cluster_rows = F, cluster_cols = F,scale = "row", main="All data scaled, contains NA")

#take psoriasis H, t, and controlls from the data set.


ctrls=data_all[,1:25]
ph=data_all[,26:59]
pt=data_all[,60:94]

#PHvsC
data=cbind(ctrls,ph)
t=data[,1:59]
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PH")
cont.vals =c("C")
attribute2[attribute %in% c(case.vals)] = "Group1"
attribute2[attribute %in% c(cont.vals)] = "Group2"
attribute = factor(attribute2)
colnames(t)=coldata
head(t)
mm = model.matrix(~ attribute - 1)
colnames(mm) = make.names(colnames(mm))
colnames(mm)
contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
fit = lmFit(t, mm)
fit = contrasts.fit(fit, contrast)
fit = eBayes(fit)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
if(nrow(tt) > 0){
  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
}

print(thr.p.value)
print(tt) 
ph_c=tt
genes=rownames(ph_c)
data_heatm=t[rownames(t)%in%genes,]
colnames(data_heatm)=cols
data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pheatmap(data_heatm, cluster_rows = T, cluster_cols = F, scale = "row", main="Psoriasis. Differentially expressed in PHvsC")
tt[,1]=tt[,1]*(-1)
Gene=rownames(tt)
tt=cbind(Gene,tt)
down=tt[tt[,2]<0,]
up=tt[tt[,2]>0,] 
write.table(up,file="~/Documents/KAI/cytoskin/after_filt/31102014/PHvsCup.txt",sep="\t",row.names=F,quote=FALSE)
write.table(down,file="~/Documents/KAI/cytoskin/after_filt/31102014/PHvsCdown.txt",sep="\t",row.names=F,quote=FALSE)