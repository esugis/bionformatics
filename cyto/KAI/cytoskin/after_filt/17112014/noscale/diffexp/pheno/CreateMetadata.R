#create metadata
library(pheatmap)

#read metadata for controlls
c=read.table(file="~/Documents/KAI/cytoskin/controlls_all.txt",sep = "\t", header=T)
c.meta=c[,1:23]
head(c.meta)
rownames(c.meta)=as.character(unname(unlist(c.meta[,1])))
colnames(c.meta)[23]=c("Psoriasis_in_family")
c.meta=t(c.meta[2:23])
rownames(c.meta)[3:4]=c("Height","Weight")
#read metadata for psoriasis
p=read.delim2(file="~/Documents/KAI/cytoskin/psoriasis_all.txt", sep="\t", header=F, dec = ",")
p=p[1:37,]
p.meta=p[,1:37]
head(p.meta)
colnames(p.meta)=as.character(unname(unlist(p.meta[2,])))
p.meta=p.meta[3:37,]
rownames(p.meta)=as.character(unname(unlist(p.meta[,1])))
p.meta=t(p.meta[,2:37])

#read metadata for Vitiligo
v=read.table(file="~/Documents/KAI/cytoskin/vitiliigo_all.txt", sep="\t", header=T)
v=v[1:68,]
v.meta=v[,1:68]
head(v.meta)
colnames(v.meta)=as.character(unname(unlist(v.meta[1,])))
v.meta=v.meta[2:19,]
rownames(v.meta)=as.character(unname(unlist(v.meta[,1])))
v.meta=t(v.meta[,2:68])
rownames(v.meta)[33]=c("Kobner_fenomen")

####merge metadata
c.p.meta=merge(c.meta,p.meta,by="row.names", all=T)
c.p.meta[c.p.meta=="ei"]<-"no"
rownames(c.p.meta)=c.p.meta[,1]
c.p.meta=c.p.meta[,2:length(colnames(c.p.meta))]
meta=merge(c.p.meta,v.meta,by="row.names", all=T)
rownames(meta)=meta[,1]
meta=meta[,2:length(colnames(meta))]
write.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/pheno/metadata.txt",meta, sep="\t", quote=F )      





