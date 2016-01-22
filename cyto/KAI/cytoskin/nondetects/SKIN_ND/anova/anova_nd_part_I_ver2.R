library("limma")
library(stringr)
library(pheatmap)


load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/anova/Disease_onset_anova.RData")
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
     data=data_skin_cp_nd
     data=data[,c(1:71, 73:87)]
     datat=t(data)
     #### PL vs PNL
     p.test=data[,24:86] #
     p=p.test
     cols=colnames(p)
     rows=rownames(p)
     colnames(p)=cols
     rownames(p)=rows
     attribute=colnames(p)
     attribute=str_replace(attribute,"0.*","")
     coldata=colnames(p)
     coldata=str_replace(coldata,"0.*","")
     attribute = factor(attribute)
     attribute2 = rep("Muu", length(attribute))
     # Case-control variables and corresponding group values 
     case.vals=c("PL")
     cont.vals =c("PNL")
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
     ph_pt=tt
     genes=rownames(ph_pt)
     ph_pt_heatm=p[rownames(p)%in%genes,]
     colnames(ph_pt_heatm)=cols
     ph_pt_heatm=ph_pt_heatm[ , ! apply( ph_pt_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
     pheatmap(ph_pt_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in psoriasis lesional skin (PL) vs nonlesional skin (PNL) after imputation")
     pl_pnl=rownames(tt)
     save(pl_pnl,file="pl_vs_pnl.RData")
     Gene=rownames(tt)
     tt=cbind(Gene,tt)
     down=tt[tt[,2]<0,]
     up=tt[tt[,2]>0,] 
     #write.table(up,file="PLvsPNLup_filt.txt",sep="\t",row.names=F,quote=FALSE)
     #write.table(down,file="PLvsPNLdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
     #PhvsPt=Gene#NEW
     # [1] "IL1F6"   "S100A9"  "S100A8"  "PI3"     "IL17A"   "FOXP3"   "IL17F"   "CCL27"   "CXCL1"   "LCN2"   
     #[11] "IL22"    "IFNGR"   "IL8b"    "CXCL10"  "IFNG"    "CCL2"    "CTLA4"   "IL1b"    "IL1RN"   "CCL20"  
     #[21] "CXCL2"   "IL22RA1" "AIM2"    "IL26"    "PYCARD"  "TRGC1"   "IL10"    "NLRP1"   "CCL5"    "TNF"    
     #[31] "MICB"  
     #OLD
     #[1] "IL17A"   "PI3"     "S100A8"  "S100A9"  "CCL27"   "FOXP3"   "IL1F6"   "IL1b"    "CXCL2"   "IL8b"    "CCL20"  
     #[12] "AIM2"    "CCL2"    "CXCL10"  "IFNG"    "PYCARD"  "IL1RN"   "IL22RA1" "TNF"     "IFNGR"   "IL22RA2" "IFIH1"  
     #[23] "NLRP1"
     
     
     
     ########### PL vs C
     t=data[,1:55]
     cols=colnames(t)
     attribute=colnames(t)
     attribute=str_replace(attribute,"0.*","")
     coldata=colnames(t)
     coldata=str_replace(coldata,"0.*","")
     attribute = factor(attribute)
     attribute2 = rep("Muu", length(attribute))
     # Case-control variables and corresponding group values 
     case.vals=c("PL")
     cont.vals =c("C")
     attribute2[attribute %in% c(case.vals)] = "Group1"
     attribute2[attribute %in% c(cont.vals)] = "Group2"
     attribute = factor(attribute2)
     colnames(t)=coldata
     head(t)
     mm = model.matrix(~ attribute - 1)
     colnames(mm) = make.names(colnames(mm))
     colnames(mm)
     contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
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
     pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in lesional skin in psoriatic patients (PL) vs controls (C) after imputation")
     tt[,1]=tt[,1]*(-1)
     pl_c=rownames(tt)
     save(pl_c,file="pl_vs_c.RData")
     Gene=rownames(tt)
     tt=cbind(Gene,tt)
     down=tt[tt[,2]<0,]
     up=tt[tt[,2]>0,] 
     #write.table(up,file="PLvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
     #write.table(down,file="PlvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
     
     #NEW
     # [1] "IL1F6"   "PI3"     "IL17A"   "S100A9"  "IL17F"   "IL22"    "IL1RN"   "S100A8"  "FOXP3"   "CXCL1"  
     #[11] "CCL20"   "CXCL10"  "IFNG"    "IL22RA1" "CTLA4"   "TNF"     "LCN2"    "CCL2"    "IL10"    "AIM2"   
     #[21] "MICB"    "CCL5"    "CXCL2"   "CCL27"   "EOMES"   "IL8b"    "IL26"    "IL22RA2" "TRGC1"   "IL1b"   
     #[31] "PYCARD"  "IL20RA" 
     #OLD
     #[1] "IL10"    "IL17F"   "IL17A"   "PI3"     "S100A9"  "S100A8"  "IL1F6"   "LCN2"    "CCL20"   "CXCL10"  "IFNG"   
     #[12] "CXCL1"   "CCL27"   "CXCL8"   "CTLA4"   "IL1RN"   "CXCL2"   "CCL2"    "IL1b"    "FOXP3"   "AIM2"    "IL22RA1"
     #[23] "TNF"     "DEFB1"   "CCL5"    "IFIH1"   "IL26"    "EOMES"  
     
     ############# PNL vs C
     
     t=data[, c(1:23,56:86)]
     cols=colnames(t)
     attribute=colnames(t)
     attribute=str_replace(attribute,"0.*","")
     coldata=colnames(t)
     coldata=str_replace(coldata,"0.*","")
     attribute = factor(attribute)
     attribute2 = rep("Muu", length(attribute))
     # Case-control variables and corresponding group values 
     case.vals=c("PNL")
     cont.vals =c("C")
     attribute2[attribute %in% c(case.vals)] = "Group1"
     attribute2[attribute %in% c(cont.vals)] = "Group2"
     attribute = factor(attribute2)
     colnames(t)=coldata
     head(t)
     mm = model.matrix(~ attribute - 1)
     colnames(mm) = make.names(colnames(mm))
     colnames(mm)
     contrast = makeContrasts(attributeGroup1 - attributeGroup2, levels = colnames(mm))
     fit = lmFit(t, mm)
     fit = contrasts.fit(fit, contrast)
     fit = eBayes(fit)
     adj.method="fdr" 
     thr.p.value=0.05
     thr.lfc=1
     tt = topTable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
     #if(nrow(tt) > 0){
     #  tt = tt[, c("ID", "ID", "logFC", "P.Value", "adj.P.Val")]
     #}
     
     print(thr.p.value)
     print(tt) 
     pt_c=tt
     genes=rownames(pt_c)
     data_heatm=t[rownames(t)%in%genes,]
     #data_heatm=as.data.frame(t(data_heatm))
     colnames(data_heatm)=cols
     data_heatm=data_heatm[ , ! apply( data_heatm , 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
     pheatmap(data_heatm, cluster_rows = F, cluster_cols = F, scale = "row", main="Differentially expressed genes in nonlesional skin in psoriatic parients (PNL) vs controls (C) after imputation")
     tt[,1]=tt[,1]*(-1)
     pnl_c=rownames(tt)
     save(pnl_c,file="pnl_c.RData")
     Gene=rownames(tt)
     tt=cbind(Gene,tt)
     down=tt[tt[,2]<0,]
     up=tt[tt[,2]>0,] 
     #write.table(up,file="PNLvsCup_filt.txt",sep="\t",row.names=F,quote=FALSE)
     #write.table(down,file="PNLvsCdown_filt.txt",sep="\t",row.names=F,quote=FALSE)
     #PNLvsC=Gene
     #NEW
     #[1] "IL22RA1" "IL1RN"   "IL1F6"   "TNF"     "EOMES"   "CCL20"   "CXCL1"   "CTLA4"   "CXCL10"  "IFNG"
     #OLD
     #[1] "IL1F6"   "CXCL10"  "IFNG"    "IL22RA1" "CCL20"   "DEFB1" 
     ####selected genes
     selected.genes=unique(c(pl_pnl,pl_c,pnl_c))
     length(selected.genes)#34
     save(selected.genes,file="selected.genes.RData")
     
################################# Density plots of gene expression in selected genes in PL and PNL  ##############
annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t=cbind(annotation.pso.h.t, interactions=interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset))
ann=cbind(annotation.pso.h.t, sample=rownames(annotation.pso.h.t))
data.skin=data_skin_cp_nd
summary(t(data_skin_cp_nd))
data=as.data.frame(data.skin)
dt=t(data)
#dt=dt[24:87,]

dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
dt=dt[!colnames(dt)%in%"IL26"]
colnames(dt)[42]="sample"
dt[,42]=as.character(dt[,42])
dt[,42]= str_replace(dt[,42],"0.*","")
#dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt)
colnames(dt.merged.melt)[3]="expression"
#distribution in all genes
tmp <- dt.merged.melt
save(dt,file="dt.RData")
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=sample, colour=sample,fill=sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed data after imputation of non-detects.")

#distribution in selected genes
tmp <- subset(dt.merged.melt,variable%in%selected.genes)
#save(dt,file="dt.RData")
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=sample, colour=sample,fill=sample), alpha=0.5)+
  ggtitle("Density plot of gene expression of selected genes in PL, PNL and C groups")

