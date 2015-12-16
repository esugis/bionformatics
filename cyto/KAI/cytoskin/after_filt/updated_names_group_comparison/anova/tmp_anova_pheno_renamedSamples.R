#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/updated_names_PL_PLN/anova")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
load(file = "~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/anova/data_log_pso_remaned_samples.RData")
rownames(data_log)[24]="CXCL8"
summary(data_log)
data=data_log
colnames(data)=str_replace(colnames(data),"PSL","PL")#sick
colnames(data)=str_replace(colnames(data),"PSNL","PNL")#healthy
save(data,file="data_log_pso_remaned_samples.RData")
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="Gene expression in control group and psoriatic patients")

##################Find diff expressed PL+PNL vs C
t=data
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PL","PNL")
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
tt_log = toptable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt_log_all = toptable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)
dim(tt_log)#22,5
selected_genes_log=rownames(tt_log)
selected_genes_log_all= rownames(tt_log_all)
save(tt_log,tt_log_all,selected_genes_log,selected_genes_log_all, file = "PSvsC_diffexp_log.RData")


#####################################Find  diff expressed genes in PSO vs C

###correct Z-scored data. Cahnge the sample name. remove C024, rename ILb8 to CXCL8
load("~/Documents/KAI/cytoskin/after_filt/21112014/data_all_filt_nontransformed.RData")
rownames(data_all_filt)[24]="CXCL8"
data_all_filt=data_all_filt[, c(1:22,24:118)]
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"C","C")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"PH","PL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"PT","PNL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"VH","VL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"VT","VNL")
save(data_all_filt, file = "data_zscores_filt_remaned_samples.RData")
data_raw=data_all_filt[,1:86]
data=scale_rows(data_raw)
summary(data)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controlls and Psoriasis")

#################################### Wilcoxon test Z-score data #########################################
p.val=c()

for(i in 1: nrow(data)){
  u=wilcox.test(data[i,24:86],data[i,1:23])$p.value  
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
dim(wlx_genes) #26
#V1    p.val.adj
#1     AIM2 2.058208e-04
#3     CCL2 5.376484e-05
#4    CCL20 3.004238e-05
#5    CCL27 2.301644e-03
#6     CCL5 2.058208e-04
#7    CTLA4 4.906690e-06
#8    CXCL1 1.959311e-04
#9   CXCL10 4.906690e-06
#10   CXCL2 1.342933e-03
#11   DEFB1 1.265431e-02
#12   EOMES 1.547919e-02
#13   FOXP3 5.077555e-03
#14   IFIH1 1.959311e-04
#16    IFNG 4.906690e-06
#18    IL1b 3.567789e-02
#19   IL1F6 3.131922e-08
#20   IL1RN 4.432504e-07
#22 IL22RA1 1.959311e-04
#24   CXCL8 2.656471e-04
#26    LCN2 7.279780e-06
#27    MICB 4.025851e-03
#29     PI3 3.004238e-05
#30  PYCARD 1.265431e-02
#31  S100A8 4.574261e-04
#32  S100A9 5.376484e-05
#33     TNF 1.497132e-03

selected_genes_wlx=genes[,1]
save(wlx_genes,selected_genes_wlx, file = "selected_genes_wlx.RData")
selected_genes_log_all[!selected_genes_log_all%in%selected_genes_wlx]
#[1] "IFNGR"


##########################################Continue with log transformed data 

#pheatmap(data_log, cluster_rows = F, cluster_cols = F, scale = "none", main="Gene expression in control group and psoriatic patients")
#summary(data_log)

data=as.data.frame(data)
dt=t(data)
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[36]="sample"
dt[,36]=as.character(dt[,36])
#dt[,36]=str_replace(dt[,36],"PL*","PL")# ne rabotaet
#dt[,36]=str_replace(dt[,36],"PNL*","PNL")#ne rabotaet
dt=dt[24:86,]
save(dt, file="dt.RData")

##create correct annotation
################################################ Disease onset ##################
###peredelano iz vitiliigo, Ne nashla original 
p=read.table(file="~/Documents/KAI/cytoskin/psoriasis_all.txt", sep="\t", header=T)
p.meta=p[,1:37]
head(p.meta)
colnames(p.meta)=as.character(unname(unlist(p.meta[1,])))
p.meta=p.meta[2:37,]
rownames(p.meta)=as.character(unname(unlist(p.meta[,1])))
p.meta=t(p.meta[,2:37])
annotation=t(p.meta)
annotationh=annotation
annotationt=annotation
rownames(annotationh)=str_replace(rownames(annotation),"P","PL")
rownames(annotationt)=str_replace(rownames(annotation),"P","PNL")
annotationh=annotationh[1:35,]
annotationt=annotationt[1:35,]


p=data[,24:86]#take only psoriasis samples
save(p,file="psoriasis_for_anova.RData")

dim(p)#[1] 35 63
################################################ Disease Onset ############################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Age_of_disease_onset")]
annotation.h[annotation.h[,2]>=40, 2]<-"late"
annotation.h[annotation.h[,2]< 40,2]<-"early"
annotation.h=annotation.h[order(annotation.h[,2]),]

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Age_of_disease_onset")]
annotation.t[annotation.t[,2]>=40, 2]<-"late"
annotation.t[annotation.t[,2]<40, 2]<-"early"
annotation.t=annotation.t[order(annotation.t[,2]),]

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"


#Merge annotation with the data, and from the merged data frame subset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"

##
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)

############################################ Draw data disctributions in the differentially expressed genes
####### PSL,PSNL
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed data.")


############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2 <- tmp2[,c(4,5,6)]
library(plyr)
outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_dis_onset_pso_fam_comp=outDf[outDf$Groups%in%c("PL.late-PL.early","PNL.late-PNL.early"),]
save(pso_dis_onset_pso_fam_comp, file = "Disease_onset_anova.RData")
pso_dis_onset_pso_fam_comp_signif=pso_dis_onset_pso_fam_comp[pso_dis_onset_pso_fam_comp$p_adj<0.05,]
write.table(pso_dis_onset_pso_fam_comp_signif, file="pso_dis_onset_pso_fam_comp_signif.txt", sep="\t", quote=F)
##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and disease onset: < 40 years(early)/ > 40 years(late)
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_dis_onset_pso_fam_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
  
############################################################################################
################################################# Psoriasis in the family ##################

#annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Psoriasis_in_family")]
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Psoriasis_in_family")]

annotation.h[,2]=str_replace(annotation.h[,2],"au.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"si.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"fa.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"mo.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"so.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"bro.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"chi.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"2_.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"no","sporadic")
annotation.h=annotation.h[order(annotation.h[,2]),]

#annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Psoriasis_in_family")]

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Psoriasis_in_family")]
annotation.t[,2]=str_replace(annotation.t[,2],"au.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"si.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"fa.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"mo.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"so.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"bro.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"chi.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"2_.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"no","sporadic")
annotation.t=annotation.t[order(annotation.t[,2]),]

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)
annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])

annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Psoriasis_in_family))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
####### ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_fam_groups_comp=outDf[outDf$Groups%in%c("PL.sporadic-PL.familial","PNL.sporadic-PNL.familial"),]
save(pso_fam_groups_comp, file = "Psoriasis_in_family_anova.RData")
pso_fam_groups_comp_signif=pso_fam_groups_comp[pso_fam_groups_comp$p_adj<0.05,]
write.table(pso_fam_groups_comp_signif, file="pso_fam_groups_comp_signif.txt", sep="\t", quote=F)

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and presence of psoriasis in the family: yes/no.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_fam_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
          

######################################### Psoriatic artritis##############################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Psoriatic_arthritis")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Psoriatic_arthritis")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)
annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])

annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Psoriatic_arthritis))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Psoriatic_arthritis <- unname(tmp2$Psoriatic_arthritis)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_arthritis_groups_comp=outDf[outDf$Groups%in%c("PL.yes-PL.no","PNL.yes-PNL.no"),]
pso_arthritis_groups_comp_signif=pso_arthritis_groups_comp[pso_arthritis_groups_comp$p_adj<0.05,]
write.table(pso_arthritis_groups_comp_signif, file="pso_arthritis_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_arthritis_groups_comp, file = "Psoriatic_arthritis_in_family_anova.RData")

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and presence of psoriatic arthritis: yes/no.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_arthritis_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

########################################### Nail_involvment ##################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)


annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])

annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Nail_involvment))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Nail_involvment <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_nail_inv_groups_comp=outDf[outDf$Groups%in%c("PL.yes-PL.no","PNL.yes-PNL.no"),]
pso_nail_inv_groups_comp_signif=pso_nail_inv_groups_comp[pso_nail_inv_groups_comp$p_adj<0.05,]
write.table(pso_nail_inv_groups_comp_signif, file="pso_nail_inv_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_nail_inv_groups_comp, file = "Psoriasis_nail_involvement_anova.RData")

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and cases of nail involvement: yes/no.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_nail_inv_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  xlim(-2.5, 8)+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
          

###################### PASI ########################################################################

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","PASI_(activity_score)")]
annotation.h[,2]=str_replace(annotation.h[,2],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,2]=as.numeric(as.vector(annotation.h[,2]))
annotation.h[,2]=cut(annotation.h[,2], c(0, 10, 20, Inf))

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","PASI_(activity_score)")]
annotation.t[,2]=str_replace(annotation.t[,2],",",".")
annotation.t=data.frame(annotation.t)
annotation.t[,2]=as.numeric(as.vector(annotation.t[,2]))
annotation.t[,2]=cut(annotation.t[,2], c(0, 10, 20, Inf))

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
colnames(annotation.pso.h.t)[3]="PASI"
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$PASI))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$PASI <- unname(tmp2$PASI)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_pasi_groups_comp=outDf[outDf$Groups%in%c("PL.(10,20]-PL.(0,10]","PL.(20,Inf]-PL.(10,20]", "PL.(20,Inf]-PL.(0,10]","PNL.(10,20]-PNL.(0,10]","PNL.(20,Inf]-PNL.(0,10]","PNL.(20,Inf]-PNL.(10,20]"),]
pso_pasi_groups_comp_signif=pso_pasi_groups_comp[pso_pasi_groups_comp$p_adj<0.05,]
write.table(pso_pasi_groups_comp_signif, file="pso_pasi_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_pasi_groups_comp, file = "Psoriasis_PASI_anova.RData")

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and PASI score: low (0,10], medium (10,20], high (20,Inf].
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_pasi_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
     
################################# Age >40, age<40 and psoriasis  in the family ########################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"

annotation.h[,2]=str_replace(annotation.h[,2],"au.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"si.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"fa.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"mo.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"so.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"bro.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"chi.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"2_.*","familial")
annotation.h[,2]=str_replace(annotation.h[,2],"no","sporadic")
annotation.h=annotation.h[order(annotation.h[,2]),]

annotation.t[,2]=str_replace(annotation.t[,2],"au.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"si.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"fa.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"mo.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"so.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"bro.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"chi.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"2_.*","familial")
annotation.t[,2]=str_replace(annotation.t[,2],"no","sporadic")
annotation.t=annotation.t[order(annotation.t[,2]),]

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"


#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$Psoriasis_in_family <- unname(tmp2$Psoriasis_in_family)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_dis_onset_pso_fam_comp=outDf[outDf$Groups%in%c("PL.early.sporadic-PL.early.familial","PL.late.sporadic-PL.late.familial","PNL.early.sporadic-PNL.early.familial","PNL.late.sporadic-PNL.late.familial"),]
pso_dis_onset_pso_fam_comp_signif=pso_dis_onset_pso_fam_comp[pso_dis_onset_pso_fam_comp$p_adj<0.05,]
write.table(pso_dis_onset_pso_fam_comp_signif, file="pso_dis_onset_pso_fam_comp_signif.txt", sep="\t", quote=F)

save(pso_dis_onset_pso_fam_comp, file = "Psoriasis_dis_onset_pso_fam_anova.RData")
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL),age of disease onset: <40(early), >40(late) and cases of psoriasis in family: yes/no.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_dis_onset_pso_fam_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
        
################################### Age of disease onset <40 yr, >= 40 yr, PASI ######################

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("PASI_(activity_score)","Age_of_disease_onset")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.h[annotation.h[,2]>=40, 2]<-"late"
annotation.h[annotation.h[,2]< 40,2]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "PASI_(activity_score)")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.t[annotation.t[,2]>=40, 2]<-"late"
annotation.t[annotation.t[,2]<40, 2]<-"early"

annotation.h[,1]=str_replace(annotation.h[,1],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,1]=as.numeric(as.vector(annotation.h[,1]))
annotation.h[,1]=cut(annotation.h[,1], c(0, 10, 20, Inf))
annotation.h=annotation.h[order(annotation.h[,1]),]

annotation.t[,1]=str_replace(annotation.t[,1],",",".")
annotation.t=data.frame(annotation.t)
annotation.t[,1]=as.numeric(as.vector(annotation.t[,1]))
annotation.t[,1]=cut(annotation.t[,1], c(0, 10, 20, Inf))
annotation.t=annotation.t[order(annotation.t[,1]),]

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t[,2]=str_replace (annotation.pso.h.t[,2],"10,20","10_20")
annotation.pso.h.t[,2]=str_replace (annotation.pso.h.t[,2],"0,10","0_10")
annotation.pso.h.t[,2]=str_replace (annotation.pso.h.t[,2],"20,Inf","20_Inf")
annotation.pso.h.t[,2]=str_replace (annotation.pso.h.t[,2],"\\(","")
annotation.pso.h.t[,2]=str_replace (annotation.pso.h.t[,2],"\\]","")
colnames(annotation.pso.h.t)[2]="PASI"


annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$PASI))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
colnames(tmp)[6]="Gene"
ggplot(tmp, aes(x=Gene,y=expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$PASI <- unname(tmp2$PASI)
tmp2 <- tmp2[,c(4,5,6)]
outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
#pso_dis_onset_pasi_comp=outDf[outDf$Groups%in%c("PL.late.20_Inf-PL.late.10_20","PL.late.20_Inf-PL.late.0_10","PL.late.10_20-PL.late.0_10","PL.early.20_Inf-PL.early.10_20","PL.early.20_Inf-PL.early.0_10","PL.early.10_20-PL.early.0_10",
#                                                "PNL.late.20_Inf-PNL.late.10_20","PNL.late.20_Inf-PNL.late.0_10","PNL.late.10_20-PNL.late.0_10","PNL.early.20_Inf-PNL.early.10_20","PNL.early.20_Inf-PNL.early.0_10","PNL.early.10_20-PNL.early.0_10"),]
#because of small number od samples in PL.late, PNL.late all subgroups, additional plot is produced for PL.early, PNL.early all groups comparisons.
pso_dis_onset_pasi_comp_early=outDf[outDf$Groups%in%c("PL.early.20_Inf-PL.early.10_20","PL.early.20_Inf-PL.early.0_10","PL.early.10_20-PL.early.0_10",
                                               "PNL.early.20_Inf-PNL.early.10_20","PNL.early.20_Inf-PNL.early.0_10","PNL.early.10_20-PNL.early.0_10"),]

pso_dis_onset_pasi_comp_signif_early=pso_dis_onset_pasi_comp_early[pso_dis_onset_pasi_comp_early$p_adj<0.05,]
write.table(pso_dis_onset_pasi_comp_signif_early, file="pso_dis_onset_pasi_comp_signif_early.txt", sep="\t", quote=F)

save(pso_dis_onset_pasi_comp_early, file = "Psoriasis_dis_onset_pasi_anova_early.RData")

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL),age of disease onset: <40(early),>40(late), and PASI score: low (0,10], medium (10,20], high (20,Inf].
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_dis_onset_pasi_comp_early, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

#expression of some groups is missing therefore we have gaps on the plots. See the example of AIM2 gene below.
#> tmp_aim2=tmp2[tmp2$variable%in%c("AIM2"),]
#> tmp_aim2[tmp_aim2$interactions%in%c("PH.late.20_Inf"),]
#interactions variable expression
#1  PH.late.20_Inf     AIM2         NA
#14 PH.late.20_Inf     AIM2         NA
#32 PH.late.20_Inf     AIM2         NA

#################################### Skin phototype ###########################################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Skin_phototype")]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Skin_phototype")]
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample,annotation.pso.h.t$Skin_phototype))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
############################################# ANOVA part on the selected genes in log-transformed data ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Skin_phototype <- unname(tmp2$Skin_phototype)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_skin_phototype_comp=outDf[outDf$Groups%in%c("PL.III-PL.II","PNL.III-PNL.II"),]
pso_skin_phototype_comp_signif=pso_skin_phototype_comp[pso_skin_phototype_comp$p_adj<0.05,]
write.table(pso_skin_phototype_comp_signif, file="pso_skin_phototype_comp_signif.txt", sep="\t", quote=F)

save(pso_skin_phototype_comp, file = "Psoriasis_pso_skin_phototype_anova.RData")
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and skin phototype: I/II/III.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_skin_phototype_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
  
          
###############################Duration of the  psoriasis ##################

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Duration_of_disease(years)")]
annotation.h[,2]=str_replace(annotation.h[,2],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,2]=as.numeric(as.vector(annotation.h[,2]))
annotation.h[,2]=cut(annotation.h[,2], c(0, 5, 15, Inf))

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Duration_of_disease(years)")]
annotation.t[,2]=str_replace(annotation.t[,2],",",".")
annotation.t=data.frame(annotation.t)
annotation.t[,2]=as.numeric(as.vector(annotation.t[,2]))
annotation.t[,2]=cut(annotation.t[,2], c(0, 5, 15, Inf))

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
colnames(annotation.pso.h.t)[3]="Duration_of_psoriasis"
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Duration_of_psoriasis))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)


tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Duration_of_psoriasis <- unname(tmp2$Duration_of_psoriasis)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_durations_comp=outDf[outDf$Groups%in%c("PL.(15,Inf]-PL.(0,5]","PL.(5,15]-PL.(0,5]","PL.(5,15]-PL.(0,5]","PNL.(15,Inf]-PNL.(0,5]","PNL.(5,15]-PNL.(0,5]","PNL.(5,15]-PNL.(0,5]"),]
pso_durations_comp_signif=pso_durations_comp[pso_durations_comp$p_adj<0.05,]
write.table(pso_durations_comp_signif, file="pso_durations_comp_signif.txt", sep="\t", quote=F)

save(pso_durations_comp, file = "Psoriasis_durations_comp_anova.RData")

#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and duration of psoriasis: <5,5-15,>15 years.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_durations_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

################################# Disease onset >40, <40 and nail involvment ########################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Nail_involvment))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Nail_involvment <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_dis_onset_nail_inv_comp=outDf[outDf$Groups%in%c("PL.early.yes-PL.early.no","PL.late.yes-PL.late.no","PNL.early.yes-PNL.early.no","PNL.late.yes-PNL.late.no"),]
pso_dis_onset_nail_inv_comp_signif=pso_dis_onset_nail_inv_comp[pso_dis_onset_nail_inv_comp$p_adj<0.05,]
write.table(pso_dis_onset_nail_inv_comp_signif, file="pso_dis_onset_nail_inv_comp_signif.txt", sep="\t", quote=F)

save(pso_dis_onset_nail_inv_comp, file = "Psoriasis_dis_onset_nail_inv_anova.RData")
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), age of disease onset: <40(early), >40(late), and nail involvement: yes/no".
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_dis_onset_nail_inv_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
          