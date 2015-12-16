#Metaanalysis psoriasis
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")
load("~/Documents/KAI/cytoskin/after_filt/21112014/selected_genes_wlx.RData")

selected_genes_wlx
data_log=data_log2_filt[,1:87]
colnames(data_log)=str_replace(colnames(data_log),"PH","PS")
colnames(data_log)=str_replace(colnames(data_log),"PT","PH")
pheatmap(data_log, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controls and Psoriasis")
summary(data_log)
data=data_log
data=as.data.frame(data)
dt=t(data)
dt=dt[25:87,]
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[36]="sample"
dt[,36]=as.character(dt[,36])

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
rownames(annotationh)=str_replace(rownames(annotation),"P","PH")
rownames(annotationt)=str_replace(rownames(annotation),"P","PT")
annotationh=annotationh[1:35,]
annotationt=annotationt[1:35,]

load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")
data_log=data_log2_filt[,1:87]#filter out vitiligo samples
p=data_log[,25:87]
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


#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"

outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_dis_onset_pso_fam_comp=outDf
save(pso_dis_onset_pso_fam_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Disease_onset_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Disease onset:late (> 40 years), early (< 40 years)")

#ggsave('/Users/nikolaeva/testplot.png', height = 11, width = 8.5)




############################################################################################
################################################# Psoriasis in the family ##################

#annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Psoriasis_in_family")]
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Gender","Psoriasis_in_family")]

annotation.h[,2]=str_replace(annotation.h[,2],"au.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"si.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"fa.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"mo.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"so.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"bro.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"chi.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"2_.*","yes")
annotation.h=annotation.h[order(annotation.h[,2]),]

#annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Psoriasis_in_family")]

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Gender","Psoriasis_in_family")]
annotation.t[,2]=str_replace(annotation.t[,2],"au.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"si.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"fa.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"mo.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"so.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"bro.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"chi.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"2_.*","yes")
annotation.t=annotation.t[order(annotation.t[,2]),]
#rownames(annotationh)=str_replace(rownames(annotation),"P","PH")
#rownames(annotationt)=str_replace(rownames(annotation),"P","PT")
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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_fam_groups_comp=outDf
save(pso_fam_groups_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_in_family_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Psoriasis in the family: yes/no")

####################################### ##Psoriatic artritis##############################
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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_arthritis_groups_comp=outDf


save(pso_arthritis_groups_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriatic_arthritis_in_family_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy);Psoriatic arthritis: yes/no")



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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_nail_inv_groups_comp=outDf


save(pso_nail_inv_groups_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasisc_nail_involvement_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Nail involvement: yes/no")


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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_pasi_groups_comp=outDf


save(pso_pasi_groups_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_PASI_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy);PASI score: low (0,10], medium (10,20], high (20,Inf]")

################################# Age >40, age<40 and psoriasis  in the family ########################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"old"
annotation.h[annotation.h[,1]< 40,1]<-"young"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"old"
annotation.t[annotation.t[,1]<40, 1]<-"young"


annotation.h[,2]=str_replace(annotation.h[,2],"au.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"si.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"fa.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"mo.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"so.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"bro.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"chi.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"2_.*","yes")
annotation.h=annotation.h[order(annotation.h[,2]),]


annotation.t[,2]=str_replace(annotation.t[,2],"au.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"si.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"fa.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"mo.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"so.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"bro.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"chi.*","yes")
annotation.t[,2]=str_replace(annotation.t[,2],"2_.*","yes")
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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_dis_onset_pso_fam_comp=outDf


save(pso_dis_onset_pso_fam_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_dis_onset_pso_fam_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Disease onset:old(>40), young(<40), Psoriasis in family: yes/no ")

#ggsave('/Users/nikolaeva/testplot.png', height = 11, width = 8.5)


################################### Age of disease onset <40 yr, >= 40 yr, PASI ######################

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("PASI_(activity_score)","Age_of_disease_onset")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.h[annotation.h[,2]>=40, 2]<-"old"
annotation.h[annotation.h[,2]< 40,2]<-"young"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "PASI_(activity_score)")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.t[annotation.t[,2]>=40, 2]<-"old"
annotation.t[annotation.t[,2]<40, 2]<-"young"

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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
pso_dis_onset_pasi_comp=outDf


save(pso_dis_onset_pasi_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_dis_onset_pasi_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
 # theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Disease onset:old(>40), young(<40),PASI score: (0,10]/(10,20]/(20,Inf]")
ggsave('~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/pso_phpt_dis_onset_pasi.png', height = 20, width = 11)

#expression of some groups is missing therefore we have gaps on the plots. See the example of AIM2 gene below.
#> tmp_aim2=tmp2[tmp2$variable%in%c("AIM2"),]
#> tmp_aim2[tmp_aim2$interactions%in%c("PH.old.20_Inf"),]
#interactions variable expression
#1  PH.old.20_Inf     AIM2         NA
#14 PH.old.20_Inf     AIM2         NA
#32 PH.old.20_Inf     AIM2         NA

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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
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
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_skin_phototype_comp=outDf


save(pso_skin_phototype_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_pso_skin_phototype_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  # theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Skin phototype: I/II/III")

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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)


tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Duration_of_psoriasis <- unname(tmp2$Duration_of_psoriasis)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_durations_comp=outDf


save(pso_durations_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_durations_comp_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  # theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Duration_of_psoriasis: <5,5-15,>15 years")


################################# Age >40, age<40 and nail involvment ########################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"old"
annotation.h[annotation.h[,1]< 40,1]<-"young"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"old"
annotation.t[annotation.t[,1]<40, 1]<-"young"

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
tmp <- subset(dt.merged.melt,variable%in%selected_genes_wlx)

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Nail_involvment <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PH","PS")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
outDf[,6]=str_replace(outDf[,6],"PT","PH")
pso_dis_onset_nail_inv_comp=outDf


save(pso_dis_onset_nail_inv_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/Psoriasis_dis_onset_nail_inv_comp_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  # theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PS (psoriasis sick),PH (psoriasis healthy); Nail involvement: yes/no")

ggsave('~/Documents/KAI/cytoskin/after_filt/21112014/anova_pheno_ver1/pso_dis_onset_nail_inv_log.png', height = 20, width = 11)
