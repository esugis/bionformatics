#Comparison of the groups defined on the meeting 04.12

#1.early onset.family'+'
# late onset.family'-'
#means that we compare 
#PL.early.yes
#PNL.early.yes
#PL.late.no
#PNL.late.no
setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/updated_names_PL_PLN/anova")
#load the gene expression  data 
load(file="psoriasis_for_anova.RData")
#?load(file="dt.RData")
load(file="diffexp_genes_PL_PLN_vs_C.RData")
library(stringr)
library(reshape)
#load data
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))

#extract annonations
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
rownames(annotationh)=str_replace(rownames(annotation),"P","PSL")
rownames(annotationt)=str_replace(rownames(annotation),"P","PSNL")
annotationh=annotationh[1:35,]
annotationt=annotationt[1:35,]


p=data_log[,24:86]
################################# Age >40, age<40 and psoriasis  in the family ########################
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"


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

#Keep only the following groups
#PL.early.yes
#PNL.early.yes
#PL.late.no
#PNL.late.no
keep=c("PL.early.yes","PNL.early.yes", "PL.late.no","PNL.late.no")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)

#make anova + turkey test
### ANOVA part on the selected genes ##############################
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

PS_typeI_typeII_PL_PN=outDf[outDf$Groups%in%c("PL.typeII-PL.typeI","PNL.typeII-PNL.typeI"),]
#PS_EarlyFam_vs_LateNoFam$Groups=str_replace(PS_EarlyFam_vs_LateNoFam$Groups,"early.yes","typeI")
#PS_EarlyFam_vs_LateNoFam$Groups=str_replace(PS_EarlyFam_vs_LateNoFam$Groups,"early.yes","typeI")
#PS_EarlyFam_vs_LateNoFam$Groups=str_replace(PS_EarlyFam_vs_LateNoFam$Groups,"late.no","typeII")
#PS_EarlyFam_vs_LateNoFam$Groups=str_replace(PS_EarlyFam_vs_LateNoFam$Groups,"late.no","typeII")
#PS_EarlyFam_vs_LateNoFam=PS_EarlyFam_vs_LateNoFam[order(PS_EarlyFam_vs_LateNoFam$variable,PS_EarlyFam_vs_LateNoFam$p_adj),]
PS_typeI_typeII_PL_PN_signif=PS_typeI_typeII_PL_PN[PS_typeI_typeII_PL_PN$p_adj<0.05,]
save(PS_typeI_typeII_PL_PN,PS_typeI_typeII_PL_PN_signif, file = "PS_typeI_typeII_PL_PN_anova.RData")
write.table(PS_typeI_typeII_PL_PN_signif, file="PS_typeI_typeII_PL_PNL.txt", sep="\t", quote=F)
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL),and the combination of early disease onset(<40 years) with cases of psoriasis in the family; late disease onset (>40 years) with no cases of psoriasis in family.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(PS_typeI_typeII_PL_PN, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

####!!!!!!!!esli ne delit na PSL i PSNL, a vziat tolko type1(early.yes) i type2(late.no)

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"

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
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"

#Keep only the following groups
#PSL.early.yes
#PSNL.early.yes
#PSL.late.no
#PSNL.late.no
keep=c("early.yes","late.no")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)

#make anova + turkey test
### ANOVA part on the selected genes ##############################
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
outDf$Groups=str_replace(outDf$Groups,"early.yes","typeI")
outDf$Groups=str_replace(outDf$Groups,"late.no","typeII")
PS_type1_type2=outDf
PS_type1_type2=PS_type1_type2[order(PS_type1_type2$p_adj),]
save(PS_type1_type2, file = "PS_type1_type2_anova.RData")
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL),and the combination of early disease onset(<40 years) with cases of psoriasis in the family; late disease onset (>40 years) with no cases of psoriasis in family.
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

#############Early onset. High PASI (>20) Late onset. Low PASI (<10)################
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
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$PASI))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
keep=c("early.20_Inf","late.0_10")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
tmp=tmp[tmp$interactions%in%keep,]
#plot groups 
colnames(tmp)[6]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(Expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$PASI <- unname(tmp2$PASI)
tmp2 <- tmp2[,c(4,5,6)]
outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
outDf$Groups=str_replace(outDf$Groups,"20_Inf","severe")
outDf$Groups=str_replace(outDf$Groups,"0_10","mild")
ps_early.high_late.low=outDf
ps_early.high_late.low=ps_early.high_late.low[order(ps_early.high_late.low$p_adj),]
ps_early.high_late.low=outDf
ps_early.high_late.low_signif=ps_early.high_late.low[ps_early.high_late.low$p_adj<0.05,]
save(ps_early.high_late.low,ps_early.high_late.low_signif, file = "PS_early.severe_late.mild_anova.RData")
write.table(ps_early.high_late.low_signif, file="ps_early.severe_late.mild_signif.txt", quote=F, sep="\t")
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL),age of disease onset: <40(early),>40(late), and PASI score: low (0,10], medium (10,20], high (20,Inf].
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

############################### PL.early.severe PNL.early.severe;PL.late.severe  PNL.late.severe
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
keep=c("PL.early.20_Inf","PNL.early.20_Inf","PL.late.20_Inf","PNL.late.20_Inf")
ann_filt=ann[ann$interactions%in%keep,]
#ann_filt$interactions=str_replace(ann_filt$interactions,"20_Inf","severe")
#ann_filt$interactions=str_replace(ann_filt$interactions,"20_Inf","severe")
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)

colnames(dt.merged.melt)[7]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[6]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()

#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(Expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$PASI <- unname(tmp2$PASI)
tmp2 <- tmp2[,c(4,5,6)]
outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

PL_PNL_early_late_severe=outDf[outDf$Groups%in%c("PNL.late.20_Inf-PNL.early.20_Inf","PL.late.20_Inf-PL.early.20_Inf"),]
PL_PNL_early_late_severe_signif=PL_PNL_early_late_severe[PL_PNL_early_late_severe$p_adj<0.05,]
save(PL_PNL_early_late_severe,PL_PNL_early_late_severe_signif, file = "PL_PNL_early_late_severe_anova.RData")
write.table(PL_PNL_early_late_severe_signif, file="PL_PNL_early_late_severe_signif.txt", quote=F, sep="\t")

ggplot(PL_PNL_early_late_severe, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

########################## Type I + pos.PS Arthritis Type I + neg PS Arthritis; Type II + pos. PS Arthritis Type II + neg. PS Arthritis 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.h=annotation.h[,c(1,3,2)]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"
annotation.t=annotation.t[,c(1,3,2)]
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
annotation.pso.h.t[,4]=as.character(annotation.pso.h.t[,4])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family,annotation.pso.h.t$Psoriatic_arthritis))
colnames(annotation.pso.h.t)[5]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[6]="sample"
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"typeI.no","typeI.PSA-")
ann$interactions=str_replace(ann$interactions,"typeII.no","typeII.PSA-")
ann$interactions=str_replace(ann$interactions,"typeI.yes","typeI.PSA+")
ann$interactions=str_replace(ann$interactions,"typeII.yes","typeII.PSA+")
keep=c("typeI.PSA-","typeII.PSA-","typeI.PSA+","typeII.PSA+" )
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[8]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[7]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()

#make anova + turkey test
### ANOVA part on the selected genes ##############################
#Aov example for many genes
myAov <- function(df){
  aov_list <- aov(Expression~interactions, data=df)
  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
  aov_out$groups <- rownames(aov_out)
  return(aov_out)
}

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$Psoriasis_in_family <- unname(tmp2$Psoriasis_in_family)
tmp2$Psoriatic_arthritis <- unname(tmp2$Psoriatic_arthritis)
tmp2 <- tmp2[,c(5,6,7)]

outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

PS_type1_type2_psa=outDf
PS_type1_type2_psa=PS_type1_type2_psa[order(PS_type1_type2_psa$p_adj),]
PS_type1_type2_psa_signif=PS_type1_type2_psa[PS_type1_type2_psa$p_adj<0.05,]
save(PS_type1_type2_psa,PS_type1_type2_psa_signif, file = "PS_type1_type2_psa_anova.RData")
write.table(PS_type1_type2_psa_signif, file="PS_type1_type2_psa_signif.txt", quote=F, sep="\t")

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))



########################## Type I + NI(+) Type I + NI(-); Type II + NI(+) Type II + NI(-) 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.h=annotation.h[,c(1,3,2)]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"
annotation.t=annotation.t[,c(1,3,2)]
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
annotation.pso.h.t[,4]=as.character(annotation.pso.h.t[,4])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family,annotation.pso.h.t$Nail_involvment))
colnames(annotation.pso.h.t)[5]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[6]="sample"
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"typeI.no","typeI.NI-")
ann$interactions=str_replace(ann$interactions,"typeII.no","typeII.NI-")
ann$interactions=str_replace(ann$interactions,"typeI.yes","typeI.NI+")
ann$interactions=str_replace(ann$interactions,"typeII.yes","typeII.NI+")
keep=c("typeI.NI-","typeII.NI-","typeI.NI+","typeII.NI+" )
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[8]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[7]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  #scale_colour_manual(values = c("red","blue", "green","black"))+
  #geom_point(aes(shape = interactions),size = 4)+
 # geom_boxplot(aes(fill = factor(interactions)))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()



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
tmp2$Psoriatic_arthritis <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(5,6,7)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

PS_type1_type2_ni=outDf
PS_type1_type2_ni=PS_type1_type2_ni[order(PS_type1_type2_ni$p_adj),]
PS_type1_type2_ni_signif=PS_type1_type2_ni[PS_type1_type2_ni$p_adj<0.05,]
save(PS_type1_type2_ni,PS_type1_type2_ni_signif, file = "PS_type1_type2_ni_anova.RData")
write.table(PS_type1_type2_ni_signif, file="PS_type1_type2_ni_signif.txt", quote=F, sep="\t")

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))


#######  Type I + PSA(+) + PASI moderate (10-20), Type I + PSA(+) + PASI mild ( <10), Type II + PSA(+) + PASI moderate (10-20), Type II + PSA(+) + PASI mild ( <10) 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis","PASI_(activity_score)")]
annotation.h=annotation.h[,c(2,4,3,1)]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis","PASI_(activity_score)")]
annotation.t=annotation.t[,c(2,4,3,1)]
annotation.t=annotation.t[order(annotation.t[,1]),]
annotation.t[annotation.t[,1]>=40, 1]<-"late"
annotation.t[annotation.t[,1]<40, 1]<-"early"


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


annotation.h[,4]=str_replace(annotation.h[,4],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,4]=as.numeric(as.vector(annotation.h[,4]))
annotation.h[,4]=cut(annotation.h[,4], c(0, 10, 20, Inf))
annotation.h=annotation.h[order(annotation.h[,4]),]

annotation.t[,4]=str_replace(annotation.t[,4],",",".")
annotation.t=data.frame(annotation.t)
annotation.t[,4]=as.numeric(as.vector(annotation.t[,4]))
annotation.t[,4]=cut(annotation.t[,4], c(0, 10, 20, Inf))
annotation.t=annotation.t[order(annotation.t[,4]),]

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t[,5]=str_replace (annotation.pso.h.t[,5],"10,20","10_20")
annotation.pso.h.t[,5]=str_replace (annotation.pso.h.t[,5],"0,10","0_10")
annotation.pso.h.t[,5]=str_replace (annotation.pso.h.t[,5],"20,Inf","20_Inf")
annotation.pso.h.t[,5]=str_replace (annotation.pso.h.t[,5],"\\(","")
annotation.pso.h.t[,5]=str_replace (annotation.pso.h.t[,5],"\\]","")
colnames(annotation.pso.h.t)[5]="PASI"

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t[,4]=as.character(annotation.pso.h.t[,4])
annotation.pso.h.t[,5]=as.character(annotation.pso.h.t[,5])
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family,annotation.pso.h.t$Psoriatic_arthritis,annotation.pso.h.t$PASI))
colnames(annotation.pso.h.t)[6]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[7]="sample"
ann$interactions=str_replace(ann$interactions,"early.yes","typeI")
ann$interactions=str_replace(ann$interactions,"late.no","typeII")
ann$interactions=str_replace(ann$interactions,"typeI.yes","typeI.PSA+")
ann$interactions=str_replace(ann$interactions,"typeII.yes","typeII.PSA+")
ann$interactions=str_replace(ann$interactions,"10_20","moderate")
ann$interactions=str_replace(ann$interactions,"0_10","mild")
keep=c("typeI.PSA+.moderate","typeII.PSA+.moderate","typeI.PSA+.mild","typeII.PSA+.mild")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[9]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_log_all)
tmp=tmp[tmp$interactions%in%keep,]
#plot groups 
colnames(tmp)[8]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression, color=interactions))+
  geom_point()+
  geom_boxplot()+
  theme_bw()

#make anova + turkey test
### ANOVA part on the selected genes ##############################
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
tmp2$Psoriatic_arthritis <- unname(tmp2$Psoriatic_arthritis)
tmp2 <- tmp2[,c(5,6,7)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

PS_type1_type2_psa=outDf
PS_type1_type2_psa=PS_type1_type2_psa[order(PS_type1_type2_psa$p_adj),]
PS_type1_type2_psa_signif=PS_type1_type2_psa[PS_type1_type2_psa$p_adj<0.05,]
save(PS_type1_type2_psa,PS_type1_type2_psa_signif, file = "PS_type1_type2_psa_anova.RData")
write.table(PS_type1_type2_psa_signif, file="PS_type1_type2_psa_signif.txt", quote=F, sep="\t")

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))



