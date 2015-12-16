#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND/anova")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
library(plyr)
load(file="data_for_anova.RData")#load selected_genes_blood,dt,data_tmp_or, annotation.pso from anova_blood_nd_part_I.R

##create correct annotation
################################################ Disease onset ##################
p=data_tmp_or[,71:173]
dim(p) 
################################################ Disease Onset ############################
annotation.do=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Sample","Gender","Age_of_disease_onset")]
annotation.do$Age_of_disease_onset=as.numeric(as.character(annotation.do$Age_of_disease_onset))
annotation.do$Gender=as.character(annotation.do$Gender)
str(annotation.do)
annotation.do[annotation.do[,3]>=40, 3]<-"late"
annotation.do[annotation.do[,3]<40,3]<-"early"
annotation.do[annotation.do[,3]==5,3]<-"early"
annotation.do[annotation.do[,3]==6,3]<-"early"
annotation.do=annotation.do[order(annotation.do[,3]),]

annotation.do=cbind(annotation.do, interactions=interaction(annotation.do$Sample, annotation.do$Age_of_disease_onset))
ann=cbind(annotation.do, sample=rownames(annotation.do))

#Merge annotation with the data, and from the merged data frame subset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"

tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)

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
pso_dis_onset_comp=outDf[outDf$Groups%in%c("72h.late-72h.early","SEB.late-SEB.early","PMA.late-PMA.early"),]
save(pso_dis_onset_comp, file = "Disease_onset_anova.RData")
pso_dis_onset_comp_signif=pso_dis_onset_comp[pso_dis_onset_comp$p_adj<0.05,]
write.table(pso_dis_onset_comp_signif, file="pso_dis_onset_comp_signif.txt", sep="\t", quote=F)
##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html
#95% confidence intervals visualizing the results of pairwise comparisons of the mean levels of gene expression
#grouped by skin state: lesional(PSL)/nonlesional(PSNL), and disease onset: < 40 years(early)/ > 40 years(late)
#The intervals that do not overlap with zero point correspond to statistically significant differences at P < 0.05.

ggplot(pso_dis_onset_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_pso_dis_onset_pso_fam_comp.png", width=10, height=10)
############################################################################################
################################################# Psoriasis in the family ##################
annotation.pif=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Gender","Psoriasis_in_family","Sample")]

annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"au.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"si.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"fa.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"mo.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"so.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"bro.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"chi.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"2_.*","familial")
annotation.pif$Psoriasis_in_family=str_replace(annotation.pif$Psoriasis_in_family,"no","sporadic")
annotation.pif=annotation.pif[order(annotation.pif$Psoriasis_in_family),]

annotation.pif=as.data.frame(annotation.pif)
annotation.pif[,1]=as.character(annotation.pif[,1])
annotation.pif[,2]=as.character(annotation.pif[,2])
annotation.pif[,3]=as.character(annotation.pif[,3])

annotation.pif=cbind(annotation.pif, interactions=interaction(annotation.pif$Sample, annotation.pif$Psoriasis_in_family))
ann=cbind(annotation.pif, rownames(annotation.pif))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
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
pso_fam_groups_comp=outDf[outDf$Groups%in%c("72h.sporadic-72h.familial","SEB.sporadic-SEB.familial","PMA.sporadic-PMA.familial"),]
save(pso_fam_groups_comp, file = "Psoriasis_in_family_anova.RData")
pso_fam_groups_comp_signif=pso_fam_groups_comp[pso_fam_groups_comp$p_adj<0.05,]
write.table(pso_fam_groups_comp_signif, file="pso_fam_groups_comp_signif.txt", sep="\t", quote=F)
ggplot(pso_fam_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_Psoriasis_in_family.png", width=10, height=10)

######################################### Psoriatic artritis##############################
annotation.pa=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Gender","Psoriatic_arthritis","Sample")]
annotation.pa=annotation.pa[order(annotation.pa$"Psoriatic_arthritis"),]

annotation.pa=as.data.frame(annotation.pa)
annotation.pa[,1]=as.character(annotation.pa[,1])
annotation.pa[,2]=as.character(annotation.pa[,2])
annotation.pa[,3]=as.character(annotation.pa[,3])

annotation.pa=cbind(annotation.pa, interactions= interaction(annotation.pa$Sample, annotation.pa$Psoriatic_arthritis))
ann=cbind(annotation.pa, sample=rownames(annotation.pa))

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
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
pso_arthritis_groups_comp=outDf[outDf$Groups%in%c("72h.yes-72h.no","SEB.yes-SEB.no","PMA.yes-PMA.no"),]
pso_arthritis_groups_comp_signif=pso_arthritis_groups_comp[pso_arthritis_groups_comp$p_adj<0.05,]
write.table(pso_arthritis_groups_comp_signif, file="pso_arthritis_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_arthritis_groups_comp, file = "Psoriatic_arthritis_in_family_anova.RData")
ggplot(pso_arthritis_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_psoriatic_arthritis_in_family.png", width=10, height=10)

########################################### Nail_involvment ##################
annotation.ni=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Gender","Nail_involvment","Sample")]
annotation.ni=annotation.ni[order(annotation.ni$Nail_involvment),]

annotation.ni=as.data.frame(annotation.ni)
annotation.ni[,1]=as.character(annotation.ni[,1])
annotation.ni[,2]=as.character(annotation.ni[,2])
annotation.ni[,3]=as.character(annotation.ni[,3])

annotation.ni=cbind(annotation.ni, interactions=interaction(annotation.ni$Sample, annotation.ni$Nail_involvment))
ann=cbind(annotation.ni, sample=rownames(annotation.ni))

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
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
pso_nail_inv_groups_comp=outDf[outDf$Groups%in%c("72h.yes-72h.no","SEB.yes-SEB.no","PMA.yes-PMA.no"),]
pso_nail_inv_groups_comp_signif=pso_nail_inv_groups_comp[pso_nail_inv_groups_comp$p_adj<0.05,]
write.table(pso_nail_inv_groups_comp_signif, file="pso_nail_inv_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_nail_inv_groups_comp, file = "Psoriasis_nail_involvement_anova.RData")
ggplot(pso_nail_inv_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  xlim(-6, 10)+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_Psoriasis_nail_involvement.png", width=10, height=10)

###################### PASI ########################################################################

annotation.pasi=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Gender","PASI_.activity_score.","Sample")]
annotation.pasi$PASI_.activity_score.=str_replace(annotation.pasi$PASI_.activity_score.,",",".")
annotation.pasi=data.frame(annotation.pasi)
annotation.pasi$PASI_.activity_score.=as.numeric(as.vector(annotation.pasi$PASI_.activity_score.))
annotation.pasi$PASI_.activity_score.=cut(annotation.pasi$PASI_.activity_score., c(0, 10, 20, Inf))

annotation.pasi=as.data.frame(annotation.pasi)
annotation.pasi[,1]=as.character(annotation.pasi[,1])
annotation.pasi[,2]=as.character(annotation.pasi[,2])
annotation.pasi[,3]=as.character(annotation.pasi[,3])
colnames(annotation.pasi)[3]="PASI"
annotation.pasi=cbind(annotation.pasi, interaction(annotation.pasi$Sample, annotation.pasi$PASI))
colnames(annotation.pasi)[4]="interactions"
ann=cbind(annotation.pasi, rownames(annotation.pasi))
colnames(ann)[5]="sample"
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
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
pso_pasi_groups_comp=outDf[outDf$Groups%in%c("72h.(10,20]-72h.(0,10]","72h.(20,Inf]-72h.(10,20]", "72h.(20,Inf]-72h.(0,10]","SEB.(10,20]-SEB.(0,10]","SEB.(20,Inf]-SEB.(0,10]","SEB.(20,Inf]-SEB.(10,20]",
                                             "PMA.(10,20]-PMA.(0,10]","PMA.(20,Inf]-PMA.(0,10]","PMA.(20,Inf]-PMA.(10,20]"),]
pso_pasi_groups_comp_signif=pso_pasi_groups_comp[pso_pasi_groups_comp$p_adj<0.05,]
write.table(pso_pasi_groups_comp_signif, file="pso_pasi_groups_comp_signif.txt", sep="\t", quote=F)

save(pso_pasi_groups_comp, file = "Psoriasis_PASI_anova.RData")
ggplot(pso_pasi_groups_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_Psoriasis_PASI.png", width=10, height=10)

############################Duration of the  psoriasis ##################

annotation.dp=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Gender","Duration_of_disease.years.","Sample")]
annotation.dp$Duration_of_disease.years.=str_replace(annotation.dp$Duration_of_disease.years.,",",".")
annotation.dp=data.frame(annotation.dp)
annotation.dp$Duration_of_disease.years.=as.numeric(as.vector(annotation.dp$Duration_of_disease.years.))
annotation.dp$Duration_of_disease.years.=cut(annotation.dp$Duration_of_disease.years., c(0, 5, 15, Inf))

annotation.dp=as.data.frame(annotation.dp)
annotation.dp[,1]=as.character(annotation.dp[,1])
annotation.dp[,2]=as.character(annotation.dp[,2])
annotation.dp[,3]=as.character(annotation.dp[,3])
colnames(annotation.dp)[3]="Duration_of_psoriasis"
annotation.dp=cbind(annotation.dp, interaction(annotation.dp$Sample, annotation.dp$Duration_of_psoriasis))
colnames(annotation.dp)[4]="interactions"
ann=cbind(annotation.dp, rownames(annotation.dp))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Duration_of_psoriasis <- unname(tmp2$Duration_of_psoriasis)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_durations_comp=outDf[outDf$Groups%in%c("72h.(15,Inf]-72h.(0,5]","72h.(5,15]-72h.(0,5]","72h.(5,15]-72h.(0,5]","SEB.(15,Inf]-SEB.(0,5]","SEB.(5,15]-SEB.(0,5]","SEB.(5,15]-SEB.(0,5]",
                                           "PMA.(15,Inf]-PMA.(0,5]","PMA.(5,15]-PMA.(0,5]","PMA.(5,15]-PMA.(0,5]"),]
pso_durations_comp_signif=pso_durations_comp[pso_durations_comp$p_adj<0.05,]
write.table(pso_durations_comp_signif, file="pso_durations_comp_signif.txt", sep="\t", quote=F)

save(pso_durations_comp, file = "Psoriasis_durations_comp_anova.RData")

ggplot(pso_durations_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_Psoriasis_durations_comp.png", width=10, height=10)


################################# Age >40, age<40 and psoriasis  in the family ########################
annotation.do.pin=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Age_of_disease_onset", "Psoriasis_in_family", "Sample")]
annotation.do.pin=annotation.do.pin[order(annotation.do.pin$Age_of_disease_onset),]
annotation.do.pin$Age_of_disease_onset=as.numeric(as.character(annotation.do.pin$Age_of_disease_onset))
annotation.do.pin$Psoriasis_in_family=as.character(annotation.do.pin$Psoriasis_in_family)

annotation.do.pin[annotation.do.pin[,2] >=40, 2]<-"late"
annotation.do.pin[annotation.do.pin[,2]<40,2]<-"early"
annotation.do.pin[annotation.do.pin[,2]==5, 2]<-"early"
annotation.do.pin[annotation.do.pin[,2] ==6,2]<-"early"


annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"au.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"si.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"fa.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"mo.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"so.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"bro.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"chi.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"2_.*","familial")
annotation.do.pin$Psoriasis_in_family=str_replace(annotation.do.pin$Psoriasis_in_family,"no","sporadic")
annotation.do.pin=annotation.do.pin[order(annotation.do.pin$Psoriasis_in_family),]

annotation.do.pin=as.data.frame(annotation.do.pin)
annotation.do.pin[,1]=as.character(annotation.do.pin[,1])
annotation.do.pin[,2]=as.character(annotation.do.pin[,2])
annotation.do.pin[,3]=as.character(annotation.do.pin[,3])
annotation.do.pin=cbind(annotation.do.pin, interaction(annotation.do.pin$Sample, annotation.do.pin$Age_of_disease_onset, annotation.do.pin$Psoriasis_in_family))
colnames(annotation.do.pin)[4]="interactions"
ann=cbind(annotation.do.pin, rownames(annotation.do.pin))
colnames(ann)[5]="sample"


#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
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
pso_dis_onset_pso_fam_comp=outDf[outDf$Groups%in%c("72h.early.sporadic-72h.early.familial","72h.late.sporadic-72h.late.familial",
                                                   "SEB.early.sporadic-SEB.early.familial","SEB.late.sporadic-SEB.late.familial",
                                                   "PMA.early.sporadic-PMA.early.familial","PMA.late.sporadic-PMA.late.familial"),]
pso_dis_onset_pso_fam_comp_signif=pso_dis_onset_pso_fam_comp[pso_dis_onset_pso_fam_comp$p_adj<0.05,]
write.table(pso_dis_onset_pso_fam_comp_signif, file="pso_dis_onset_pso_fam_comp_signif.txt", sep="\t", quote=F)

save(pso_dis_onset_pso_fam_comp, file = "Psoriasis_dis_onset_pso_fam_anova.RData")
ggplot(pso_dis_onset_pso_fam_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_Psoriasis_dis_onset_pso_fam.png", width=10, height=10)
################################### Age of disease onset <40 yr, >= 40 yr, PASI ######################

annotation.do.pasi=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("PASI_.activity_score.","Age_of_disease_onset","Sample")]

annotation.do.pasi$Age_of_disease_onset=as.numeric(as.character(annotation.do.pasi$Age_of_disease_onset))

annotation.do.pasi[annotation.do.pasi[,3] >=40, 3]<-"late"
annotation.do.pasi[annotation.do.pasi[,3]<40,3]<-"early"
annotation.do.pasi[annotation.do.pasi[,3]==5, 3]<-"early"
annotation.do.pasi[annotation.do.pasi[,3] ==6,3]<-"early"

annotation.do.pasi[,2]=str_replace(annotation.do.pasi[,2],",",".")
annotation.do.pasi=data.frame(annotation.do.pasi)
annotation.do.pasi[,2]=as.numeric(as.vector(annotation.do.pasi[,2]))
annotation.do.pasi[,2]=cut(annotation.do.pasi[,2], c(0, 10, 20, Inf))
annotation.do.pasi=annotation.do.pasi[order(annotation.do.pasi[,2]),]


annotation.do.pasi[,2]=str_replace (annotation.do.pasi[,2],"10,20","10_20")
annotation.do.pasi[,2]=str_replace (annotation.do.pasi[,2],"0,10","0_10")
annotation.do.pasi[,2]=str_replace (annotation.do.pasi[,2],"20,Inf","20_Inf")
annotation.do.pasi[,2]=str_replace (annotation.do.pasi[,2],"\\(","")
annotation.do.pasi[,2]=str_replace (annotation.do.pasi[,2],"\\]","")
colnames(annotation.do.pasi)[2]="PASI"


annotation.do.pasi=as.data.frame(annotation.do.pasi)
annotation.do.pasi[,1]=as.character(annotation.do.pasi[,1])
annotation.do.pasi[,2]=as.character(annotation.do.pasi[,2])
annotation.do.pasi[,3]=as.character(annotation.do.pasi[,3])
annotation.do.pasi=cbind(annotation.do.pasi, interaction(annotation.do.pasi$Sample, annotation.do.pasi$Age_of_disease_onset, annotation.do.pasi$PASI))
colnames(annotation.do.pasi)[4]="interactions"
ann=cbind(annotation.do.pasi, rownames(annotation.do.pasi))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)
colnames(tmp)[6]="Gene"
ggplot(tmp, aes(x=Gene,y=expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
ggsave("Jitter_Psoriasis_dis_onset_pasi.png", width=10, height=10)
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
outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

#because of small number od samples in PL.late, PNL.late all subgroups, additional SEBot is produced for PL.early, PNL.early all groups comparisons.
pso_dis_onset_pasi_comp_early=outDf[outDf$Groups%in%c("72h.early.20_Inf-72h.early.10_20","72h.early.20_Inf-72h.early.0_10","72h.early.10_20-72h.early.0_10",
                                                      "SEB.early.20_Inf-SEB.early.10_20","SEB.early.20_Inf-SEB.early.0_10","SEB.early.10_20-SEB.early.0_10",
                                                      "PMA.early.20_Inf-PMA.early.10_20","PMA.early.20_Inf-PMA.early.0_10","PMA.early.10_20-PMA.early.0_10"),]

pso_dis_onset_pasi_comp_signif_early=pso_dis_onset_pasi_comp_early[pso_dis_onset_pasi_comp_early$p_adj<0.05,]
write.table(pso_dis_onset_pasi_comp_signif_early, file="pso_dis_onset_pasi_comp_signif_early.txt", sep="\t", quote=F)

save(pso_dis_onset_pasi_comp_early, file = "Psoriasis_dis_onset_pasi_anova_early.RData")

ggplot(pso_dis_onset_pasi_comp_early, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_Psoriasis_dis_onset_pasi.png", width=10, height=10)

################################# Disease onset >40, <40 and nail involvment ########################
annotation.do.ni=annotation.pso[rownames(annotation.pso)%in%colnames(p)[1:103],colnames(annotation.pso)%in%c("Age_of_disease_onset", "Nail_involvment","Sample")]
annotation.do.ni=annotation.do.ni[order(annotation.do.ni[,2]),]
annotation.do.ni$Age_of_disease_onset=as.numeric(as.character(annotation.do.ni$Age_of_disease_onset))
annotation.do.ni[annotation.do.ni[,2]>=40, 2]<-"late"
annotation.do.ni[annotation.do.ni[,2]< 40,2]<-"early"
annotation.do.ni[annotation.do.ni[,2]==5, 2]<-"early"
annotation.do.ni[annotation.do.ni[,2] ==6,2]<-"early"

annotation.do.ni=as.data.frame(annotation.do.ni)
annotation.do.ni[,1]=as.character(annotation.do.ni[,1])
annotation.do.ni[,2]=as.character(annotation.do.ni[,2])
annotation.do.ni[,3]=as.character(annotation.do.ni[,3])
annotation.do.ni=cbind(annotation.do.ni, interaction(annotation.do.ni$Sample, annotation.do.ni$Age_of_disease_onset, annotation.do.ni$Nail_involvment))
colnames(annotation.do.ni)[4]="interactions"
ann=cbind(annotation.do.ni, rownames(annotation.do.ni))
colnames(ann)[5]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes_blood)

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Nail_involvment <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(4,5,6)]

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"
pso_dis_onset_nail_inv_comp=outDf[outDf$Groups%in%c("72h.early.yes-72h.early.no","72h.late.yes-72h.late.no",
                                                    "SEB.early.yes-SEB.early.no","SEB.late.yes-SEB.late.no",
                                                    "PMA.early.yes-PMA.early.no","PMA.late.yes-PMA.late.no"),]
pso_dis_onset_nail_inv_comp_signif=pso_dis_onset_nail_inv_comp[pso_dis_onset_nail_inv_comp$p_adj<0.05,]
write.table(pso_dis_onset_nail_inv_comp_signif, file="pso_dis_onset_nail_inv_comp_signif.txt", sep="\t", quote=F)

save(pso_dis_onset_nail_inv_comp, file = "Psoriasis_dis_onset_nail_inv_anova.RData")

ggplot(pso_dis_onset_nail_inv_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_Psoriasis_dis_onset_nail_inv.png", width=10, height=10)
