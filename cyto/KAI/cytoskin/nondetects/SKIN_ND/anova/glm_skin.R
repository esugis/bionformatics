#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND/anova")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
library(plyr)
library(glmulti)
#load the gene expression  data 
load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
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

##create correct annotation
################################################ Disease onset ##################
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

################################################ Disease Onset ############################
annotation.h=data.frame(annotationh,stringsAsFactors=FALSE)
annotation.h[annotation.h$Age_of_disease_onset>=40, 9]<-"late"
annotation.h[annotation.h$Age_of_disease_onset< 40,9]<-"early"
annotation.h=annotation.h[order(annotation.h$Age_of_disease_onset),]
#Psoriasis_in_family
annotation.h[,27]=str_replace(annotation.h[,27],"2_.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"au.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"si.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"fa.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"mo.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"so.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"bro.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"chi.*","familial")
annotation.h[,27]=str_replace(annotation.h[,27],"no","sporadic")
annotation.h=annotation.h[order(annotation.h[,27]),]
######################################### Psoriatic artritis##############################

########################################### Nail_involvment ##################

###################### PASI ########################################################################
annotation.h[,7]=str_replace(annotation.h[,7],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,7]=as.numeric(as.vector(annotation.h[,7]))
annotation.h[,7]=cut(annotation.h[,7], c(0, 10, 20, Inf))


#################################### Skin phototype ###########################################

###############################Duration of the  psoriasis ##################
annotation.h[,10]=str_replace(annotation.h[,10],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,10]=as.numeric(as.vector(annotation.h[,10]))
annotation.h[,10]=cut(annotation.h[,10], c(0, 5, 15, Inf))

################################# merge with expression ########################
data=as.data.frame(data.skin)
p=data[,24:87]
dt=t(data)
dt=cbind(sample=row.names(dt),dt)
dt=dt[24:87,]
annotation.h=annotation.h[rownames(annotation.h)%in%colnames(p)[1:64],]
annotation.t=annotation.h
rownames(annotation.t)=str_replace(rownames(annotation.t),"PL","PNL")

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"PL0.*","PL")
Sample=str_replace(Sample,"PNL0.*","PNL")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t,stringsAsFactors=FALSE)
#annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Nail_involvment))
#colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[38]="sample"

#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged[, 39:80] <- sapply(dt.merged[, 39:80], as.character)
dt.merged[, 39:80] <- sapply(dt.merged[, 39:80], as.numeric)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[40]="expression"
load(file="selected.genes.RData")

tmp <- subset(dt.merged.melt,variable%in%selected.genes)

res <- glmulti(OUTPUT ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat,
               level=1, crit="aicc”)






outDf <- ddply(tmp, "variable", function(tmp)  myLm (expression ~ Age_of_disease_onset + Sample +
                                                       PASI_.activity_score.+ Duration_of_disease.years. +
                                                       Skin_phototype + Psoriatic_arthritis +
                                                       Psoriasis_in_family +Nail_involvment, tmp))
outDf <- ddply(tmp, "variable", function(tmp) lm (data=tmp, expression ~ Age_of_disease_onset + Sample +
                                                    PASI_.activity_score.+ Duration_of_disease.years. +
                                                    Skin_phototype + Psoriatic_arthritis +
                                                    Psoriasis_in_family +Nail_involvment))

#for one gene

lmList <- lm(expression ~ Age_of_disease_onset + Sample +
               PASI_.activity_score.+ Duration_of_disease.years. +
               Skin_phototype + Psoriatic_arthritis +
               Psoriasis_in_family +Nail_involvment, data=tmp[tmp$variable%in%"AIM2",])
lmOut <- data.frame(t(lmList$coefficients))
names(lmOut) <- c("intercept","Age_of_disease_onset_coef","Sample_coef",
                  #                    "PASI_.activity_score._coef","Duration_of_disease.years._coef","Skin_phototype_coef",
                  "Psoriatic_arthritis_coef","Psoriasis_in_family_coef","Nail_involvment_coef")








ggplot(pso_dis_onset_nail_inv_comp, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_Psoriasis_dis_onset_nail_inv.png", width=15, height=10)