#Metaanalysis psoriasis
setwd("~/Documents/KAI/cytoskin/nondetects/SKIN_ND/anova")
#load log2transformed data
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
library("limma")
library(plyr)
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
#tmp=data.frame(tmp,stringsAsFactors=FALSE) 
#colnames(tmp)[40]="expression"
#tmp2 <- tmp[,-c(1)]
#tmp2$Sample <- unname(tmp2$Sample)
#tmp2$Nail_involvment <- unname(tmp2$Nail_involvment)
#tmp2 <- tmp2[,c(4,5,6)]
lmskin <- lm(data=tmp, expression ~ Age_of_disease_onset+ PASI_.activity_score.+Duration_of_disease.years.+Skin_phototype+Psoriatic_arthritis+Psoriasis_in_family+ Nail_involvment)
sum(lmskin$coefficients)
confint(lmskin, level=0.95) # CIs for model parameters 
fitted(lmskin) # predicted values

residuals(lmskin) # residuals

anova(lmskin) # anova table 
vcov(lmskin) # covariance matrix for model parameters 
influence(lmskin)
lmOut <- data.frame(t(lmskin$coefficients))
summary(lmskin)
plot(lmskin)
#pl 
tmp_pl=tmp[tmp$Sample%in%"PL",]
lmskinpl <- lm(data=tmp_pl, expression ~ Age_of_disease_onset+ PASI_.activity_score.+Duration_of_disease.years.+Skin_phototype+Psoriatic_arthritis+Psoriasis_in_family+ Nail_involvment)
summary(lmskinpl)
tmp_pnl=tmp[tmp$Sample%in%"PNL",]
lmskinpnl <- lm(data=tmp_pnl, expression ~ Age_of_disease_onset+ PASI_.activity_score.+Duration_of_disease.years.+Skin_phototype+Psoriatic_arthritis+Psoriasis_in_family+ Nail_involvment)
summary(lmskinpnl)
# Influential Observations
# added vainstariable plots 
av.Plots(lmskin)
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(tmp)-length(lmskin$coefficients)-2)) 
plot(lmskin, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(lmskin,  id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

###
myLm <- function(formula, df){
  lmList <- lm(formula, data=df)
  lmOut <- data.frame(t(lmList$coefficients))
  names(lmOut) <- c("intercept","Age_of_disease_onset_coef","Sample_coef",
 #                    "PASI_.activity_score._coef","Duration_of_disease.years._coef","Skin_phototype_coef",
                    "Psoriatic_arthritis_coef","Psoriasis_in_family_coef","Nail_involvment_coef")
  return(lmOut)
}

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






#################tree 
library(rpart)
# grow tree 
tmp2=tmp[which(complete.cases(tmp[40])),]
tmp2_pl=tmp2[tmp2$Sample%in%"PL",]
#tmp2=tmp2[,colnames(tmp2)%in%c("sample","Sample","Age_of_disease_onset", "PASI_.activity_score.","Duration_of_disease.years.","Skin_phototype","Psoriatic_arthritis","Psoriasis_in_family", "Nail_involvment","expression")]
fit <- rpart(expression ~ Age_of_disease_onset+ + PASI_.activity_score.+Duration_of_disease.years.+Skin_phototype+Psoriatic_arthritis+Psoriasis_in_family+ Nail_involvment, 
             method="anova", data=tmp2_pl)
fit <- rpart(expression ~ Psoriatic_arthritis, method="anova", data=tmp2)
printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

# create additional plots 
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(fit) # visualize cross-validation results    

# plot tree 
plot(fit, uniform=TRUE, 
     main="Regression Tree ")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postcript plot of tree 
# prune the tree 
pfit<- prune(fit, cp=0.01160389) # from cptable   

# plot the pruned tree 
plot(pfit, uniform=TRUE, 
     main="Pruned Regression Tree")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)






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
ggsave("ANOVA_Psoriasis_dis_onset_nail_inv.png", width=15, height=10)