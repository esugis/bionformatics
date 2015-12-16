
setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/anova")
library("limma")


#####################################################################################################################################
############################################ Select diffexp genes in log-transformed data for further analysis  ############################
##Since the changes can be already in non-leasional skin, we first select the genes that are diff expressed in (PSL+PSNL) vs C

load(file = "data_log_pso_remaned_samples.RData")
rownames(data_log)[24]="CXCL8"
pheatmap(data_log, cluster_rows = F, cluster_cols = F, scale = "none", main="Gene expression in control group and psoriatic patients")
summary(data_log)
data=data_log
##diff expressed PSL+PSNL vs C
t=data
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PSL","PSNL")
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

############################################# ANOVA part on the selected genes in log-transformed data ##############################

#load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/tmp_melted_dt.RData")
#tmp <- subset(melted_dt,variable%in%selected_genes_log_all)
#tmp2 <- tmp[,-c(1)]
#tmp2$Sample <- unname(tmp2$Sample)
#tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
#tmp2 <- tmp2[,c(3,4,5)]

##Aov example for many genes
#myAov <- function(df){
#  aov_list <- aov(expression~interactions, data=df)
#  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
#  aov_out$groups <- rownames(aov_out)
#  return(aov_out)
#}

##test <- subset(tmp2, variable=="WIPI1")
##myAov(test)
#outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
#colnames(outDf)[5] <- "p_adj"
#dis_onset_groups_comp=outDf
#save(dis_onset_groups_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/Psoriasis_disease_onset_group_comparison.RData")

###faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

#ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
#  geom_point()+
#  theme_bw()+
#  facet_wrap(~variable)+
#  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
#  ggtitle("Comparison of gene expression in different groups using log-transformed data")



#######################################################################################################
###################Test on nonlogtransformed data converted to Z scores ###############################
library(stringr)
library(reshape)
library(ggplot2)
library(pheatmap)
load("~/Documents/KAI/cytoskin/after_filt/21112014/data_all_filt_nontransformed.RData")
rownames(data_all_filt)[24]="CXCL8"
data_all_filt=data_all_filt[, c(1:22,24:118)]
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"C","C")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"PH","PSL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"PT","PSNL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"VH","VL")
colnames(data_all_filt)=str_replace(colnames(data_all_filt),"VT","VNL")
save(data_all_filt, file = "data_zscores_filt_remaned_samples.RData")
data_raw=data_all_filt[,1:86]
pheatmap(data_raw, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controlls and Psoriasis")

pheatmap(data_raw, cluster_rows = F, cluster_cols = F, scale = "row", main="Data Controlls and Psoriasis")


summary(data_raw)
data=scale_rows(data_raw)
summary(data)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controlls and Psoriasis")

###diff expressed PH+PT vs C
#t=data
#cols=colnames(t)
#attribute=colnames(t)
#attribute=str_replace(attribute,"0.*","")
#coldata=colnames(t)
#coldata=str_replace(coldata,"0.*","")
#attribute = factor(attribute)
#attribute2 = rep("Muu", length(attribute))
## Case-control variables and corresponding group values 
#case.vals=c("PH","PT")
#cont.vals =c("C")
#attribute2[attribute %in% c(case.vals)] = "Group1"
#attribute2[attribute %in% c(cont.vals)] = "Group2"
#attribute = factor(attribute2)
#colnames(t)=coldata
#head(t)
#mm = model.matrix(~ attribute - 1)
#colnames(mm) = make.names(colnames(mm))
#colnames(mm)
#contrast = makeContrasts(attributeGroup2 - attributeGroup1, levels = colnames(mm))
#fit = lmFit(t, mm)
#fit = contrasts.fit(fit, contrast)
#fit = eBayes(fit)
#adj.method="fdr" 
#thr.p.value=0.05
#thr.lfc=1
#tt_z_lfc = toptable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
#tt_all_z=toptable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)#diff expressed genes filtered only based on p-val <0.05
#print(thr.p.value)
#print(tt) 
#selected_genes_z=rownames(tt_lfc)#7
#selected_genes_all_z=rownames(tt_all_z)
#length(selected_genes_all_z)#28
#save(tt_z_lfc,tt_all_z,selected_genes_z,selected_genes_all_z, file = "~/Documents/KAI/cytoskin/after_filt/21112014/diffexp_zscores.RData")


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

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t=cbind(annotation.pso.h.t, interaction(annotation.pso.h.t$Sample, annotation.pso.h.t$Age_of_disease_onset))
colnames(annotation.pso.h.t)[4]="interactions"
ann=cbind(annotation.pso.h.t, rownames(annotation.pso.h.t))
colnames(ann)[5]="sample"
data=as.data.frame(data)
dt=t(data)
dt=dt[25:87,]
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[36]="sample"
dt[,36]=as.character(dt[,36])
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="expression"

################################# ANOVA Z-score data ##############
tmp_z <- subset(dt.merged.melt,variable%in%selected_genes_wlx)
tmp2 <- tmp_z[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2 <- tmp2[,c(4,5,6)]
colnames(tmp2)[1]="interactions"

outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
dis_onset_groups_comp2=outDf
save(dis_onset_groups_comp2, file = "~/Documents/KAI/cytoskin/after_filt/21112014/Psoriasis_disease_onset_group_comparison2.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+geom_point()+theme_bw()+
  facet_wrap(~variable)+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of gene expression in different groups. Gene expression data was converted to Z-scores")



###########################################################################################################################
################# Compare differentially expressed genes in log-transformed and Z-score transformed data data ##############

#length(selected_genes_log)#21 limma diffexp genes for logtransformed data using lfc parameter
length(selected_genes_log_all)#27 limma diffexp genes for logtransformed data using only p-val parameter filt

length(selected_genes_wlx)#27 wilcoxon test diffexp genes for z-score data
length(selected_genes_all_z)#28 limma diffexp genes for z-score data

#Comparion of diffexp genes obtained in log transformed data and z-score data
selected_genes_wlx[selected_genes_wlx%in%selected_genes_log_all]#27 genes=> selected_genes_log_all contains the same genes as selected_genes_wlx

#[1] "AIM2"    "CCL2"    "CCL20"   "CCL27"   "CCL5"    "CTLA4"   "CXCL1"   "CXCL10"  "CXCL2"   "DEFB1"   "EOMES"   "FOXP3"   "IFIH1"  
#[14] "IFNG"    "IFNGR"   "IL1b"    "IL1F6"   "IL1RN"   "IL22RA1" "IL8b"    "LCN2"    "MICB"    "PI3"     "PYCARD"  "S100A8"  "S100A9" 
#[27] "TNF"

#### Comparison for diffexp genes in z-score data obtained by limma and wilcoxon test
selected_genes_wlx[!(selected_genes_wlx%in%selected_genes_all_z)]#[1] "CCL20"
selected_genes_all_z[!(selected_genes_all_z%in%selected_genes_wlx)]#[1] "IFNAR1" "IL20RA"

#### Q? Shall I take overlap of the results found by both methods or prefer one of them?
#### As wilcoxon test and limma for log transformed data gave the same result,
#### this set of diffexp genes will be used in the analysis

#############take this 23.03.2015
############################## Show the comparisons of distributions for raw nonscaled, Z-score converted and log transformed#####
#Plot distributions for log-transformed data in PH and Pt conditions
ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed data.")

ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_histogram(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Histogram of gene expression in log-transformed data..")

##########################skip for now,not sure if needed#########
#ggplot(tmp, aes(x=expression))+
#  facet_wrap(~variable)+
#  theme_bw()+
#  geom_density(aes(group=interactions, colour=interactions,fill=interactions), alpha=0.5)+
#  ggtitle("Density plot of gene expression in log-transformed data. Data distribution is shown in all coditions.")



##Plot distributions for Z-scored data in PH and Pt conditions
#ggplot(tmp_z, aes(x=expression))+
#  facet_wrap(~variable)+
#  theme_bw()+
#  # log from negative value scale_x_log10()+
#  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
#  ggtitle("Density plot of gene expression in data converted to Z-scores.")

#ggplot(tmp_z, aes(x=expression))+
#  facet_wrap(~variable)+
#  theme_bw()+
#  # log from negative value scale_x_log10()+
#  geom_histogram(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
#  #aes(y = ..density..)+
#  ggtitle("Histogram of gene expression in data converted to Z-scores.")



##Plot distributions for raw data in PH and Pt conditions
#data_raw=data_all_filt[,1:87]
#data_raw=as.data.frame(data_raw)
#dt_raw=t(data_raw)
#dt_raw=dt_raw[25:87,]
#dt_raw=as.data.frame(dt_raw)
#dt_raw=cbind(dt_raw, rownames(dt_raw))
#colnames(dt_raw)[36]="sample"
#dt_raw[,36]=as.character(dt_raw[,36])
#dt.merged_raw=merge(ann,dt_raw,by="sample",all=T)
#dt.merged.melt_raw=melt(dt.merged_raw)
#colnames(dt.merged.melt_raw)[7]="expression"
#tmp_raw <- subset(dt.merged.melt_raw,variable%in%selected_genes_wlx)


#ggplot(tmp_raw, aes(x=expression))+
#  facet_wrap(~variable)+ 
#  geom_density(aes(group=interactions, colour=interactions,fill=interactions), alpha=0.5, na.rm=TRUE)+
#  scale_x_log10()+
#  theme_bw()+
#  ggtitle("Density plot of gene expression in raw data. Data distribution is shown in all coditions.")

##plot(density(tmp_raw[tmp_raw$variable%in%c("CCL27"),7 ]))





######################################################################################
###################################example with smaller number og genes#############
#targets$samples <- rownames(targets)
#dft$samples <- rownames(dft)
#merged_dt <- merge(targets, dft, by = "samples",all.x=TRUE)
#save(dft,merged_dt,melted_dt,tmp, tmp2, file = "~/Documents/KAI/cytoskin/after_filt/21112014/tmp_melted_dt.RData")
#library(reshape2)
#melted_dt <- melt(merged_dt, id.vars=1:3, value.name="expression")
#melted_dt$interactions <- interaction(melted_dt$Sample,melted_dt$Age_of_disease_onset)

##model_lm2 <- lm(data=melted_dt, expression~interactions)
##model_lm <- lm(data=melted_dt, expression ~ Sample*Age_of_disease_onset)

#tmp <- subset(melted_dt,variable=="TRGC1"| variable=="AIM2"| variable=="WIPI1")
#tmp2 <- tmp[,-c(1)]
#tmp2$Sample <- unname(tmp2$Sample)
#tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
#tmp2 <- tmp2[,c(3,4,5)]

##Aov example for many genes
#myAov <- function(df){
#  aov_list <- aov(expression~interactions, data=df)
#  aov_out <- as.data.frame(TukeyHSD(aov_list, "interactions")$interactions)
#  aov_out$groups <- rownames(aov_out)
#  return(aov_out)
#}

##test <- subset(tmp2, variable=="WIPI1")
##myAov(test)
#outDf <- ddply(tmp2, "variable", function(tmp2)  myAov(tmp2))
#colnames(outDf)[5] <- "p_adj"


###faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

#ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+geom_point()+theme_bw()+
#  facet_wrap(~variable)+
#  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))



###aov example for 1 gene
##tmp_aov <- aov(data=tmp2, expression~interactions)
##ac <- as.data.frame(TukeyHSD(tmp_aov, "interactions",ordered=T)$interactions)
##ac$groups <- rownames(ac)

###linear model
##lm3 <- lm(data=tmp, expression ~ Age_of_disease_onset)
##sum(lm3$coefficients)#nondiseased, young



#myLm <- function(formula, df){
#  lmList <- lm(formula, data=df)
#  lmOut <- data.frame(t(lmList$coefficients))
#  #names(lmOut) <- c("intercept","Age_of_disease_onset_coef","Sample_coef")
#  return(lmOut)
#}

#outDf <- ddply(tmp2, "variable", function(tmp2)  myLm(expression ~ Age_of_disease_onset + Sample, tmp2))