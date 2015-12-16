#beginning is in anova_nd_blood_part_II.R
################################# 1.Age >40, age<40 and psoriasis  in the family ########################
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
#Keep only the following groups
#72h.early.yes
#SEB.early.yes
#PMA.early.yes
#72h.late.no
#SEB.late.no
#PMA.late.no
keep=c("72h.early.familial","SEB.early.familial","PMA.early.familial","72h.late.sporadic","SEB.late.sporadic","PMA.late.sporadic")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
ann$interactions=str_replace(ann$interactions,"early.familial","typeI")
ann$interactions=str_replace(ann$interactions,"early.familial","typeI")
ann$interactions=str_replace(ann$interactions,"late.sporadic","typeII")
ann$interactions=str_replace(ann$interactions,"late.sporadic","typeII")

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

PS_typeI_typeII=outDf[outDf$Groups%in%c("72h.typeII-72h.typeI","SEB.typeII-SEB.typeI","PMA.typeII-PMA.typeI"),]
PS_typeI_typeII_signif=PS_typeI_typeII[PS_typeI_typeII$p_adj<0.05,]
save(PS_typeI_typeII,PS_typeI_typeII_signif, file = "PS_typeI_typeII_anova.RData")
write.table(PS_typeI_typeII_signif, file="PS_typeI_typeII.txt", sep="\t", quote=F)

ggplot(PS_typeI_typeII, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_PS_typeI_typeII.png", width=10, height=10)


#############2. Early onset. High PASI (>20) Late onset. Low PASI (<10)################
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
#annotation.do.pasi=cbind(annotation.do.pasi, interaction(annotation.do.pasi$Sample, annotation.do.pasi$Age_of_disease_onset, annotation.do.pasi$PASI))
annotation.do.pasi=cbind(annotation.do.pasi, interaction(annotation.do.pasi$Age_of_disease_onset, annotation.do.pasi$PASI))

colnames(annotation.do.pasi)[4]="interactions"
ann=cbind(annotation.do.pasi, rownames(annotation.do.pasi))
colnames(ann)[5]="sample"

#keep=c("72h.early.20_Inf","72h.late.0_10",
#       "SEB.early.20_Inf","SEB.late.0_10",
#       "PMA.early.20_Inf","PMA.late.0_10")
keep=c("early.20_Inf","late.0_10")
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt


#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)
colnames(dt.merged.melt)[7]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes)
tmp=tmp[tmp$interactions%in%keep,]
#plot groups 
colnames(tmp)[6]="Gene"
tmp$interactions=str_replace(tmp$interactions,"20_Inf","severe")
tmp$interactions=str_replace(tmp$interactions,"0_10","mild")
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
ggsave("jitter_ps_early.severe_late.mild_2groups.png", width=10, height=10)

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

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_ps_early.severe_late.mild_2groups.png", width=10, height=10)

############################### 3. PL.early.severe PNL.early.severe;PL.late.severe  PNL.late.severe

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
keep=c("72h.early.20_Inf","SEB.early.20_Inf","PMA.early.20_Inf",
       "72h.late.20_Inf","SEB.late.20_Inf","PMA.late.20_Inf")
      
ann_filt=ann[ann$interactions%in%keep,]
ann=ann_filt
#Merge annotation with the data andsubset differentially expressed genes
dt.merged=merge(ann,dt,by="sample",all=T)
dt.merged.melt=melt(dt.merged)

colnames(dt.merged.melt)[7]="Expression"
tmp <- subset(dt.merged.melt,variable%in%selected_genes)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[6]="Gene"
tmp$interactions=str_replace(tmp$interactions,"20_Inf","severe")
tmp$interactions=str_replace(tmp$interactions,"0_10","mild")
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
ggsave("jitter_early_late_severe.png", width=10, height=10)
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
tmp2=tmp2[!tmp2$Gene%in%c("IL17F"),]
outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))

colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

early_late_severe=outDf[outDf$Groups%in%c("SEB.late.severe-SEB.early.severe","PMA.late.severe-PMA.early.severe","72h.late.severe-72h.early.severe"),]
early_late_severe_signif=early_late_severe[early_late_severe$p_adj<0.05,]
save(early_late_severe,early_late_severe_signif, file = "early_late_severe_anova.RData")
write.table(early_late_severe_signif, file="early_late_severe_signif.txt", quote=F, sep="\t")

ggplot(early_late_severe, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))
ggsave("ANOVA_early_late_severe.png", width=10, height=10, dpi=100)


########################## 4. Type I + pos.PS Arthritis Type I + neg PS Arthritis; Type II + pos. PS Arthritis Type II + neg. PS Arthritis 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:64],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.h=annotation.h[,c(1,3,2)]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:64],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis")]
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
Sample=str_replace(Sample,"PL0.*","PL")
Sample=str_replace(Sample,"PNL0.*","PNL")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

annotation.pso.h.t=as.data.frame(annotation.pso.h.t)
annotation.pso.h.t[,1]=as.character(annotation.pso.h.t[,1])
annotation.pso.h.t[,2]=as.character(annotation.pso.h.t[,2])
annotation.pso.h.t[,3]=as.character(annotation.pso.h.t[,3])
annotation.pso.h.t[,4]=as.character(annotation.pso.h.t[,4])
annotation.pso.h.t=cbind(annotation.pso.h.t,interactions= interaction(annotation.pso.h.t$Age_of_disease_onset, annotation.pso.h.t$Psoriasis_in_family,annotation.pso.h.t$Psoriatic_arthritis))

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
tmp <- subset(dt.merged.melt,variable%in%selected_genes)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[7]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
ggsave("Gene_expression_PS_type1_type2_psa.png", width=20, height=10, dpi=100)

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
str(PS_type1_type2_psa)
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
ggsave("ANOVA_PS_type1_type2_psa.png", width=15, height=10, dpi=100)


########################## 5 Type I + NI(+) Type I + NI(-); Type II + NI(+) Type II + NI(-) 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:64],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"
annotation.h=annotation.h[,c(1,3,2)]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:64],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Nail_involvment")]
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
Sample=str_replace(Sample,"PL0.*","PL")
Sample=str_replace(Sample,"PNL0.*","PNL")
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
tmp <- subset(dt.merged.melt,variable%in%selected_genes)
tmp=tmp[tmp$interactions%in%keep,]
colnames(tmp)[7]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw()
ggsave("Gene_exp_PS_type1_type2_ni.png", width=, height=10, dpi=100)

tmp2 <- tmp[,-c(1)]
tmp2$Sample <- unname(tmp2$Sample)
tmp2$Age_of_disease_onset <- unname(tmp2$Age_of_disease_onset)
tmp2$Psoriasis_in_family <- unname(tmp2$Psoriasis_in_family)
tmp2$Psoriatic_arthritis <- unname(tmp2$Nail_involvment)
tmp2 <- tmp2[,c(5,6,7)]

outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
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
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))

ggsave("ANOVA_PS_type1_type2_ni.png", width=15, height=10, dpi=100)

#######  Type I + PSA(+) + PASI moderate (10-20), Type I + PSA(+) + PASI mild ( <10), Type II + PSA(+) + PASI moderate (10-20), Type II + PSA(+) + PASI mild ( <10) 
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:64],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis","PASI_(activity_score)")]
annotation.h=annotation.h[,c(2,4,3,1)]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"late"
annotation.h[annotation.h[,1]< 40,1]<-"early"

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:64],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family","Psoriatic_arthritis","PASI_(activity_score)")]
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
Sample=str_replace(Sample,"PL0.*","PL")
Sample=str_replace(Sample,"PNL0.*","PNL")
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
tmp <- subset(dt.merged.melt,variable%in%selected_genes)
tmp=tmp[tmp$interactions%in%keep,]
#plot groups 
colnames(tmp)[8]="Gene"
ggplot(tmp, aes(x=Gene,y=Expression,color=interactions))+
  geom_jitter(aes(shape = interactions),position = position_jitter(width = .21),size = 2.5)+
  theme_bw() + ggsize(width=12, height=10, dpi=100)

ggsave("Gene_exp_PS_type1_type2_psa_mild_moderate.png", width=20, height=10)

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
tmp2 <- tmp2[,c(6,7,8)]

outDf <- ddply(tmp2, "Gene", function(tmp2)  myAov(tmp2))
colnames(outDf)[5] <- "p_adj"
colnames(outDf)[2] <- "Difference"
colnames(outDf)[6] <- "Groups"

PS_type1_type2_psa=outDf
PS_type1_type2_psa=PS_type1_type2_psa[order(PS_type1_type2_psa$p_adj),]
PS_type1_type2_psa_signif=PS_type1_type2_psa[PS_type1_type2_psa$p_adj<0.05,]
#save(PS_type1_type2_psa,PS_type1_type2_psa_signif, file = "PS_type1_type2_psa_mild_moderate_anova.RData")
#write.table(PS_type1_type2_psa_signif, file="PS_type1_type2_psa_pasi_mild_moderate_signif.txt", quote=F, sep="\t")

ggplot(outDf, aes(x=Difference,y=Groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=9))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))


#ggsave("ANOVA_exp_PS_type1_type2_psa_pasi_mild_moderate.png", width=15, height=10, dpi=100)

