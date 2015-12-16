library("limma")
library(stringr)
library(pheatmap)
setwd("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/diffexp/pheno")
load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")
data_filt=data_log2_filt
##Controlls
ctrls=data_filt[,1:24]

##PSO.h
ph=data_filt[,25:56]

##PSO.t
pt=data_filt[,57:87]

##VIT.h
vh=data_filt[,88:103]

##VIT.t
vt=data_filt[,104:118]

###########################Put Controlls, Pso, Vit in one matrix##########
data=cbind(ctrls,ph,pt,vh,vt)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "row", main="All data after standardization")
summary(data)

#!!!!!Change metada!!!!


#annotation
annotation=t(p.meta)
annotationh=annotation
annotationt=annotation
rownames(annotationh)=str_replace(rownames(annotation),"P","PH")
rownames(annotationt)=str_replace(rownames(annotation),"P","PT")


#p.test=data[,25:87]
#p=p.test
#cols=colnames(p)
#colnames(p)=cols

#Test skin phototype influence
#annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Skin_phototype")]
#annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Skin_phototype")]

annotation.h=annotationh[rownames(annotationh)%in%colnames(data)[25:56],colnames(annotationh)%in%c("Gender","Skin_phototype")]
annotation.t=annotationt[rownames(annotationt)%in%colnames(data)[57:87],colnames(annotationt)%in%c("Gender","Skin_phototype")]


rownames(annotationh)=str_replace(rownames(annotation),"P","PH")
rownames(annotationt)=str_replace(rownames(annotation),"P","PT")
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)



#targets=data.frame(annotation.pso.h.t[,c(1,3)])
#pheno=paste(targets$Sample, targets$Skin_phototype, sep=".")
#pheno=factor(pheno, levels=pheno)
#design <- model.matrix(~0 +pheno)
#colnames(design)
#colnames(design)=str_replace(colnames(design),"pheno","")
##colnames(design) <- levels(pheno)
#fit <- lmFit(p, design)
##1. Type I
#contrast = makeContrasts(PHvsPTinI=PH.I-PT.I, levels=colnames(design))
#fit2 <- contrasts.fit(fit, contrast)
#fit2 <- eBayes(fit2)
#topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

##results <- decideTests(fit2)
##vennDiagram(results)
##2. type II
#contrast = makeContrasts(PHvsPTinII=PH.II-PT.II, levels=colnames(design))
#fit2 <- contrasts.fit(fit, contrast)
#fit2 <- eBayes(fit2)
#topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

##3 Type III
#contrast = makeContrasts(PHvsPTinIII=PH.III-PT.III, levels=colnames(design))
#fit2 <- contrasts.fit(fit, contrast)
#fit2 <- eBayes(fit2)
#topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#cont.dif <- makeContrasts(
#  DifII_I =(PH.II-PT.II)-(PH.I-PT.I),
#  DifIII_II=(PH.III-PT.III)-(PH.II-PT.II),
#  DiffIII_I=(PH.III-PT.III)-(PH.I-PT.I),
#  levels=design)
#  fit2 <- contrasts.fit(fit, cont.dif)
#  fit2 <- eBayes(fit2)
#topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

### Age >40, age<40
#annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Age_of_disease_onset")]
annotation.h=annotationh[rownames(annotationh)%in%colnames(data)[25:56],colnames(annotationh)%in%c("Gender","Age_of_disease_onset")]

annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.h[annotation.h[,2]>=40, 2]<-"old"
annotation.h[annotation.h[,2]< 40,2]<-"young"
#annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Age_of_disease_onset")]
annotation.t=annotationt[rownames(annotationt)%in%colnames(data)[57:87],colnames(annotationt)%in%c("Gender","Age_of_disease_onset")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.t[annotation.t[,2]>=40, 2]<-"old"
annotation.t[annotation.t[,2]<40, 2]<-"young"
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

targets=data.frame(annotation.pso.h.t[,c(1,3)])
pheno=paste(targets$Sample, targets$Age_of_disease_onset, sep=".")
pheno=factor(pheno, levels=unique(pheno))
design <- model.matrix(~0 +pheno)
colnames(design)<-levels(pheno)
#colnames(design)=str_replace(colnames(design),"pheno","")
fit <- lmFit(p, design)

#contrasts of interest
 cont.matrix <- makeContrasts(
   OldvsYounginPT=PT.old-PT.young,
   OldvsYounginPH=PH.old-PH.young,
   Diff=(PH.old-PH.young)-(PT.old-PT.young),
   levels=design)
 fit2 <- contrasts.fit(fit, cont.matrix)
 fit2 <- eBayes(fit2)

adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt1=toptable(fit2,coef=1,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt2=toptable(fit2,coef=2,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt3=toptable(fit2,coef=3,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)


###Results
#> tt1
#logFC        t     P.Value  adj.P.Val         B
#CTLA4 1.63004 3.227419 0.002143678 0.03751437 -1.512352
#> tt2
#logFC        t      P.Value    adj.P.Val          B
#NLRP1 1.180474 4.935930 7.075138e-06 0.0002476298  3.5341099
#IFIH1 1.497466 4.100469 1.300241e-04 0.0015169473  0.7648948
#CTLA4 1.446181 3.177619 2.477128e-03 0.0111941166 -1.9672288
#FOXP3 1.012937 3.145971 2.558655e-03 0.0111941166 -2.0290130
#KLRK1 1.128399 2.964326 4.536615e-03 0.0158781516 -2.5145683
#AIM2  1.113410 2.806603 7.315827e-03 0.0232776299 -2.8931975
#> tt3
#data frame with 0 columns and 0 rows


#design <- model.matrix(~SampleType*AgeOfOnset)
#fit <- lmFit(p, design)
#temporarily swap PH i PT
targets2=rbind(targets[1:32,], targets[33:63,])
SampleType=factor(targets2$Sample, levels=c("PT","PH"))
AgeOfOnset=factor(targets2$Age_of_disease_onset, levels=unique(targets2$Age_of_disease_onset))
contrasts(SampleType) =contr.sum(2)
contrasts(AgeOfOnset) = contr.sum(2)
design = model.matrix(~Sample+AgeOfOnset)
fit =lmFit(p, design)
cont.matrix = cbind(OldvsYounginPT=c(0,0,-2,-2),OldvsYounginPH=c(0,0,-2,2),Diff=c(0,0,0,4))
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
adj.method="fdr" 
thr.p.value=0.05
thr.lfc=1
tt11=toptable(fit2,coef=1,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt21=toptable(fit2,coef=2,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt31=toptable(fit2,coef=3,p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
#the results are identical with the first approach
#> tt11
#logFC        t     P.Value  adj.P.Val         B
#CTLA4 1.63004 3.320099 0.001632745 0.02857303 -1.281748
#> tt21
#logFC        t      P.Value   adj.P.Val          B
#NLRP1 1.180474 4.973978 6.167086e-06 0.000215848  3.6611535
#IFIH1 1.497466 4.008937 1.763370e-04 0.002057265  0.4928487
#CTLA4 1.446181 3.075089 3.323016e-03 0.014140848 -2.2108772
#FOXP3 1.012937 3.062472 3.263384e-03 0.014140848 -2.2302183
#KLRK1 1.128399 2.703659 9.195786e-03 0.032185251 -3.0778234
#AIM2  1.113410 2.590086 1.280639e-02 0.040747595 -3.3299323
#> tt31
#data frame with 0 columns and 0 rows



######Psoriasis in the family

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Psoriasis_in_family")]

annotation.h[,2]=str_replace(annotation.h[,2],"au.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"si.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"fa.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"mo.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"so.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"bro.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"chi.*","yes")
annotation.h[,2]=str_replace(annotation.h[,2],"2_.*","yes")
annotation.h=annotation.h[order(annotation.h[,2]),]

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Psoriasis_in_family")]
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

targets=data.frame(annotation.pso.h.t[,c(1,3)])
pheno=paste(targets$Sample, targets$Psoriasis_in_family, sep=".")
#pheno=factor(pheno, levels=pheno)
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")
fit <- lmFit(p, design)
#1.No psoriasis in the family
contrast = makeContrasts(PHvsPTinNO=PH.no-PT.no, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#results <- decideTests(fit2)
#vennDiagram(results)
#2. Patients with psoritic family history
contrast = makeContrasts(PHvsPTinYES=PH.yes-PT.yes, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#3. whhich genes respond diﬀerently in patients with different phenotype in healthy and diseased condititons.
cont.dif <- makeContrasts(
  Difyes_no =(PH.yes-PT.yes)-(PH.no-PT.no),
  levels=design)
fit2 = contrasts.fit(fit, cont.dif)
fit2 = eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)



##Psoriatic artritis

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Psoriatic_arthritis")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Psoriatic_arthritis")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

targets=data.frame(annotation.pso.h.t[,c(1,3)])
pheno=paste(targets$Sample, targets$Psoriatic_arthritis, sep=".")
#pheno=factor(pheno, levels=pheno)
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")
fit <- lmFit(p, design)
#1.No psoriatic arthritis in the family
contrast = makeContrasts(PHvsPTinNO=PH.no-PT.no, levels=colnames(design))
fit2 = contrasts.fit(fit, contrast)
fit2 = eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#results <- decideTests(fit2)
#vennDiagram(results)
#2. Patients with psoritic family history
contrast = makeContrasts(PHvsPTinYES=PH.yes-PT.yes, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#3. whhich genes respond diﬀerently in patients with different phenotype in healthy and diseased condititons.
cont.dif <- makeContrasts(
  Difyes_no =(PH.yes-PT.yes)-(PH.no-PT.no),
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)



### Nail_involvment
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

targets=data.frame(annotation.pso.h.t[,c(1,3)])
pheno=paste(targets$Sample, targets$Nail_involvment, sep=".")
#pheno=factor(pheno, levels=pheno)
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")
fit <- lmFit(p, design)
#1.No nail involvment
contrast = makeContrasts(PHvsPTinNO=PH.no-PT.no, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#results <- decideTests(fit2)
#vennDiagram(results)
#2. Nail involvment
contrast = makeContrasts(PHvsPTinYES=PH.yes-PT.yes, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#3. whhich genes respond diﬀerently in patients with different phenotype in healthy and diseased condititons.
cont.dif <- makeContrasts(
  Difyes_no =(PH.yes-PT.yes)-(PH.no-PT.no),
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)



#PASI 

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Gender","PASI_(activity_score)")]
annotation.h[,2]=str_replace(annotation.h[,2],",",".")
annotation.h=data.frame(annotation.h)
annotation.h[,2]=as.numeric(as.vector(annotation.h[,2]))
annotation.h[,2]=cut(annotation.h[,2], c(0, 10, 20, Inf))

annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Gender","PASI_(activity_score)")]
annotation.t[,2]=str_replace(annotation.t[,2],",",".")
annotation.t=data.frame(annotation.t)
annotation.t[,2]=as.numeric(as.vector(annotation.t[,2]))
annotation.t[,2]=cut(annotation.t[,2], c(0, 10, 20, Inf))

annotation.pso.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.pso.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.pso.h.t=cbind(Sample,annotation.pso.h.t)

targets=data.frame(annotation.pso.h.t[,c(1,3)])
pheno=paste(targets$Sample, targets$PASI_.activity_score., sep=".")
#pheno=factor(pheno, levels=pheno)
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")
colnames(design)=c("PH.10","PH.10.20","PH.20","PT.10","PT.10.20","PT.20")
fit <- lmFit(p, design)

#1. Type I
contrast = makeContrasts(PHvsPTin010=PH.10-PT.10, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000 )
#tt = topTable(fit2, coef = 1, p.value = 0.05, lfc = 1, adjust.method ="fdr", number = 10000000) ##no diff expressed genes in healthy patients with different phototypes
#results <- decideTests(fit2)
#vennDiagram(results)
#2. type II
contrast = makeContrasts(PHvsPTin1020=PH.10.20-PT.10.20, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr",p.value = 0.05, number = 10000000 )
#tt = topTable(fit2, coef = 1, p.value = 0.05, lfc = 1, adjust.method ="fdr", number = 10000000) ##
#3 Type III
contrast = makeContrasts(PHvsPTin20=PH.20-PT.20, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr",p.value = 0.05, number = 10000000 )
#4. Diff of differences
cont.dif <- makeContrasts(
  Dif1020_10 =(PH.10.20-PT.10.20)-(PH.10-PT.10),
  Dif20_1020=(PH.20-PT.20)-(PH.10.20-PT.10.20),
  Diff20_10=(PH.20-PT.20)-(PH.10-PT.10),
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)


### Age >40, age<40 and psoriasis  in the family
annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
annotation.h=annotation.h[order(annotation.h[,1]),]
annotation.h[annotation.h[,1]>=40, 1]<-"old"
annotation.h[annotation.h[,1]< 40,1]<-"young"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family")]
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

targets=data.frame(annotation.pso.h.t)
pheno=paste(targets$Sample,targets$Psoriasis_in_family,  targets$Age_of_disease_onset, sep=".")
#pheno=factor(pheno, levels=pheno)
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")

fit <- lmFit(p, design)
#1. no pso in the family Young
contrast = makeContrasts(PHvsPTin.no.young=PH.no.young-PT.no.young, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#2. no pso in the family Old
contrast = makeContrasts(PHvsPTin.no.old=PH.no.old-PT.no.old, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#3. yes pso in the family Young
contrast = makeContrasts(PHvsPTin.yes.young=PH.yes.young-PT.yes.young, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#results <- decideTests(fit2)
#vennDiagram(results)

#4. yes pso in the family Old
contrast = makeContrasts(PHvsPTin.yes.old=PH.yes.old-PT.yes.old, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)


#5. whhich genes respond diﬀerently in patients with different phenotype in healthy and diseased condititons. No psoriasis in the family
cont.dif <- makeContrasts(
  Dif.no.old_young =(PH.no.old-PT.no.old)-(PH.no.young-PT.no.young),
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#tt = topTable(fit2, coef = 1, p.value = 0.05, lfc = 1, adjust.method ="fdr", number = 10000000)####no diff expressed genes

#6. whhich genes respond diﬀerently in patients with different phenotype in healthy and diseased condititons.Pso in the family
cont.dif <- makeContrasts(
  Dif.yes.old_young =(PH.yes.old-PT.yes.old)-(PH.yes.young-PT.yes.young),
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
#tt = topTable(fit2, coef = 1, p.value = 0.05, lfc = 1, adjust.method ="fdr", number = 10000000)####no diff expressed genes


#### Age of disease onset <40 yr, >= 40 yr, PASI 

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:34],colnames(annotationh)%in%c("PASI_(activity_score)","Age_of_disease_onset")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.h[annotation.h[,2]>=40, 2]<-"old"
annotation.h[annotation.h[,2]< 40,2]<-"young"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[35:69],colnames(annotationt)%in%c("Age_of_disease_onset", "PASI_(activity_score)")]
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


targets=data.frame(annotation.pso.h.t)
colnames(targets)[2]=c("PASI")


pheno=paste(targets$Sample, targets$PASI,targets$Age_of_disease_onset, sep=".")
design <- model.matrix(~0 +pheno)
colnames(design)
colnames(design)=str_replace(colnames(design),"pheno","")

fit <- lmFit(p, design)

#1. Type PASI 0-10, early diasease onset
contrast = makeContrasts(PHvsPTin010young=PH.0_10.young-PT.0_10.young, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#2. Type PASI 0-10, late diasease onset
contrast = makeContrasts(PHvsPTin010old=PH.0_10.old-PT.0_10.old, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)


#3. Type PASI 10-20, early diasease onset
contrast = makeContrasts(PHvsPTin1020young=PH.10_20.young-PT.10_20.young, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#4. Type PASI 10-20, late diasease onset
contrast = makeContrasts(PHvsPTin1020old=PH.10_20.old-PT.10_20.old, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#5. Type PASI 20-inf, early diasease onset
contrast = makeContrasts(PHvsPTin20Infyoung=PH.20_Inf.young-PT.20_Inf.young, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#6. Type PASI 20-inf, late diasease onset
contrast = makeContrasts(PHvsPTin20Infold=PH.20_Inf.old-PT.20_Inf.old, levels=colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)

#####to be finished!!!!!
cont.dif <- makeContrasts(
  Dif.0_10.old.young =(PH.0_10.old-PT.0_10.old)-(PH.0_10.young-PT.0_10.young),
  Dif.10_20.old.young=(PH.10_20.old-PT.10_20.old)-(PH.10_20.young-PT.10_20.young),
  Diff.20_Inf.old.young=(PH.20_Inf.old-PT.20_Inf.old)-(PH.20_Inf.young-PT.20_Inf.young),
  
  levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="fdr", p.value = 0.05, number = 10000000)
