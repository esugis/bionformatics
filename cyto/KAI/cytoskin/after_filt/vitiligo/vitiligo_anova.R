
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

#####################################################################################################################################
############################################ Select diffexp genes in log-transformed data for further analysis  ############################
##Since the changes can be already in non-leasional skin, we first select the genes that are diff expressed in (VH+VT)vsK

#test on log-transformed data without conversion to Z-scores
load("~/Documents/KAI/cytoskin/after_filt/17112014/noscale/data_log2_filt.RData")
data_log=data_log2_filt[,c(1:24,88:118)]
pheatmap(data_log, cluster_rows = F, cluster_cols = F, scale = "none", main="Gene expression in control and vitiligo group")
summary(data_log)
data=data_log
v=data_log2_filt[,88:118]
data=as.data.frame(data)
dt=t(data)
dt=dt[25:55,]
dt=as.data.frame(dt)
dt=cbind(dt, rownames(dt))
colnames(dt)[36]="sample"
dt[,36]=as.character(dt[,36])




##diff expressed PH+PT vs C
t=data
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("VH","VT")
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
dim(tt_log)# 9 genes
#logFC         t      P.Value    adj.P.Val          B
#CTLA4   -2.3329434 -7.179217 4.845044e-09 1.550414e-07 10.4882717
#WIPI1    0.5829671  4.173715 1.033530e-04 1.285631e-03  0.5926660
#CXCL10  -2.3084534 -4.102455 1.607039e-04 1.285631e-03  0.3042263
#IFNG    -2.3084534 -4.102455 1.607039e-04 1.285631e-03  0.3042263
#IFIH1   -0.7069381 -3.623446 6.344414e-04 4.060425e-03 -1.1150293
#IL22RA1 -1.2042291 -3.353904 1.524216e-03 8.129151e-03 -1.8910440
#TNF     -0.7994326 -3.245193 1.964674e-03 8.981365e-03 -2.2009334
#KLRK1   -0.7040773 -2.755170 8.062006e-03 3.224802e-02 -3.4566030
#CCL5    -0.6901946 -2.641730 1.062270e-02 3.776960e-02 -3.7559446
selected_genes_log=rownames(tt_log)
selected_genes_log_all= rownames(tt_log_all)
save(tt_log,tt_log_all,selected_genes_log,selected_genes_log_all, file = "~/Documents/KAI/cytoskin/after_filt/vitiligo/diffexp_cvtvh_log.RData")




############## Read in metadata###########
v=read.table(file="~/Documents/KAI/cytoskin/vitiliigo_all.txt", sep="\t", header=T)
v=v[1:68,]
v.meta=v[,1:68]
head(v.meta)
colnames(v.meta)=as.character(unname(unlist(v.meta[1,])))
v.meta=v.meta[2:19,]
rownames(v.meta)=as.character(unname(unlist(v.meta[,1])))
v.meta=t(v.meta[,2:68])
rownames(v.meta)[33]=c("Kobner_fenomen")
annotation=t(v.meta)
annotationh=annotation
annotationt=annotation
rownames(annotationh)=str_replace(rownames(annotation),"V","VH")
rownames(annotationh)[2]="VH002"
rownames(annotationh)[8]="VH008"
rownames(annotationt)=str_replace(rownames(annotation),"V","VT")
rownames(annotationt)[2]="VT002"
rownames(annotationt)[8]="VT008"
###########################################
############################################# BSA <10%, 10-50%,>50% 
annotation.h=annotationh[rownames(annotationh)%in%colnames(v)[1:31],colnames(annotationh)%in%c("Gender","Involvment")]
annotation.h[,2]=str_replace(annotation.h[,2],"till 10%","10")
annotation.h[,2]=str_replace(annotation.h[,2],"11-50%","50")
annotation.h[,2]=str_replace(annotation.h[,2],"51-90%","90")
annotation.h[,2]=str_replace(annotation.h[,2],"51-90","90")
annotation.h=data.frame(annotation.h)
annotation.h[,2]=as.numeric(as.vector(annotation.h[,2]))
annotation.h[,2]=cut(annotation.h[,2], c(0, 10, 50, 90))

annotation.t=annotationt[rownames(annotationt)%in%colnames(v)[1:31],colnames(annotationt)%in%c("Gender","Involvment")]
annotation.t[,2]=str_replace(annotation.t[,2],"till 10%","10")
annotation.t[,2]=str_replace(annotation.t[,2],"11-50%","50")
annotation.t[,2]=str_replace(annotation.t[,2],"51-90%","90")
annotation.t[,2]=str_replace(annotation.t[,2],"51-90","90")
annotation.t=data.frame(annotation.t)
annotation.t[,2]=as.numeric(as.vector(annotation.t[,2]))
annotation.t[,2]=cut(annotation.t[,2], c(0, 10, 50, 90))

annotation.v.h.t=rbind(annotation.h, annotation.t)
Sample=rownames(annotation.v.h.t)
Sample=str_replace(Sample,"0.*","")
annotation.v.h.t=cbind(Sample,annotation.v.h.t)

annotation.v.h.t=as.data.frame(annotation.v.h.t)
annotation.v.h.t[,1]=as.character(annotation.v.h.t[,1])
annotation.v.h.t[,2]=as.character(annotation.v.h.t[,2])
annotation.v.h.t[,3]=as.character(annotation.v.h.t[,3])
colnames(annotation.v.h.t)[3]="Involvement"
annotation.v.h.t=cbind(annotation.v.h.t, interaction(annotation.v.h.t$Sample, annotation.v.h.t$Involvement))
colnames(annotation.v.h.t)[4]="interactions"
ann=cbind(annotation.v.h.t, rownames(annotation.v.h.t))
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
pso_durations_comp=outDf


save(pso_durations_comp, file = "~/Documents/KAI/cytoskin/after_filt/21112014/Psoriasis_durations_comp_group_comparison.RData")

##faset_wrap http://docs.ggplot2.org/0.9.3.1/facet_wrap.html

ggplot(outDf, aes(x=diff,y=groups,color=p_adj))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable)+
  # theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(size=6))+
  geom_errorbarh(aes(xmax = upr, xmin = lwr, height = .2))+
  ggtitle("Comparison of groups using log-transformed data.
          Groups are formed based on the combinations of factors (Skin:PH/PT, Duration_of_psoriasis: <5,5-15,>15 years")



#######################################################################################################
###################Test on nonlogtransformed data converted to Z scores ###############################

load("~/Documents/KAI/cytoskin/after_filt/21112014/data_all_filt_nontransformed.RData")

data_raw=data_all_filt[,1:87]
pheatmap(data_raw, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controlls and Psoriasis")

pheatmap(data_raw, cluster_rows = F, cluster_cols = F, scale = "row", main="Data Controlls and Psoriasis")


summary(data_raw)
data=scale_rows(data_raw)
summary(data)
pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="Data Controlls and Psoriasis")

##diff expressed PH+PT vs C
t=data
cols=colnames(t)
attribute=colnames(t)
attribute=str_replace(attribute,"0.*","")
coldata=colnames(t)
coldata=str_replace(coldata,"0.*","")
attribute = factor(attribute)
attribute2 = rep("Muu", length(attribute))
# Case-control variables and corresponding group values 
case.vals=c("PH","PT")
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
tt_z_lfc = toptable(fit, coef = 1, p.value = thr.p.value, lfc = thr.lfc, adjust.method = adj.method, number = 10000000)
tt_all_z=toptable(fit, coef = 1, p.value = thr.p.value, adjust.method = adj.method, number = 10000000)#diff expressed genes filtered only based on p-val <0.05
print(thr.p.value)
print(tt) 
#logFC         t      P.Value    adj.P.Val        B
#IL1RN  -1.025620 -4.481099 7.737715e-06 0.0001175657 3.324277
#IL1F6  -1.005809 -4.394541 1.153270e-05 0.0001175657 2.948863
#S100A9 -1.005204 -4.372169 1.277116e-05 0.0001175657 2.855508
#S100A8 -1.026526 -4.335299 1.509365e-05 0.0001175657 2.714685
#LCN2   -1.080236 -4.236083 2.351314e-05 0.0001175657 2.338201
#CXCL10 -1.012968 -3.942706 8.263202e-05 0.0003213468 1.176155
#IFNG   -1.012968 -3.942706 8.263202e-05 0.0003213468 1.176155
print(tt_all_z)#diff expressed genes filtered based on p-val <0.05. Fold change was not taken to account
#logFC         t      P.Value    adj.P.Val           B
#IL1RN   -1.0256203 -4.481099 7.737715e-06 0.0001175657  3.32427676
#IL1F6   -1.0058092 -4.394541 1.153270e-05 0.0001175657  2.94886329
#S100A9  -1.0052040 -4.372169 1.277116e-05 0.0001175657  2.85550793
#S100A8  -1.0265257 -4.335299 1.509365e-05 0.0001175657  2.71468508
#PI3     -0.9838566 -4.298627 1.779996e-05 0.0001175657  2.54142369
#LCN2    -1.0802363 -4.236083 2.351314e-05 0.0001175657  2.33820079
#DEFB1    0.9725343  4.249157 2.219066e-05 0.0001175657  2.33479497
#CXCL10  -1.0129684 -3.942706 8.263202e-05 0.0003213468  1.17615547
#IFNG    -1.0129684 -3.942706 8.263202e-05 0.0003213468  1.17615547
#TNF     -0.8415292 -3.668623 2.486131e-04 0.0008074323  0.09017637
#CCL2    -0.8384562 -3.663349 2.537644e-04 0.0008074323  0.06977511
#IL22RA1 -0.8313829 -3.624390 2.950362e-04 0.0008605224 -0.06744296
#IFIH1   -0.8153351 -3.476168 5.167267e-04 0.0013911874 -0.56605629
#CXCL2   -0.8733083 -3.329784 8.809105e-04 0.0020703652 -0.97565673
#CCL27    0.7833759  3.315309 9.276218e-04 0.0020703652 -1.09397572
#CCL5    -0.7575059 -3.309664 9.464526e-04 0.0020703652 -1.13537010
#CTLA4   -0.8662166 -3.264377 1.110857e-03 0.0022870577 -1.17720634
#CXCL1   -0.7927873 -2.863223 4.226159e-03 0.0082175306 -2.33852359
#PYCARD  -0.6205231 -2.698986 6.998787e-03 0.0128925027 -2.92474729
#FOXP3   -0.6091316 -2.572529 1.014926e-02 0.0177612049 -3.22614922
#IFNAR1   0.5625059  2.457678 1.404678e-02 0.0234112970 -3.53647177
#IFNGR   -0.5442206 -2.372514 1.773774e-02 0.0282191299 -3.73559242
#EOMES   -0.6803842 -2.245119 2.484178e-02 0.0347784924 -3.79408636
#IL20RA   0.5370517  2.324825 2.015572e-02 0.0306717420 -3.83901650
#MICB    -0.5451989 -2.292539 2.195177e-02 0.0320130028 -3.88702420
#IL8b    -0.7408070 -2.118071 3.426092e-02 0.0428261541 -3.95115561
#IL1b    -0.5516596 -2.191411 2.850795e-02 0.0383760837 -4.05966616
#AIM2    -0.5425715 -2.166420 3.036701e-02 0.0393646403 -4.11699407

selected_genes_z=rownames(tt_lfc)#7
selected_genes_all_z=rownames(tt_all_z)
length(selected_genes_all_z)#28

save(tt_z_lfc,tt_all_z,selected_genes_z,selected_genes_all_z, file = "~/Documents/KAI/cytoskin/after_filt/21112014/diffexp_zscores.RData")


#################################### Wilcoxon test Z-score data #########################################
p.val=c()

for(i in 1: nrow(data)){
  u=wilcox.test(data[i,25:87],data[i,1:24])$p.value  
  p.val=c(p.val,u)
}
p.val
p.val.adj=p.adjust(p.val, method="fdr")  
genes=cbind(rownames(data), p.val.adj)
genes=as.data.frame(genes)
genes[,1]=as.character(genes[,1])
genes[,2]=as.numeric(as.character(genes[,2]))
genes=genes[genes$p.val.adj< 0.05,]
dim(genes)
#V1    p.val.adj
#1     AIM2 2.058208e-04
#3     CCL2 8.501798e-05
#4    CCL20 3.339449e-05
#5    CCL27 2.527162e-03
#6     CCL5 1.764107e-04
#7    CTLA4 4.178859e-06
#8    CXCL1 1.414884e-04
#9   CXCL10 3.423536e-06
#10   CXCL2 6.346530e-04
#11   DEFB1 1.158874e-02
#12   EOMES 1.547919e-02
#13   FOXP3 6.875056e-03
#14   IFIH1 1.744757e-04
#16    IFNG 3.423536e-06
#17   IFNGR 4.811984e-02
#18    IL1b 3.695514e-02
#19   IL1F6 1.403962e-08
#20   IL1RN 2.305570e-07
#22 IL22RA1 8.761885e-05
#24    IL8b 2.656471e-04
#26    LCN2 5.883264e-06
#27    MICB 4.681199e-03
#29     PI3 3.709164e-05
#30  PYCARD 8.089010e-03
#31  S100A8 2.981038e-04
#32  S100A9 3.709164e-05
#33     TNF 7.980145e-04
selected_genes_wlx=genes[,1]


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


ggplot(tmp, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  geom_density(aes(group=interactions, colour=interactions,fill=interactions), alpha=0.5)+
  ggtitle("Density plot of gene expression in log-transformed data. Data distribution is shown in all coditions.")



#Plot distributions for Z-scored data in PH and Pt conditions
ggplot(tmp_z, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  # log from negative value scale_x_log10()+
  geom_density(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  ggtitle("Density plot of gene expression in data converted to Z-scores.")

ggplot(tmp_z, aes(x=expression))+
  facet_wrap(~variable)+
  theme_bw()+
  # log from negative value scale_x_log10()+
  geom_histogram(aes(group=Sample, colour=Sample,fill=Sample), alpha=0.5)+
  #aes(y = ..density..)+
  ggtitle("Histogram of gene expression in data converted to Z-scores.")

#Plot distributions for raw data in PH and Pt conditions
data_raw=data_all_filt[,1:87]
data_raw=as.data.frame(data_raw)
dt_raw=t(data_raw)
dt_raw=dt_raw[25:87,]
dt_raw=as.data.frame(dt_raw)
dt_raw=cbind(dt_raw, rownames(dt_raw))
colnames(dt_raw)[36]="sample"
dt_raw[,36]=as.character(dt_raw[,36])
dt.merged_raw=merge(ann,dt_raw,by="sample",all=T)
dt.merged.melt_raw=melt(dt.merged_raw)
colnames(dt.merged.melt_raw)[7]="expression"
tmp_raw <- subset(dt.merged.melt_raw,variable%in%selected_genes_wlx)


ggplot(tmp_raw, aes(x=expression))+
  facet_wrap(~variable)+ 
  geom_density(aes(group=interactions, colour=interactions,fill=interactions), alpha=0.5, na.rm=TRUE)+
  scale_x_log10()+
  theme_bw()+
  ggtitle("Density plot of gene expression in raw data. Data distribution is shown in all coditions.")

#plot(density(tmp_raw[tmp_raw$variable%in%c("CCL27"),7 ]))


