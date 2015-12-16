library("limma")
library(stringr)
library(pheatmap)
library(reshape2)
setwd("~/Documents/KAI/cytoskin/BLOOD_SKIN_COMP/")
###########################Put Controlls, Pso, Vit in one matrix##########
#load(file = "~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data_log2_filt.RData")
#data=data_log2_filt
#pheatmap(data, cluster_rows = F, cluster_cols = F, scale = "none", main="All data after standardization")
#controlls
DataControllsSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.control_newfilt.txt", sep="\t", header=T)
###Psofiasis#
DataPsoriasisSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.psoriasis_newfilt.txt", sep="\t", header=T )      
###Vitiliigo#
DataVitiliigoSkin=read.table(file="~/Documents/KAI/cytoskin/after_filt/17112014/removedC024/data.vitiliigo_newfilt.txt", sep="\t", header=T )      
ctrls=DataControllsSkin
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:43]
ctrls <- apply(ctrls[,2:43],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisSkin)# psoriasis
rows=rownames(pso[2:43,])
cols=as.character(pso[1,])
ps <- apply(pso[2:43,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiliigoSkin)# vitiliigo
rows=rownames(vit[2:43,])
cols=as.character(vit[1,])
vi <- apply(vit[2:43,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

data.all=cbind(ctrls, ps,vi)
summary(data.all)

data.log=log2(data.all)
data.skin=data.log

colnames(data.skin)=str_replace(colnames(data.skin),"PH","PL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"PT","PNL")#healthy
colnames(data.skin)=str_replace(colnames(data.skin),"VH","VL")#sick
colnames(data.skin)=str_replace(colnames(data.skin),"VT","VNL")#healthy
rownames(data.skin)=gsub("IL8b","CXCL8", rownames(data.skin))
save(data.skin, file="data.skin.Rdata")

pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL, PNL, VL and VNL groups",
         filename="data_log_NAs.png")

##take only PSL samples

psl=data.skin[,25:58]
pheatmap(psl, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in psoriatic lesional skin",
         filename="PSL_log_NAs.png")




#For interleukins Il17A IL17F skin blood  draw and calculate the values for correlation profiles in PL and PMA, SEB, unstimulated

psl_il=psl[rownames(psl)%in%c("IL17A", "IL17F"),]#data for correlation comparison skin

####Load blood sample' data
#controlls
load(file="~/Documents/KAI/cytoskin/BLOOD/datacontrol_QR.RData")
###Psofiasis#
load(file="~/Documents/KAI/cytoskin/BLOOD/datapsoriasis_QR.RData")     
###Vitiliigo#
load(file="~/Documents/KAI/cytoskin/BLOOD/datavitiligo_QR.RData")
###PCA plot of the data

ctrls=DataControllsBlood_QR
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:16]
ctrls <- apply(ctrls[,2:16],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DataPsoriasisBlood_QR)# psoriasis
rows=rownames(pso[2:16,])
cols=as.character(pso[1,])
ps <- apply(pso[2:16,],2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#VIT.
vit=t(DataVitiligoBlood_QR)# vitiliigo
rows=rownames(vit[2:16,])
cols=as.character(vit[1,])
vi <- apply(vit[2:16,],2,as.numeric)
colnames(vi)=cols
rownames(vi)=rows

#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[16,]=gsub("C.*_","",ctrlstst[16,])
ctrlstst=ctrlstst[,order(ctrlstst[16,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[17,]=gsub("C","",ctrlstst[17,])
ctrlstst[17,]=gsub("_.*","",ctrlstst[17,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,17]=as.numeric(as.character(ctrlstst[,17]))
C_72h=ctrlstst[ctrlstst[,16]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,17]),]
C_PMA=ctrlstst[ctrlstst[,16]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,17]),]
C_SEB=ctrlstst[ctrlstst[,16]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,17]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:15,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
ctrls.log=log2(ctrls.data)
rownames(ctrls.log)=genes
rownames(ctrls.log)
pheatmap(ctrls.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in controls (72h, PMA, SEB)")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[16,]=gsub("P.*_","",pstst[16,])
pstst=pstst[,order(pstst[16,])]
pstst=rbind(pstst, colnames(pstst))
pstst[17,]=gsub("P","",pstst[17,])
pstst[17,]=gsub("_.*","",pstst[17,])
pstst=as.data.frame(t(pstst))
pstst[,17]=as.numeric(as.character(pstst[,17]))
P_72h=pstst[pstst[,16]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,17]),]
P_PMA=pstst[pstst[,16]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,17]),]
P_SEB=pstst[pstst[,16]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,17]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:15,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))
ps.log=log2(ps.data)
rownames(ps.log)=genes
rownames(ps.log)
pheatmap(ps.log, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in blood samples of psoriasis patients (72h, PMA, SEB)")
ps.blood.il=ps.log[rownames(ps.log)%in%c("IL17A", "IL17F"),]#data cor correlation comparson R


#plot the data in PL and blood  72h, seb, pma
skin_blood=cbind(psl[rownames(psl)%in%rownames(ps.log),], ps.log[rownames(ps.log)%in%rownames(psl),])
pheatmap(skin_blood, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in blood samples of psoriasis patients (PL,72h, PMA, SEB)")


#######Correlation of the interleukins' profiles IL17A, IL17F in PL and SEB, PMA, 72h
head(psl_il)
head(ps.blood.il)
ps.blood.il.72h=ps.blood.il[,1:35]
ps.blood.il.pma=ps.blood.il[,36:70]
ps.blood.il.seb=ps.blood.il[,71:105]
colnames(ps.blood.il.seb)=gsub("_SEB","",colnames(ps.blood.il.seb))
colnames(ps.blood.il.seb)=gsub("0","",colnames(ps.blood.il.seb))
colnames(ps.blood.il.72h)=gsub("_72h","",colnames(ps.blood.il.72h))
colnames(ps.blood.il.72h)=gsub("0","",colnames(ps.blood.il.72h))
colnames(ps.blood.il.pma)=gsub("_PMA","",colnames(ps.blood.il.pma))
colnames(ps.blood.il.pma)=gsub("0","",colnames(ps.blood.il.pma))
colnames(psl_il)=gsub("PL","P",colnames(psl_il))
colnames(psl_il)=gsub("0","",colnames(psl_il))
length(ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17A"), ])#35
length(psl_il[rownames(psl_il)%in%c("IL17A"),])#34

#http://www.r-bloggers.com/more-on-exploring-correlations-in-r/
IL17A=rbind(psl_il[rownames(psl_il)%in%c("IL17A"),],
            ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL17A)=c("PL", "SEB","PMA", "72h")
cor(t(IL17A),method="spearman",use="pairwise.complete.obs")
#      PL        SEB        PMA      72h
#PL   1.0000000 -0.3218391 -0.3264957 0.5
#SEB -0.3218391  1.0000000  0.1052419 0.2
#PMA -0.3264957  0.1052419  1.0000000 0.5
#72h  0.5000000  0.2000000  0.5000000 1.0
cor(t(IL17A),method="pearson",use="pairwise.complete.obs")
#PL         SEB         PMA        72h
#PL   1.0000000 -0.25264095 -0.39970445 0.76386257
#SEB -0.2526409  1.00000000  0.03239979 0.06812561
#PMA -0.3997044  0.03239979  1.00000000 0.18859489
#72h  0.7638626  0.06812561  0.18859489 1.00000000

cor.test.spearman <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="spearman")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
cor.test.pearson <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="pearson")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}


#pearson
cor.test.pearson(t(IL17A))
#    PL       SEB        PMA       72h
#PL  0.00000000 0.1946100 0.04306423 0.4466010
#SEB 0.19461001 0.0000000 0.86263066 0.9318744
#PMA 0.04306423 0.8626307 0.00000000 0.8792134
#72h 0.44660104 0.9318744 0.87921339 0.0000000
pt=melt(cor.test.pearson(t(IL17A)))
pcor=melt(cor(t(IL17A),method="pearson",use="pairwise.complete.obs"))
ptcor=cbind(pcor, pt[,3])
colnames(ptcor)=c("group1", "group2", "cor", "pval")
ptcor=ptcor[-c(5,9,10,13,14,15),]

#spearman
cor.test.spearman(t(IL17A))
#    PL        SEB       PMA        72h
#PL  0.00000000 0.09522312 0.1038560 1.00000000
#SEB 0.09522312 0.00000000 0.5717242 0.91666667
#PMA 0.10385595 0.57172422 0.0000000 1.00000000
#72h 1.00000000 0.91666667 1.0000000 0.08333333
st=melt(cor.test.spearman(t(IL17A)))
scor=melt(cor(t(IL17A),method="spearman",use="pairwise.complete.obs"))
stcor=cbind(scor, st[,3])
colnames(stcor)=c("group1", "group2", "cor", "pval")
stcor=stcor[-c(5,9,10,13,14,15),]


chart.Correlation <-
  function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...)
  { # @author R Development Core Team
    # @author modified by Peter Carl
    # Visualization of a Correlation Matrix. On top the (absolute) value of the
    # correlation plus the result of the cor.test as stars. On botttom, the
    # bivariate scatterplots, with a fitted line
    
    x = checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      
      test <- cor.test(x,y, method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    f <- function(t) {
      dnorm(t, mean=mean(x), sd=sd.xts(x) )
    }
    hist.panel = function (x, ...) {
      par(new = TRUE)
      hist(x,
           col = "light gray",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD")
      lines(density(x, na.rm=TRUE),
            col = "red",
            lwd = 1)
      #lines(f, col="blue", lwd=1, lty=1) how to add gaussian normal overlay?
      rug(x)
    }
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=hist.panel, method=method, ...)
    else
      pairs(x, gap=0, lower.panel=panel.smooth, upper.panel=panel.cor, method=method, ...) 
  }


scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
# plot the data
###shall i csale the rows patient wise?
library(PerformanceAnalytics)
chart.Correlation((t(IL17A)),method="spearman", main="Spearman correlation between expression profiles of IL17A PL, SEB, PMA, 72h ")#check this one
chart.Correlation(t((IL17A)),method="pearson", main="Pearson correlation between expression profiles of IL17A PL, SEB, PMA, 72h ")#check this one

summary(scale_rows(t(IL17A)))
#pheatmap(IL17A, cluster_rows = F, cluster_cols = F,scale = "none",clustering_distance_rows="spearman", main="correlation")
p#lot(ptest[1,], ptest[2,], main="Scatterplot IL17A, psoriasis blood SEB stimulated, psoriasis lesional skin")
#abline(lm(ptest[2,]~ptest[1,]), col="red")

######IL17F correlations
IL17F=rbind(psl_il[rownames(psl_il)%in%c("IL17F"),],
            ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL17F)=c("PL", "SEB","PMA", "72h")

cor(t(IL17F),method="spearman",use="pairwise.complete.obs")
#PL         SEB         PMA         72h
#PL   1.00000000  0.01135191 -0.19285714  1.00000000
#SEB  0.01135191  1.00000000  0.36343733 -0.08571429
#PMA -0.19285714  0.36343733  1.00000000  0.08571429
#72h  1.00000000 -0.08571429  0.08571429  1.00000000

cor(t(IL17F),method="pearson",use="pairwise.complete.obs")
#PL         SEB         PMA        72h
#PL   1.00000000 -0.07853554 -0.27037733 1.00000000
#SEB -0.07853554  1.00000000  0.33756829 0.14505312
#PMA -0.27037733  0.33756829  1.00000000 0.09676006
#72h  1.00000000  0.14505312  0.09676006 1.00000000

#pearson
cor.test.pearson(t(IL17F))#not enought finit observations
cor(t(IL17F[1:3,]),method="pearson",use="pairwise.complete.obs")
#PL         SEB        PMA
#PL   1.00000000 -0.07853554 -0.2703773
#SEB -0.07853554  1.00000000  0.3375683
#PMA -0.27037733  0.33756829  1.0000000
pt=melt(cor.test.pearson(t(IL17F[1:3,])))
pcor=melt(cor(t(IL17F[1:3,]),method="pearson",use="pairwise.complete.obs"))
ptcor=cbind(pcor, pt[,3])
colnames(ptcor)=c("group1", "group2", "cor", "pval")
ptcor=ptcor[-c(4,7,8),]

#spearman
cor.test.spearman(t(IL17F))

st=melt(cor.test.spearman(t(IL17F)))
scor=melt(cor(t(IL17F),method="spearman",use="pairwise.complete.obs"))
stcor=cbind(scor, st[,3])
colnames(stcor)=c("group1", "group2", "cor", "pval")
stcor=stcor[-c(5,9,10,13,14,15),]
# plot the data

chart.Correlation(t(IL17F),method="spearman", main="Spearman correlation between expression profiles of IL17F PL, SEB, PMA, 72h ")#check this one
chart.Correlation(t(IL17F[1:3,]),method="pearson", main="Pearson correlation between expresson profiles of IL17F PL, SEB, PMA, 72h ")#check this one

