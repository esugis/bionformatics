library("limma")
library(stringr)
library(pheatmap)
library(reshape2)
setwd("~/Documents/KAI/cytoskin/nondetects/BLOOD_BACT_SKIN_COMP/")
###########################Put Controlls, Pso, Vit in one matrix##########
#load(file="~/Documents/KAI/cytoskin/nondetects/SKIN_ND/data_skin_cp_nd.RData")
data.skin=data_skin_cp_nd
pheatmap(data.skin, cluster_rows = F, cluster_cols = F,scale = "none", main="Gene expression in C, PL and PNL groups after imputation")

##take only PSL samples
psl=data.skin[,24:56]
pheatmap(psl, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in psoriatic lesional skin",
         filename="PSL_ND.png")

#For interleukins Il17A IL17F skin blood  draw and calculate the values for correlation profiles in PL and PMA, SEB, unstimulated
il=c("IL17A","IL17F", "IL22","IL4", "IL10" )
psl_il=psl[rownames(psl)%in%il,]#data for correlation comparison skin

####Load blood sample' data
load(file="~/Documents/KAI/cytoskin/nondetects/BLOOD_ND_bact/DDCTBlood_bact_F.RData")

ctrls=DDCTBlood_bact_F[1:72,]
rows=rownames(ctrls)
cols=colnames(ctrls)
ctrls <- apply(ctrls,2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO
pso=t(DDCTBlood_bact_F[72:177,])# psoriasis
rows=rownames(pso)
cols=colnames(pso)
ps <- apply(pso,2,as.numeric)
colnames(ps)=cols
rownames(ps)=rows


#split ctrls, ps and vi into groups 72h, PMA, SEB
#CTRL
ctrlstst=rbind(ctrls, colnames(ctrls))
ctrlstst[15,]=gsub("C.*_","",ctrlstst[15,])
ctrlstst=ctrlstst[,order(ctrlstst[15,])]
ctrlstst=rbind(ctrlstst, colnames(ctrlstst))
ctrlstst[16,]=gsub("C","",ctrlstst[16,])
ctrlstst[16,]=gsub("_.*","",ctrlstst[16,])
ctrlstst=as.data.frame(t(ctrlstst))
ctrlstst[,16]=as.numeric(as.character(ctrlstst[,16]))
C_72h=ctrlstst[ctrlstst[,15]%in%c("72h"),]
C_72h=C_72h[order(C_72h[,16]),]
C_PMA=ctrlstst[ctrlstst[,15]%in%c("PMA"),]
C_PMA=C_PMA[order(C_PMA[,16]),]
C_SEB=ctrlstst[ctrlstst[,15]%in%c("SEB"),]
C_SEB=C_SEB[order(C_SEB[,16]),]
ctrlstst_ordered=rbind(C_72h,C_PMA,C_SEB)
ctrls.data=t(ctrlstst_ordered)

ctrls.data=ctrls.data[1:14,]
genes=rownames(ctrls.data)
ctrls.data=as.matrix(apply(ctrls.data,2, as.numeric))
ctrls.log=log2(ctrls.data)
rownames(ctrls.log)=genes
rownames(ctrls.log)
pheatmap(ctrls.log, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in controls (72h, PMA, SEB). Blood B2mg after imputation")


#PSO
pstst=rbind(ps, colnames(ps))
pstst[15,]=gsub("P.*_","",pstst[15,])
pstst=pstst[,order(pstst[15,])]
pstst=rbind(pstst, colnames(pstst))
pstst[16,]=gsub("P","",pstst[16,])
pstst[16,]=gsub("_.*","",pstst[16,])
pstst=as.data.frame(t(pstst))
pstst[,16]=as.numeric(as.character(pstst[,16]))
P_72h=pstst[pstst[,15]%in%c("72h"),]
P_72h=P_72h[order(P_72h[,16]),]
P_PMA=pstst[pstst[,15]%in%c("PMA"),]
P_PMA=P_PMA[order(P_PMA[,16]),]
P_SEB=pstst[pstst[,15]%in%c("SEB"),]
P_SEB=P_SEB[order(P_SEB[,16]),]
pstst_ordered=rbind(P_72h,P_PMA,P_SEB)
ps.data=t(pstst_ordered)

ps.data=ps.data[1:14,]
genes=rownames(ps.data)
ps.data=as.matrix(apply(ps.data,2, as.numeric))

rownames(ps.data)=genes
rownames(ps.data)
pheatmap(ps.data, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in blood samples (B2mg) of psoriasis patients (72h, PMA, SEB) after imputation")
ps.blood.il=ps.data[rownames(ps.data)%in%il,]#data cor correlation comparson R


#plot the data in PL and blood  72h, seb, pma
skin_blood=cbind(psl[rownames(psl)%in%rownames(ps.data),], ps.data[rownames(ps.data)%in%rownames(psl),])
pheatmap(skin_blood, cluster_rows = F, cluster_cols = F,scale = "row", main="Gene expression in skin and blood samples (B2mg) of psoriasis patients (PL,72h, PMA, SEB) after imputation")


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
length(psl_il[rownames(psl_il)%in%c("IL17A"),])#33

#http://www.r-bloggers.com/more-on-exploring-correlations-in-r/
IL17A=rbind(psl_il[rownames(psl_il)%in%c("IL17A"),],
            ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL17A"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL17A)=c("PL", "SEB","PMA", "72h")
cor(t(IL17A),method="spearman",use="pairwise.complete.obs")
#             PL        SEB         PMA        72h
#PL   1.00000000  0.3308824 -0.07352941 -1.0000000
#SEB  0.33088235  1.0000000 -0.12338710 -0.2142857
#PMA -0.07352941 -0.1233871  1.00000000  0.7142857
#72h -1.00000000 -0.2142857  0.71428571  1.0000000

#cor(t(IL17A),method="pearson",use="pairwise.complete.obs")
#            PL         SEB         PMA        72h
#PL   1.0000000  0.25044693 -0.11869777 -1.0000000
#SEB  0.2504469  1.00000000 -0.09707931 -0.2887634
#PMA -0.1186978 -0.09707931  1.00000000  0.4173481
#72h -1.0000000 -0.28876341  0.41734811  1.0000000

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

#spearman
cor.test.spearman(t(IL17A))
#               PL       SEB        PMA          72h
#PL  1.115071e-05 0.1944042 0.78015452 1.000000e+00
#SEB 1.944042e-01 0.0000000 0.50687275 6.190972e-01
#PMA 7.801545e-01 0.5068728 0.00000000 8.809524e-02
#72h 1.000000e+00 0.6190972 0.08809524 4.960317e-05
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
chart.Correlation((t(IL17A)),method="spearman", main="Spearman correlation between expression profiles of IL17A in PL, SEB(Bact), PMA(Bact), 72h(Bact) ")#check this one
#chart.Correlation(t((IL17A)),method="pearson", main="Pearson correlation between expression profiles of IL17A in PL, SEB, PMA, 72h ")#check this one

#summary(scale_rows(t(IL17A)))
#pheatmap(IL17A, cluster_rows = F, cluster_cols = F,scale = "none",clustering_distance_rows="spearman", main="correlation")
#p#lot(ptest[1,], ptest[2,], main="Scatterplot IL17A, psoriasis blood SEB stimulated, psoriasis lesional skin")
#abline(lm(ptest[2,]~ptest[1,]), col="red")

######IL17F correlations
IL17F=rbind(psl_il[rownames(psl_il)%in%c("IL17F"),],
            ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL17F"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL17F)=c("PL","SEB","PMA", "72h")

cor(t(IL17F),method="spearman",use="pairwise.complete.obs")
#           PL       SEB       PMA       72h
#PL  1.0000000 0.3142857 0.4000000        NA
#SEB 0.3142857 1.0000000 0.2399774 0.1333333
#PMA 0.4000000 0.2399774 1.0000000 0.0500000
#72h        NA 0.1333333 0.0500000 1.0000000
#spearman
pval=cor.test.spearman(t(IL17F[1:3,])) #not enough finite observations
#PL       SEB          PMA
#PL  4.960317e-05 0.5638889 7.500000e-01
#SEB 5.638889e-01 0.0000000 2.807348e-01
#PMA 7.500000e-01 0.2807348 7.268764e-07 
padj=p.adjust(pval, method = "fdr")
#[1] 1.488095e-04 7.250000e-01 7.500000e-01 7.250000e-01 0.000000e+00 5.053226e-01 7.500000e-01
#[8] 5.053226e-01 3.270944e-06

chart.Correlation(t(IL17F[1:3,]),method="spearman", main="Spearman correlation between expression profiles of IL17F in PL, SEB(Bact), PMA(Bact)")#check this one

#########"IL22" correlation
IL22=rbind(psl_il[rownames(psl_il)%in%c("IL22"),],
            ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL22"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL22"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
            ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL22"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL22)=c("PL","SEB","PMA", "72h")

cor(t(IL22),method="spearman",use="pairwise.complete.obs")
#            PL        SEB        PMA       72h
#PL   1.0000000 -0.9428571 -0.9000000 0.4285714
#SEB -0.9428571  1.0000000  0.7312821 0.1828967
#PMA -0.9000000  0.7312821  1.0000000 0.0754386
#72h  0.4285714  0.1828967  0.0754386 1.0000000
chart.Correlation(t(IL22),method="spearman", main="Spearman correlation between expression profiles of IL22 in PL, SEB(Bact), PMA(Bact), 72h(Bact)")#check this one
#########IL4" is not present in psl_il
#########"IL10" correlation
IL10=rbind(psl_il[rownames(psl_il)%in%c("IL10"),],
           ps.blood.il.seb[rownames(ps.blood.il.seb)%in%c("IL10"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
           ps.blood.il.pma[rownames(ps.blood.il.pma)%in%c("IL10"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ],
           ps.blood.il.72h[rownames(ps.blood.il.72h)%in%c("IL10"),colnames(ps.blood.il.seb)%in%colnames(psl_il) ])
rownames(IL10)=c("PL","SEB","PMA", "72h")

cor(t(IL10),method="spearman",use="pairwise.complete.obs")
#              PL        SEB         PMA        72h
#PL   1.00000000  0.3350649 -0.01503759 -0.2706767
#SEB  0.33506494  1.0000000 -0.19072581 -0.1314516
#PMA -0.01503759 -0.1907258  1.00000000  0.4743842
#72h -0.27067669 -0.1314516  0.47438424  1.0000000

chart.Correlation(t(IL10),method="spearman", main="Spearman correlation between expression profiles of IL10 in PL, SEB(Bact), PMA(Bact), 72h(Bact)")#check this one
