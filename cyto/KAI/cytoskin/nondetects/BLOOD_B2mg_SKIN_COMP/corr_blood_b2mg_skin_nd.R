#calculate correlations on the scaled centered data by rows. Put together skin and blood and convert to z-scores
#then extract the row IL17F and compute the correlations of gene expression of IL17F in PL,SEB,PMA,72h
load(file="data.skin.Rdata")
data.skin
load(file="~/Documents/KAI/cytoskin/BLOOD/All_data_log.RData")
data.blood=data.log
summary(data.skin)
summary(data.blood)
data=cbind(data.skin[rownames(data.skin)%in%rownames(data.blood),], data.blood[rownames(data.blood)%in%rownames(data.skin),])
data.scale=scale(data)
data.IL=data[rownames(data)%in%c("IL17A","IL17F"),]
tdata.IL=t(data.IL)
IL17A=rbind(data.IL[rownames(data.IL)%in%c("IL17A"),25:58],
            data.IL[rownames(data.IL)%in%c("IL17A"),c(198:205,207:232 )],
            data.IL[rownames(data.IL)%in%c("IL17A"),c(233:240,242:267) ],
            data.IL[rownames(data.IL)%in%c("IL17A"),c(268:275,277:302) ])
rownames(IL17A)=c("PL", "72h","PMA", "SEB")

cor(t(IL17A),method="spearman",use="pairwise.complete.obs")
