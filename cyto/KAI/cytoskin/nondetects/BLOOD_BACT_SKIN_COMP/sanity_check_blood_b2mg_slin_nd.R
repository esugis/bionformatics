load(file="data.skin.Rdata")
data.skin
load(file="~/Documents/KAI/cytoskin/BLOOD/All_data_log.RData")
data.blood=data.log
summary(data.skin)
summary(data.blood)
data=cbind(data.skin[rownames(data.skin)%in%rownames(data.blood),], data.blood[rownames(data.blood)%in%rownames(data.skin),])
head(data)
summary(data)
tdata=t(data)
data.na.omit=na.omit(tdata)
pca.all=prcomp(data.na.omit, scale=T, center=T)
summary(pca.all)
Samples=rownames(data.na.omit)
a=qplot(pca.all$x[,1],pca.all$x[,2],xlab="Principal Component 1",ylab= "Principal Component 2", label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))+theme_bw()
a1=a+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")

b=qplot(pca.all$x[,1],pca.all$x[,3],xlab="Principal Component 1",ylab= "Principal Component 3", label=Samples)+scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
b1=b+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")
c=qplot(pca.all$x[,2],pca.all$x[,3],xlab="Principal Component 2",ylab= "Principal Component 3", label=Samples) +scale_colour_manual(values=c("black","darkblue","seagreen3","red","darkorchid1"))
c1=c+geom_text(size=3,hjust=0.1,vjust=1.6,position = "identity") +theme_bw()+ theme(legend.position="none")

##install grid extra!!!
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(a)
grid.arrange(a1,b1,c1, legend, 
             ncol=2, nrow=2, widths=c(1/2,1/2))
