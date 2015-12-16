setwd("~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/updated_names_PL_PLN/anova")
#load the gene expression  data 
load(file="~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/psoriasis_for_anova.RData")
load(file="~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/dt.RData")
load(file = "~/Documents/KAI/cytoskin/after_filt/updated_names_group_comparison/PSvsC_diffexp_log.RData")

annotation.h=annotationh[rownames(annotationh)%in%colnames(p)[1:63],colnames(annotationh)%in%c("Age_of_disease_onset", "Psoriasis_in_family","PASI_(activity_score)","Psoriatic_arthritis","Nail_involvment")]
annotation.h=annotation.h[order(annotation.h[,2]),]
annotation.h[annotation.h[,2]>=40, 2]<-"late"
annotation.h[annotation.h[,2]< 40,2]<-"early"
annotation.t=annotationt[rownames(annotationt)%in%colnames(p)[1:63],colnames(annotationt)%in%c("Age_of_disease_onset", "Psoriasis_in_family", "PASI_(activity_score)","Psoriatic_arthritis","Nail_involvment")]
annotation.t=annotation.t[order(annotation.t[,2]),]
annotation.t[annotation.t[,2]>=40, 2]<-"late"
annotation.t[annotation.t[,2]<40, 2]<-"early"


annotation.h[,5]=str_replace(annotation.h[,5],"au.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"si.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"fa.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"mo.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"so.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"bro.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"chi.*","yes")
annotation.h[,5]=str_replace(annotation.h[,5],"2_.*","yes")
annotation.h=annotation.h[order(annotation.h[,5]),]


annotation.t[,5]=str_replace(annotation.t[,5],"au.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"si.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"fa.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"mo.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"so.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"bro.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"chi.*","yes")
annotation.t[,5]=str_replace(annotation.t[,5],"2_.*","yes")
annotation.t=annotation.t[order(annotation.t[,5]),]

#######PASI cleaning ########
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
colnames(annotation.pso.h.t)[2]="PASI"

#############################
save(annotation.pso.h.t, file="annotation.pso.h.t.RData")
write.table(annotation.pso.h.t, file="annotation.pso.h.t.txt", quote=F, sep="\t", row.names=F)
