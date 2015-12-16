ctrls=DataControllsSkin
rows=as.character(ctrls[,1])
cols=colnames(ctrls)[2:44]
ctrls <- apply(ctrls[,2:44],2,as.numeric)
rownames(ctrls)=rows
ctrls=t(ctrls)

#PSO.h
pso.h=t(DataPsoriasisSkin[1:34,])#haige psoriasis
rows=rownames(pso.h[2:44,])
cols=as.character(pso.h[1,])
ph <- apply(pso.h[2:44,],2,as.numeric)
colnames(ph)=cols
rownames(ph)=rows

#PSO.t
pso.t=t(DataPsoriasisSkin[35:69,])#terve psoriasis  
rows=rownames(pso.t[2:44,])
cols=as.character(pso.t[1,])
pt <- apply(pso.t[2:44,],2,as.numeric)
colnames(pt)=cols
rownames(pt)=rows

#VIT.h
vit.h=t(DataVitiliigoSkin[1:16,])#haige vitiliigo
rows=rownames(vit.h[2:44,])
cols=as.character(vit.h[1,])
vh <- apply(vit.h[2:44,],2,as.numeric)
colnames(vh)=cols
rownames(vh)=rows

#VIT.t
vit.t=t(DataVitiliigoSkin[17:32,])#terve vitiliigo
rows=rownames(vit.t[2:44,])
cols=as.character(vit.t[1,])
vt <- apply(vit.t[2:44,],2,as.numeric)
colnames(vt)=cols
rownames(vt)=rows  





require(vioplot)
require(devtools)
require(digest)
source_gist("https://gist.github.com/mbjoseph/5852613")
vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE, side="both") 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q2 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  radj <- ifelse(side == "right", 0, 1)
  ladj <- ifelse(side == "left", 0, 1)
  if (!(is.null(h))) 
    args <- c(args, h = h)
  med.dens <- rep(NA, n)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q2[i] <- quantile(data, 0.5)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    med.dat <- do.call("sm.density", 
                       c(list(data, xlim=est.xlim,
                              eval.points=med[i], display = "none")))
    med.dens[i] <- med.dat$estimate
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    med.dens[i] <- med.dens[i] * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(x = c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              y = c(base[[i]], rev(base[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - radj*boxwidth/2, 
             q1[i], 
             at[i] + ladj*boxwidth/2, 
             q3[i], col = rectCol)
        # median line segment
        lines(x = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              y = rep(med[i],3))
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), 
              c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - radj*boxwidth/2, q3[i], at[i] + 
               ladj*boxwidth/2, col = rectCol)
        lines(y = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              x = rep(med[i],3))
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
plot(x=NULL, y=NULL,
     xlim = c(0.5, 2.5), ylim=c(min(values), max(values),
     type="n", ann=FALSE, axes=F)
axis(1, at=c(1, 2, 3),  labels=c("I", "II"))
axis(2)
for (i in unique(treatment)) {
  for (j in unique(group)){
    vioplot2(values[which(treatment == i & group == j)],
             at = ifelse(i =="I", 1, 2),
              
             side = ifelse(j == 1, "left", "right"),
             col = ifelse(j == 1, "purple", "lightblue"),
             add = T)
  }
}

title("Violin plot", xlab="Skin Phototype")
legend("bottomright", fill = c("purple", "lightblue"),
       legend = c("PT", "PH"), box.lty=0)


vioplot2(values,
         at = ifelse(i =="I", 1, 2),
         
         side = ifelse(j == 1, "left", "right"),
         col = ifelse(j == 1, "purple", "lightblue"),
         add = T)
}
pso=cbind(pt,ph)

pso=pso[ , ! apply( pso, 2 , function(x) all(is.na(x)) ) ]#removes columns with complete NAs
pso=pso[!rowSums(!is.finite(pso)),]#remove nas
annotation.pso.h.t=annotation.pso.h.t[order(annotation.pso.h.t[,3]),]

target=rownames( annotation.pso.h.t)#vector according to which we will order ph_pt_heatm
pso=pso[,match(target, colnames(pso))]#order columns by vector
pso=pso[ , ! apply( pso, 2 , function(x) all(is.na(x)) ) ]
one.t=pso[,2]
one.t=as.vector(one.t)
one.h=pso[,1]
one.h=as.vector(one.h)
two.h= pso[,3:16]
two.h=as.vector(two.h)
two.t= pso[,17:29] 
two.t=as.vector(two.t)
three.h =pso[,30:47]
three.h=as.vector(three.h)
three.t=pso[,48:length(colnames(pso))]
three.t=as.vector(three.t)

values=c(one.t, one.h, two.t, two.h, three.t, three.h)
treatment=c("I","II","III")
vioplot(one.t,one.h, two.t,two.h, three.t, three.h, names=c("I terve","I haige", "II terve","II haige", "III terve", "III haige"),col="purple", ylim=c(0,2000) )
title("Violin plots of Skin Phototypes")

require(beanplot)
treatment <- rep(c("I", "II","III"), each=n.each*2)#n.each is the length of the  expression vectors 
group <- rep(c(1, 2, 1, 2), each=n.each)
beanplot(values~group*treatment,mail="Bean plot", side="both", xlab="Phototype", col=list("purple",c("lightblue", "black")), axes=F)
axis(1, at=c(1, 2),  labels=c("I", "II", "III"))
axis(2)
legend("bottomright", fill = c("purple", "lightblue"),
       legend = c("PH", "PT"), box.lty=0)
