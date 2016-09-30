require("vegan")
require("MESS")
require("flowFDA")
require("gridExtra")
require("boot")
### Diversity from FCM data (Hill numbers D0, D1 and D2)
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### n = number of replicates
Diversity <- function(x, d=4, plot=FALSE, R=999){
  require("boot")
  D2.boot <- function(x,i) 1/sum((x[i]/sum(x[i]))^2)
  D1.boot <- function(x,i) exp(-sum((x[i]/sum(x[i]))*log(x[i]/sum(x[i]))))
  x <- x@basis/apply(x@basis, 1, max)
  D=matrix(ncol=3,nrow=nrow(x))
    ### Observed richness
    D0 = apply(x,1,FUN=function(x) {
      x = round(x,d);x<-x[x!=0];sum(x!=0)
    })
    ### D1
    D1 = apply(fbasis@basis,1,FUN=function(x) {
      x = round(x,d);x<-x[x!=0];
      boot(data=x,statistic=D1.boot,R=R)
    })
    ### D2
    D2 = apply(fbasis@basis,1,FUN=function(x) {
      x = round(x,d);x<-x[x!=0];
      boot(data=x,statistic=D2.boot,R=R)
    })
    results <- data.frame(Sample_name=attr(x,"dimnames")[[1]],
                          D0, 
                          t(data.frame(lapply(D1,FUN=function(x) c(mean(x$t),sd(x$t))))),
                          t(data.frame(lapply(D2,FUN=function(x) c(mean(x$t),sd(x$t))))))
    colnames(results) = c("Sample_name","D0","D1","sd.D1","D2",
                          "sd.D2")
    rownames(results) = attr(x,"dimnames")[[1]]
    if (plot==TRUE) {
      p <- ggplot(results, aes(x=seq(1:nrow(results)),y=D2)) + geom_point(shape=16,size=4,alpha=0.7,colour="blue")+
        geom_point(colour = "grey90", size = 1.5) + labs(x="Samples",y="Phenotypic diversity - D2")+
        geom_line(colour="blue",alpha=0.4,linetype=2)+
        geom_errorbar(aes(ymin=D2-sd.D2,ymax=D2+sd.D2), width=0.25)
      print(p)
    }
    cat(paste0("Alpha diversity metrics (D1,D2) have been computed after ",R," bootstraps"))
    return(results)
}

### Evenness based on cumulative paretocurves
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### The higher the evenness the higher the value
### n = number of replicates
Evenness <- function(x,d=3,n=1,plot=FALSE){
  x<- x@basis/apply(x@basis, 1, max)
  require(MESS)
  AUC = as.numeric(matrix(nrow= length(x[,1]), ncol=1))
  for(i in 1:nrow(x)){
    AUC[i] = auc(cum_Richness(x[i,],d=d)[,1],cum_Richness(x[i,],d=d)[,2])
    AUC[i] = 1-(AUC[i]-0.5)/(0.5)
  }
  if(n>1){
    results = matrix(nrow=length(x[,1])/3, ncol=2)
    results[,1]=trip(AUC,n)[,1]
    results[,2]=trip(AUC,n)[,2]
    results=data.frame(results)
    colnames(results)=c("Evenness","sdev")
    rownames(results) = trip_col(attr(x,"dimnames")[[1]],n)
  }
  else{
    results = matrix(nrow=length(x[,1]), ncol=1)
    results[,1] = AUC
    results = data.frame(results)
    colnames(results) = c("Evenness")
    rownames(results) = attr(x,"dimnames")[[1]]
  }
  if (plot==TRUE) {
    plot(results$Evenness,pch=21,bg=adjustcolor("blue",0.7),
         col=adjustcolor("blue",0.7),cex=1.5,las=1,ylab="Evenness",xlab="Samples")
  }
  return(results)
  print("1 = maximum evenness; 0 = minimum evenness")
}

# Cumulative richness for Pareto curve in Evenness function

cum_Richness <- function(x, d=3){
  x = round(x,d)
  x = x[x!=0]
  x = x[rev(order(x))]/sum(x)
  for(i in 1: (length(x)-1)){
    x[i+1] = x[i+1] + x[i] 
  }
  y=seq(0,1,by=(1/(length(x)-1)))
  result=cbind(y,x)
  return(result)
}

### Structural organisation (So), based on Koch et al. 2014
### So describes the relative difference in abundance of indivual bins, 
### and the average abundance of all distinct bins.
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### n = number of replicates
So <- function(x,d=3,n=1,plot=FALSE){
  x<- x@basis/apply(x@basis, 1, max)
  So = apply(x,1,FUN=function(x){
    x = round(x,d);x<-x[x!=0];sum(abs(x-mean(x)))/(length(x))
  })
  if(n>1){
    results = matrix(nrow=length(x[,1])/3, ncol=2)
    results[,1]=trip(So,n)[,1]
    results[,2]=trip(So,n)[,2]
    results=data.frame(results)
    colnames(results)=c("Organisation","sdev")
    rownames(results) = trip_col(attr(x,"dimnames")[[1]],n)
  }
  else{
    results <- So
    results = data.frame(results)
    colnames(results)=c("Organisation")
    rownames(results)=attr(x,"dimnames")[[1]]
  }
  if (plot==TRUE) {
    plot(Structural.organization.fbasis$Organisation,pch=21,
         bg=adjustcolor("blue",0.7),col=adjustcolor("blue",0.7),cex=1.5,
         las=1,ylab="SO",xlab="Samples")
  }
  return(results)
}

#### Calculating coefficient of variance (CV) of fingerprint
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### n = number of replicates
CV <- function(x,d=3,n=1,plot=FALSE){
  x <- x@basis/apply(x@basis, 1, max)
  CV= as.numeric(matrix(nrow=length(x[,1]), ncol=1))
  for(i in 1:length(x[,1])){
    CV[i] = 100*sd(round(x[i,],d)[round(x[i,],d)!=0])/mean(round(x[i,],d)[round(x[i,],d)!=0])
    # print(CV[i])
  }
  if(n>1){
    results = matrix(nrow=length(x[,1])/3, ncol=2)
    results[,1]=trip(CV,n)[,1]
    results[,2]=trip(CV,n)[,2]
    results=data.frame(results)
    colnames(results)=c("CV","sdev")
    rownames(results) = trip_col(attr(x,"dimnames")[[1]],n)
  }
  else{
    results <- CV
    results = data.frame(results)
    colnames(results)=c("CV")
    rownames(results)=attr(x,"dimnames")[[1]]
  }
  if (plot==TRUE) {
    plot(Coef.var.fbasis$CV,pch=21,bg=adjustcolor("blue",0.7),
         col=adjustcolor("blue",0.7),cex=1.5,las=1,ylab="CV",xlab="Samples")
  }
  return(results)
}

### Beta-diversity (ranked based) using Non-metric Multidimensional Scaling (NMDS)
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = rounding factor for densities 
### n = number of replicates
### dist = choice of distance metric 
### k = number of MDS dimensions 
### iter = number of random starts in search of stable solution
beta.div.fcm <- function(x, d=3, n=1, dist="bray",k=2,iter=100,ord.type=c("NMDS","PCoA")){
  x <- x@basis/apply(x@basis, 1, max)
  require('vegan')
  input <- matrix(nrow=nrow(x)/n,ncol=ncol(x))
  j<-1
  if(n>1){
    for(i in seq(1,nrow(x),n)){
      if(n>1)
        input[j,] <- round(colMeans(x[(i:(i+(n-1))),]),d)
      j=j+1
    }
    rownames(input) <- rownames(x)[seq(1,nrow(x),n)]
    input.dist <- vegdist(input,method=dist)
    if(ord.type=="NMDS") mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
    else mds.fbasis <- cmdscale(input.dist, k = 2, eig = TRUE, add = TRUE)
  }
  else{
    input.dist <- vegdist(x,method=dist)
    if(ord.type=="NMDS") mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
    else mds.fbasis <- cmdscale(input.dist, k = 2, eig = TRUE, add = TRUE)
  }
  return(mds.fbasis)
}

plot.beta.fcm <- function(x, color=NA,shape=NA,labels=c("Factor 1","Factor 2"),legend.pres=NULL){
  require('ggplot2')
  legend.ops<-NULL
  if (sum(is.na(color))>0) color=rep("f1",nrow(x$points))
  if (sum(is.na(shape))>0) {
    shape=rep("f2",nrow(x$points))
    legend.ops <- FALSE
    }
  var.pcoa <- eigenvals(x)/sum(eigenvals(x))
  PcoA <- as.data.frame(x$points)
  names(PcoA)[1:2] <- c('Axis1', 'Axis2')
  PcoA <- cbind(PcoA, color, shape)
  ggplot(PcoA, aes(x = Axis1, y = Axis2, color = color, shape= shape))+ 
    geom_point(alpha = 0.7, size = 4) +geom_point(colour = "grey90", size = 1.5)+
    scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                  "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
    labs(x=paste0("Axis1 (",round(100*var.pcoa[1],1),"%)"),y=paste0("Axis2 (",round(100*var.pcoa[2],1),"%)"))+
    ggtitle("Ordination of phenotypic fingerprints")+
    labs(color=labels[1],shape=labels[2])+
    guides(color=legend.pres, shape=legend.ops)
}

dist.fcm <- function(x, d=3, n=3, dist="bray"){
  x <- x@basis/apply(x@basis, 1, max)
  require('vegan')
  input <- matrix(nrow=nrow(x)/n,ncol=ncol(x))
  j<-1
  if(n>1){
    for(i in seq(1,nrow(x),n)){
      if(n>1)
        input[j,] <- round(colMeans(x[(i:(i+(n-1))),]),d)
      j=j+1
    }
    rownames(input) <- rownames(x)[seq(1,nrow(x),n)]
    input.dist <- vegdist(input,method=dist)
  }
  else{
    input.dist <- vegdist(x,method=dist)
  }
  return(input.dist)
}

###### This function subsets the input frames based on analysis time
###### Suited for on-line time series analysis
# flowData_transformed: preprocessed data (normalized, denoised)
# create: TRUE if you make want a new folder for the new flowframes
# analysis.length: dataframe with $time defining the total analysis time for each sample
# start: vector of length n(flowData_transformed) that indicates for each sample
# at what time point it should start the function 
# (e.g., first 10 minutes are irrelevant, start = 10*60s)
time.discretization <- function(flowData_transformed,analysis.length,create=FALSE,start=0,time.interval){
  for(j in 1:length(flowData_transformed)){
    number <- max(round((round(analysis.length/time.interval,0)+1)/10,0))
    old.wd <- getwd()
    if(create) {
      dir.create(paste(strsplit(rownames(analysis.length)[j],".fcs")[[1]], paste(time.interval), sep="_"))
      setwd(paste(strsplit(rownames(analysis.length)[j],".fcs")[[1]], paste(time.interval), sep="_"))
    }
    if(is.integer(analysis.length$time[j]/time.interval)) teller <- analysis.length$time[j]/time.interval
    else teller <- round(analysis.length$time[j]/time.interval,0)+1
    if(length(start)>0) teller = teller - start[j]/time.interval
    if(max(start) == 0) {
      start<-c()
      start[1:length(flowData_transformed)]<-0
    }
    res <- 0
    for(i in 1:teller){
      bottom <- (i-1)*time.interval + res + start[j]
      top <- i*time.interval + start[j]
      time.gate <- rectangleGate(filterId = "Time discretization", "Time" = c(bottom, top), "FL1-H" = c(0, 1))
      res <- 0.1
      flowData.temp <- Subset(flowData_transformed[j],time.gate)
      flowData.temp[[1]]@description$`$VOL` <- as.numeric(as.numeric(flowData_transformed[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j]))
      print(as.numeric(as.numeric(flowData_transformed[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j])))
      write.FCS(x=flowData.temp[[1]], filename=paste(i+(10*number),time.interval,paste(rownames(analysis.length)[j]), sep="_"), what="numeric")
    }
    setwd(old.wd)
  }
}

### Some small functions for averaging replicates
trip <- function(x,n=3){
  y = c(); s = c()
  j=1
  for(i in seq(1,length(x),by=n)){
    y[j] = mean(x[i:(i+n-1)])
    s[j] = sd(x[i:(i+n-1)])
    j=j+1
  }
  y=cbind(y,s)
  return(y)
}

trip_col <- function(x,n=3){
  y = c()
  j=1
  for(i in seq(1,length(x),by=n)){
    y[j] = x[(i+n-1)]
    j=j+1
  }
  return(y)
}

################################################################################
### Function for sampling to equal nr. of cells in flowset
### Standard is minimum number of cells 
################################################################################

FCS.resample <- function(x, sample=0, replace=FALSE){
  library(easyGgplot2)
  library(devtools)
  sample_distr <- data.frame(counts=fsApply(x,FUN=function(x) nrow(x),use.exprs=TRUE))
  p1 <- ggplot2.histogram(data=sample_distr , xName='counts',
                    fill="white", color="black",
                    linetype="longdash",addMeanLine=TRUE, meanLineColor="red",
                    meanLineType="dashed", meanLineSize=1)+
    theme_bw() + labs(y="Frequency", title="Original count distribution")
  if(sample==0) sample <- min(fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE))
  ## Remove all .fcs files with less observations than the specified sample
  x <- x[fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE)>sample]
  for(i in 1:length(x)){
    exprs(x[[i]]) <- exprs(x[[i]])[sample(1:nrow(exprs(x[[i]])), sample, replace=replace),]
  }
  print(p1)
  cat(paste0("Your samples were randomly subsampled to ",sample," cells"))
  return (x)
}


