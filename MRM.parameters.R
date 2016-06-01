require("vegan")
require("MESS")
require("flowFDA")

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


### Diversity from FCM data (Hill numbers D0, D1 and D2)
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### n = number of replicates
Diversity <- function(x,d=4,n=3, plot=FALSE){
  x<- x@basis/apply(x@basis, 1, max)
  D=matrix(ncol=3,nrow=nrow(x))
    D[,1] = apply(x,1,FUN=function(x) {
      x = round(x,d);x<-x[x!=0];sum(x!=0)
    })
    D[,2] = apply(x,1,FUN=function(x){
      x = round(x,d);x<-x[x!=0];exp(-sum((x/sum(x))*log(x/sum(x))))
    })
    D[,3] = apply(x,1,FUN=function(x){
      x = round(x,d);x<-x[x!=0];1/sum((x/sum(x))^2)
    })
  if(n>1){
    D2=matrix(ncol=6,nrow=length(x[,1])/n)
    D2[,1] = trip(D[,1],n)[,1]
    D2[,2] = trip(D[,1],n)[,2]
    D2[,3] = trip(D[,2],n)[,1]
    D2[,4] = trip(D[,2],n)[,2]
    D2[,5] = trip(D[,3],n)[,1]
    D2[,6] = trip(D[,3],n)[,2]
    results = data.frame(D2)
    colnames(results) = c("D0","sd.D0","D1","sd.D1","D2",
                          "sd.D2")
    rownames(results) = trip_col(attr(x,"dimnames")[[1]],n=n)
  }
  else{
    D2=matrix(ncol=3,nrow=length(x[,1]))
    D2[,1] = D[,1]
    D2[,2] = D[,2]
    D2[,3] = D[,3]
    results = data.frame(D2)
    colnames(results) = c("D0","D1","D2")
    rownames(results) = attr(x,"dimnames")[[1]]
  }
    if (plot==TRUE) {
      plot(results$D2,pch=21,bg=adjustcolor("blue",0.7),
                        col=adjustcolor("blue",0.7),cex=1.5,las=1,ylab="Diversity",xlab="Samples")
    }
    return(results)
}

### Evenness based on cumulative paretocurves
### x = flowBasis object from fingerprint (e.g., fingerprint)
### d = threshold for denoising
### The higher the evenness the higher the value
### n = number of replicates
Evenness <- function(x,d=3,n=3,plot=FALSE){
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
So <- function(x,d=3,n=3,plot=FALSE){
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
CV <- function(x,d=3,n=3,plot=FALSE){
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
beta.div.fcm <- function(x, d=3, n=3, dist="bray",k=2,iter=100){
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
    mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
  }
  else{
    input.dist <- vegdist(x,method=dist)
    mds.fbasis <- metaMDS(input.dist,autotransform=FALSE, k,trymax=iter)
  }
  stressplot(mds.fbasis)
  return(mds.fbasis)
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

################################################################################
### Function for sampling to equal nr. of cells in flowset
### Standard is minimum number of cells 
################################################################################

FCS.resample <- function(x, sample=0,replace=FALSE){
  ## Remove all .fcs files with 0 observations
  x <- x[fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE)!=0]
  if(sample==0) sample <- min(fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE))
  for(i in 1:length(x)){
    exprs(x[[i]]) <- exprs(x[[i]])[sample(1:nrow(exprs(x[[i]])), sample, replace=replace),]
  }
  return (x)
}

