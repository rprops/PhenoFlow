#########################################################################
### Code for creating phenotypic fingerprint and its ecological indices
### Ruben Props, CMET, Ghent University
#########################################################################

require('flowFDA')
require("vegan")
require("MESS")
source("MRM.parameters.R")
set.seed(777)

### Insert here the path to your data
### Output files will be stored in this directory
path = "path"
### Import .fcs data
flowData <- read.flowSet(path = path, 
                         transformation = FALSE, pattern=".fcs")


### Select parameters (standard: two scatters and two FL) and 
### Transform data using the inverse hyperbolic sine
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                   `SSC-H`=asinh(`SSC-H`), 
                                   `FL3-H`=asinh(`FL3-H`), 
                                   `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

### Normalize data between [0,1] on average, 
### this is required for using the bw=0.01 in the fingerprint calculation
summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = mean(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))


### Create a PolygonGate for extracting the single-cell information
### Input coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8,8,14,14,3,7.5,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                   bw=0.01,normalize=function(x) x)

### Example of a fingerprint (x-axis: bins
### y-axis: normalized density ~ probability of cell being in this bin)
plot(fbasis@basis[1,],bg=adjustcolor("blue",0.7),col=adjustcolor("blue",0.7)
     ,cex=0.2)

### Calculate ecological parameters from normalized fingerprint 
### Densities will be normalized to the interval [0,1]
### n = number of replicates
### d = rounding factor
Diversity.fbasis <- Diversity(fbasis,d=3,n=1)
Evenness.fbasis  <- Evenness(fbasis,d=3,n=1)
Structural.organization.fbasis  <- So(fbasis,d=3,n=1)
Coef.var.fbasis  <- CV(fbasis,d=3,n=1)

### Export ecological data to .csv file in the chosen directory
write.csv2(file="results.metrics.csv",
           cbind(Diversity.fbasis, Evenness.fbasis,
                                          Structural.organization.fbasis,
                 Coef.var.fbasis))

#########################################################################
### Beta-diversity assessment of fingerprint
#########################################################################
beta.div <- beta.div.fcm(fbasis,n=1)
plot(beta.div)

#########################################################################
### Comparison to HNA/LNA fingerprint described by:
### Prest EI (2013). Monitoring microbiological changes in drinking water systems
### using a fast and reproducible flow cytometric method. Water Research.
#########################################################################

### Creating a rectangle gate, set correct threshold here for FL1
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/max,"FL3-H"=c(0,20)/max, 
                           filterId = "HNA bacteria")


### Check if rectangle gate is correct, if not readjust the above line
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)


### Extract the cell counts
a <- filter(flowData_transformed, rGate_HNA) 
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Exporting cell counts to .csv into chosen directory
write.csv2(file="results.counts.csv",
           data.frame(Samples=flowData_transformed@phenoData@data$name, 
                      Total.cells = TotalCount$true,HNA.cells = HNACount$true))

### Plotting Total cells, HNA cells and LNA cells
plot(TotalCount$true,pch=21,bg=adjustcolor("black",0.7),
     col=adjustcolor("black",0.7),cex=1.25,las=1,ylab="Counts (# cells)",
     xlab="Samples",ylim=c(0,max(TotalCount$true)))
points(HNACount$true,pch=21,bg=adjustcolor("red",0.7),
       col=adjustcolor("red",0.7),cex=1.25)
points((TotalCount$true-HNACount$true),pch=21,bg=adjustcolor("blue",0.7),
       col=adjustcolor("blue",0.7),cex=1.25)
legend("topright",legend=c("Total cells","HNA cells","LNA cells"), 
       pch=c(21,21,21),pt.bg=c(adjustcolor("black",0.7),adjustcolor("red",0.7),
                               adjustcolor("blue",0.7)),bty="n")

### Plotting % HNA cells vs diversity + testing for correlation (pearson)
plot(y=Diversity.fbasis$D2,x=HNACount$true/TotalCount$true,pch=21,
     bg=adjustcolor("black",0.7),col=adjustcolor("black",0.7),cex=1.25,
     las=1,ylab="Diversity (D2)",xlab="% HNA cells")
cor.test(Diversity.fbasis$D2,HNACount$true/TotalCount$true)

### Plotting % LNA cells vs diversity + testing for correlation (pearson)
plot(y=Diversity.fbasis$D2,x=(1-HNACount$true/TotalCount$true),pch=21,
     bg=adjustcolor("black",0.7),col=adjustcolor("black",0.7),cex=1.25,
     las=1,ylab="Diversity (D2)",xlab="% LNA cells")
cor.test(Diversity.fbasis$D2,(1-HNACount$true/TotalCount$true))

### Plotting diversity, % HNA and % LNA
plot(HNACount$true/TotalCount$true,pch=21,bg=adjustcolor("red",0.7),
     col=adjustcolor("red",0.7),cex=1.25,las=1,ylab="% NA cells",
     xlab="Samples",ylim=c(0,1.2))
### Plotting % LNA cells vs diversity
points((1-HNACount$true/TotalCount$true),pch=21,bg=adjustcolor("blue",0.7),
       col=adjustcolor("blue",0.7),cex=1.25,las=1,ylab='',xlab='')
par(new=TRUE)
plot(Diversity.fbasis$D2,bg=adjustcolor("black",0.7),
     col=adjustcolor("black",0.7),cex=1.25,las=1,ylab="",
     xlab="",xaxt="n",yaxt="n",type="l",lwd=2)
legend("topright",legend=c("Diversity","HNA cells","LNA cells"), 
       pch=c(21,21,21),pt.bg=c(adjustcolor("black",0.7),adjustcolor("red",0.7),
                               adjustcolor("blue",0.7)),bty="n")

