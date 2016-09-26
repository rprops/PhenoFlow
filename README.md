#PhenoFlow _(v1.1)_
##Phenotypic alpha-diversity and beta-diversity parameters for flow cytometry data of microbial communities
===============
- **Authors**: Ruben Props [Ruben.Props@UGent.be], Pieter Monsieurs, Mohamed Mysara, Lieven Clement, Nico Boon

Accompanying code for <a href="http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12607/full"> Props et al. (2016), Measuring the biodiversity of microbial communities by flow cytometry, *Methods in Ecology and Evolution*, doi:10.1111/2041-210X.12607</a>

<p align="justify">The goal of this method is to translate (longitudinal) flow cytometry data, into diversity estimates of microbial communities, such as depicted below (data available at <a href="http://datadryad.org/resource/doi:10.5061/dryad.m1c04"> doi:10.5061/dryad.m1c04 </a>), which can be compared with 16S amplicon data sets. </p>
![alt text][logo]
[logo]: https://github.com/rprops/PhenoFlow/blob/master/Animation_low_res.gif "Figure 1"

## Installation of required packages
From CRAN, the following packages need to be installed:
```R
install.packages("mclust")
install.packages("vegan")
install.packages("MESS")
install.packages("multcomp")
install.packages("KernSmooth")
install.packages("mvtnorm")
install.packages("lattice")
install.packages("survival")
install.packages("TH.data")
install.packages("gridExtra")
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
source("https://bioconductor.org/biocLite.R")
biocLite("flowViz")
```
<p align="justify">Next, install the flowFDA package from https://github.com/lievenclement/flowFDA. Or, alternatively download the package from this repository and unzip it directly into the R library folder.</p>

```R
install.packages("devtools")
library(devtools)
install_github("lievenclement/flowFDAExampleData")
install_github("lievenclement/flowFDA", build_vignettes=TRUE)
```

Also install the folowing packages for easy count histogram plots:
```R
install_github("easyGgplot2", "kassambara")
```
## Scripts and functions
Script  | Content
------------| -----------
MRM.parameters.R | Contains the functions used for calculating the phenotypic diversity
Fingerprint.R | Code for assessing the phenotypic diversity and cell counts from flow cytometry tutorial data

Functions  | Actions
------------| -----------
flowBasis | Function part of the flowFDA package (in development) which performes bivariate kernel density estimations on the phenotypic parameters
trip/trip_col | Small functions for taking replicate averages and the corresponding sample names
Diversity | Calculation of Hill diversities of order 0, 1 and 2 from fingerprint object
Evenness | Calculation of pareto evenness (Wittebolle L. et al. (2009)) from fingerprint object (1 = maximum evenness, 0 = minimum evenness)
cum_Richness | Used in Evenness function for calculating cumulative density functions
So | Calculation of Structural Organization parameter (Koch et al. (2014), Frontiers in Microbiology)
CV | Calculation of Coefficient of Variation (CV) of the fingerprint object
beta.diversity.fcm | Non-metric Multidimensional Scaling (NMDS) or PCoA of the phenotypic fingerprints
dist.fcm | Calculating distance matrix between fingerprints
time.discretization | Function for subsetting .fcs files in time intervals and exporting them as new .fcs files. Designed for the analysis of on-line experiments.
FCS.resample | Resamples sample files from flowSet object to an equal number of cells. Standard is to the minimum sample size.

## Input required by the user

- Path to .fcs files and path to where the output excel file is to be put.

- Gating strategy for isolating the cellular information and discarding the instrument/sample noise.

Full tutorial data is available at: https://flowrepository.org/id/RvFr3eLf9W5MNLBtMv8cQG41U05HcOr1pQ8pTpnFKBYfHAjTiYpbWuweKbSD3mQF
## Output
- CSV files with the diversity parameters, total cell counts, high nucleic acid (HNA) and low nucleic acid content (LNA) cell counts in the specified gate(s).  

<p align="justify"><b>Note: different flow cytometers may use different nomenclature for detector signals, e.g., FL1-H can be named FL1 log in the .fcs file. The fingerprint.R script will have to be adjusted accordingly.</b></p>


# How to use the scripts

## Load the packages and source code
Open RStudio. Make sure that your working directory is located where the scripts are saved.
```R
require('flowFDA')
require("vegan")
require("MESS")
```
The source code MRM.parameters.R contains all the necessary functions.
```R
source("MRM.parameters.R")
```
Set a fixed seed to ensure reproducible analysis
```R
set.seed(777)
```
## Data analysis
Insert the path to your data folder, for example test_data for the tutorial data.
```R
path = "test_data"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")
```
<p align="justify">At this point we select the phenotypic features of interest and transform their intensity values according to the hyperbolic arcsin. In this case we chose two fluorescent parameters and two scatter parameters in their height format (-H). Depending on the FCM, the resolution may increase by using the area values (-A) since many detectors have a higher signal resolution for area values. For transparency we store the transformed data in a new object, called <code>flowData_transformed</code>. Due to filtering of relevant parameters, we also reduce the data size of the <code>flowSet</code>. This becomes relevant for larger datasets. For example, a dataset of 200 samples of an average of 25 000 cells will require 200 - 300 MB of RAM.</p>

```R
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                   `SSC-H`=asinh(`SSC-H`), 
                                   `FL3-H`=asinh(`FL3-H`), 
                                   `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)
```
<p align="justify">Now that the data has been formatted, we need to discard all the signals detected by the FCM which correspond to instrument noise and (in)organic background. This is done by selecting the cells in a scatterplot on the primary fluorescence or scatter signals. For SYBR Green I, this is done based on the <code>FL1-H</code> and <code>FL3-H</code> parameters. For this example, an initial polygon gate (<code>polyGate1</code>) is created and adjusted based on the sample type in question. For each specific experiment, it is advised to use identical gating for each sample. The gating strategy is evaluated on the <code>xyplot</code> and adjusted if necessary. A more detailed guideline for gating or denoising aqueous microbial samples can be found <a href="http://jornades.uab.cat/workshopmrama/sites/jornades.uab.cat.workshopmrama/files/Assessing_water_quality_with_the_BD_Accuri_C6_flow_cytometer.pdf">here</a> (p.6):</p>
```R
### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
```
Here is an example of a good and bad filtering approach:
![My image](https://cloud.githubusercontent.com/assets/19682548/16420078/44b0d5ec-3d50-11e6-800d-10a1e3413ca2.png)
When the optimal gate has been chosen, the data can be denoised using the `Subset` function.  

```R
### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)
```
<p align="justify">Next, the phenotypic intensity values of each cell are normalized to the [0,1] range. This is required for using a bandwidth of 0.01 in the fingerprint calculation. Each parameter is normalized based on the maximum FL1-H intensity value over the data set.</p>

```R
summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))
```

<p align="justify">The denoised data can now be used for calculating the phenotypic fingerprint using the <code>flowBasis</code> function. Changing <code>nbin</code> increases the grid resolution of the density estimation but also steeply increases the computation time. The most determining factor for getting an accurate kernel density estimation is the number of cells counted in <code>polyGate1</code>. In general, 1,000 cells will give a good estimation of the mean alpha-diversity (D<sub>2</sub>) but in order to reduce the variance on this estimate I would reccomend a cell count of 10,000 cells or more. 

If desired, you can randomly resample your samples to the lowest sample size or any user specified sample size using the <code>FCS.resample()</code> function. Samples with a sample size that is equal to 0 or that is lower than the specified size will be discarded. </p>
```R
### Randomly resample to the lowest sample size
flowData_transformed <- FCS.resample(flowData_transformed)
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                   bw=0.01,normalize=function(x) x)
```
<p align="justify">From the phenotypic fingerprint, alpha diversity metrics can be calculated. <code>n</code> is the number of technical replicates, <code>d</code> is a rounding factor which is used to eliminate unstable density values from the dataset. Different rounding factors usually only scale the diversity estimates by a fixed factor and do not affect temporal trends or comparative analysis.</p>
```R
### Calculate ecological parameters from normalized fingerprint 
### Densities will be normalized to the interval [0,1]
### n = number of replicates
### d = rounding factor
Diversity.fbasis <- Diversity(fbasis,d=3,n=1,plot=FALSE)
Evenness.fbasis <- Evenness(fbasis,d=3,n=1,plot=FALSE)
Structural.organization.fbasis <- So(fbasis,d=3,n=1,plot=FALSE)
Coef.var.fbasis <- CV(fbasis,d=3,n=1,plot=FALSE)
```

Add the argument `plot=TRUE` in case a quick plot of the D<sub>2</sub> diversity values is desired.
![diversity example](https://cloud.githubusercontent.com/assets/19682548/17303302/17b92b4e-57ee-11e6-954a-cea341c5a152.png)

Alpha diversity analysis has completed: time to export all the data to your working directory. If you are not sure where this is, type <code>getwd()</code>.  
```R
### Export ecological data to .csv file in the chosen directory
write.csv2(file="results.metrics.csv",
           cbind(Diversity.fbasis, Evenness.fbasis,
                                          Structural.organization.fbasis,
                 Coef.var.fbasis))
```
Optionally, you can also perform a beta diversity analysis using Non-metric Multidimensional Scaling (NMDS) or PCoA from the <code>vegan</code> package. The <code>beta.div.fcm</code> function does the calculation (accepts all dissimilarity metrics of the <code>vegdist</code> function under the <code>dist</code> argument) and the <code>plot.beta.fcm</code> plots the ordination. Additional factorial arguments to be given to the plot function are: <code>color</code> and <code>shape</code> which can be used for visualization (<code>color</code> and <code>shape</code> have to be factor vectors with length equal to number of samples). The title of the legends for these factors can be specified in the <code>labels</code> argument. In case no legend has to be displayed, <code>legend.pres</code> can be put to <code>FALSE</code>.
<p align="justify"><b> Attention: this dimension reduction utilizes the density values as opposed to the count values which are used for beta-diversity analysis in community 16S data sets</b></p>

```R
### Beta-diversity assessment of fingerprint
beta.div <- beta.div.fcm(fbasis,n=1,ord.type="PCoA")
plot.beta.fcm(beta.div,legend.pres=FALSE)
```
![pcoa example](https://cloud.githubusercontent.com/assets/19682548/17302990/bc8d877a-57ec-11e6-8709-4f5d5dfee679.png)

<p align="justify">It is often also useful to know the exact cell densities of your sample. This is performed by the following code. Additionally it quantifies the amount of High Nucleic Acid (HNA) and Low Nucleic Acid (LNA) bacteria as defined by <a href="http://www.sciencedirect.com/science/article/pii/S0043135413008361">Prest et al. (2013)</a>.</p>

<p align="justify"><b>Warning: the HNA/LNA partition is only valid for data gathered on a BD C6 Accuri flow cytometer.
For other flow cytometers the threshold (FL1-H = 20 000) should be adjusted according to the appropriate reference samples.</b></p>
```R
### Creating a rectangle gate for counting HNA and LNA cells
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/max,"FL3-H"=c(0,20)/max, 
                           filterId = "HNA bacteria")
### Normalize total cell gate
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
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

### Exporting cell counts to .csv file to working directory
write.csv2(file="results.counts.csv",
           data.frame(Samples=flowData_transformed@phenoData@data$name, 
                      Total.cells = TotalCount$true,HNA.cells = HNACount$true))
```
===============
