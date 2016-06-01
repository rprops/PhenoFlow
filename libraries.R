############################################
#### Install required packages
############################################

install.packages("mclust")
install.packages("vegan")
install.packages("MESS")
install.packages("multcomp")
install.packages("KernSmooth")
install.packages("mvtnorm")
install.packages("lattice")
install.packages("survival")
install.packages("TH.data")

source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

source("https://bioconductor.org/biocLite.R")
biocLite("flowViz")

require('flowFDA')
require("vegan")
require("MESS")
