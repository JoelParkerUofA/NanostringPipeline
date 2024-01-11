################################################################################
## The purpose of this script is to load and install the required packages 
## For the workflow
################################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes most up to date version of Bioc
#BiocManager::install(version="3.15")

#BiocManager::install("NanoStringNCTools")
#BiocManager::install("GeomxTools")
#BiocManager::install("GeoMxWorkflows")

# For preprocessing
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(knitr)
library(ggplot2)
library(ggforce)
library(dplyr)
library(scales) # for percent
library(reshape2)  
library(cowplot) 


# For batch correction
library(sva)

# For dim reduction
library(umap)
library(Rtsne)
library(pheatmap)  # for pheatmap

source(paste0("R/",dir("R")))

