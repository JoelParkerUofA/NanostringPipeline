###############################################################################
###  This script is used for the following steps
###  1. Data loading (.DCC files, PKC files and sample annotation)
##   2. Quality Control
##      - SegmentQC (default = ON)
##      - ProbeQC (default = ON)
##      - Limit of Quantification QC (LOQ) (default = OFF)
##   3. Normalization
##      - Q3 Normalization
##      - Background Normalization
##      - Transcripts per million (default = ON)
##
##  4. Log transformation (default = ON, base=2)
##
##  5. Batch correction (default = ON)
###############################################################################

# Which QC metrics would you like to use?
SegmentQC = T
ProbeQC = T
LOQ = F


# Which normalization would you like to use? 
# Options at q_norm (Q3 normalization ), b_norm (background normalization)
#. or TPM (transcripts per million)

norm_method = "TPM"


# Do you want to use a log transformation? If so what base do you want to use?
log_transformation = T

base = 2



# Do you want to include batch effect adjustment?
batchEff = T



# What are the ROI and AOI names in the annotation file

ROI = "region"

AOI = "SegmentLabel"

# Batch is only required if batchEff=T
batch = "Experiment"


DCC_col = "Sample_ID"

###############################################################################
### Step 1: 
### DCC files (folder), PKC Files (folder)
### and Sample annotation (folder)
###############################################################################


# Provide path to the DCC files, PKC file and Sample annotation
DCC = "../Data/Raw/DCC"

PKC = "../Data/Raw/PKC"

SampleAnnotation = "../Data/Raw/Annotation"


DCC_col = "Sample_ID"

######################################
#### Get all paths for the raw data
######################################

# DCC
DCCFiles <- dir(DCC, pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)


# PKC
PKCFiles <- dir(PKC, pattern = ".pkc$", full.names = TRUE, recursive = TRUE)


# Sample annotation (Alternatively you could provide a direct path)
SampleAnnotationFile <- dir(SampleAnnotation, pattern = ".xlsx$",
                            full.names = TRUE, recursive = TRUE)

#########################################
### Load all data into the nanostring 
### object to store all data together.
#########################################
# 1. load data (This may take a few moments)
demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                   pkcFiles = PKCFiles,
                                   phenoDataFile = SampleAnnotationFile,
                                   phenoDataSheet = "Sheet1",
                                   phenoDataDccColName = DCC_col
)



################################################################################
### Step 2: Quality Control
###   Consists of the following subparts
###     - Count shift
###     - Segment QC
###     - Probe QC
###     - LOQ (Limit of quantification filtering
###############################################################################

##############################################################
# Shift counts to one
#############################################################
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)



############################################################
## Segment QC
############################################################
if(SegmentQC==T){
  # Set the QC parameters that we wish to use. These can be edited.  
  QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads (1000)
         percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
         percentStitched = 80,   # Minimum % of reads stitched (80%)
         percentAligned = 80,    # Minimum % of reads aligned (80%)
         percentSaturation = 50, # Minimum sequencing saturation (50%)
         minNegativeCount = 1,   # Minimum negative control counts (1)
         maxNTCCount = 9000,     # Maximum counts observed in NTC well (9000)
         minNuclei = 20,         # Minimum # of nuclei estimated (20)
         minArea = 1000)         # Minimum segment area (1000)
  
  
  # Include the QC parameters in the nanostring object 
  demoData <- setSegmentQCFlags(demoData, 
                                qcCutoffs = QC_params)   
  
  ##############################################
  ###  Create QC summary table based on
  ###  these flags.
  ##############################################
  
  QCResults <- protocolData(demoData)[["QCFlags"]]
  flag_columns <- colnames(QCResults)
  
  
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                           Warning = colSums(QCResults[, flag_columns]))
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  
  
  # Display summary of the results
  # print( knitr::kable(QC_Summary["TOTAL FLAGS", ] <-
  #  c(sum(QCResults[, "QCStatus"] == "PASS"),
  #     sum(QCResults[, "QCStatus"] == "WARNING"))))
  
  
  ################################################
  ### Visualize segment QC
  ################################################
  
  col_by <- AOI
  
  # Trimmed % 
  print( QC_histogram(sData(demoData), "Trimmed (%)", col_by, QC_params$percentTrimmed))
  
  # Stitched %
  print(QC_histogram(sData(demoData), "Stitched (%)", col_by, QC_params$percentStitched))
  
  # Aligned %
  print(QC_histogram(sData(demoData), "Aligned (%)", col_by, QC_params$percentAligned))
  
  # Saturated %
  print(QC_histogram(sData(demoData), "Saturated (%)", col_by, QC_params$percentSaturation) +
          labs(title = "Sequencing Saturation (%)",
               x = "Sequencing Saturation (%)"))
  
  # Area
  print(QC_histogram(sData(demoData), "AOISurfaceArea", col_by, QC_params$minArea,
                     scale_trans = "log10"))
  
  # Nuclei
  print(QC_histogram(sData(demoData), "AOINucleiCount", col_by, QC_params$minNuclei))
  
  
  #################################################
  ## Calculate and plot geometric mean of the 
  ## Negative probes for each segment.
  #################################################
  
  # gene annotations modules
  #  Assign the pkcs variable to the name of the gene annotation file 
  pkcs <- annotation(demoData)
  
  
  #  Create table with gene annotation file names
  modules <- gsub(".pkc", "", pkcs)
  
  
  #  Calculate the geometric means of the negative controls for each of the samples.
  negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
           assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
  
  
  # Store the geometric means in the protocol (metadata)
  protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans
  
  
  #  Explicitly copy the Negative geoMeans from sData to pData
  negCols <- paste0("NegGeoMean_", modules)
  pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
  
  
  #  Plot the geoMetric means for each sample.  
  for(ann in negCols) {
    plt <- QC_histogram(pData(demoData), ann, col_by, QC_params$minNegativeCount, scale_trans = "log10")
    print(plt)
  }
  
  
  
  ###############################
  ## Display Segment QC results
  ## and filter segments
  ###############################
  
  # Segment QC Summary
  print(kable(QC_Summary, caption = "QC Summary Table for each Segment"))
  
  # Filter only segments that passed 
  # Subset data only on the samples that passed the QC Metrics. 
  demoData <- demoData[, QCResults$QCStatus == "PASS"]
  
  
}   

####################################################################
###  Probe QC
####################################################################

if(ProbeQC == T){
  
  # Set flag values
  demoData <- setBioProbeQCFlags(demoData, 
                                 qcCutoffs = list(minProbeRatio = 0.1,
                                                  percentFailGrubbs = 20), 
                                 removeLocalOutliers = TRUE)
  
  
  # 2. Data set with raised flaggs for each probe. 
  ProbeQCResults <- fData(demoData)[["QCFlags"]]
  
  
  # 3. Table with flag summary count. 
  qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                      Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                      Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                  & !ProbeQCResults$GlobalGrubbsOutlier))
  
  print(kable(qc_df, caption = "Probes flagged or passed as outliers"))
  
  
  ##################################
  ### Filter Probes
  ##################################
  
  #Subset object to exclude all that did not pass Ratio & Global testing
  ProbeQCPassed <- 
    subset(demoData, 
           fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
             fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
  
  # Resign to demoData
  demoData <- ProbeQCPassed 
  
}    

##################################
## Create gene level count matrix
##################################


# Aggregate counts with the geometric mean. 
target_demoData <- aggregateCounts(demoData)



################################################################
###  Limit of Quantification (LOQ)
################################################################
if(LOQ==TRUE){
  ################################
  ## Calculate LOQ
  ###############################
  # Define the LOQ cutoff and the minLOQ. 
  cutoff <- 2
  minLOQ <- 2
  
  
  
  # Calculate LOQ per module tested
  LOQ <- data.frame(row.names = colnames(target_demoData))
  for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
      LOQ[, module] <-
        pmax(minLOQ,
             pData(target_demoData)[, vars[1]] * 
               pData(target_demoData)[, vars[2]] ^ cutoff)
    }
  }
  
  
  # Store the LOQ information.
  pData(target_demoData)$LOQ <- LOQ
  
  
  #####################################
  ### Create matrix with LOQ
  ### Pass/Fail information
  #####################################
  
  
  
  #  Initialize an empty vector to store the LOQ_Mat values.
  LOQ_Mat <- c()
  
  # Run for loop.
  for(module in modules) {
    ind <- fData(target_demoData)$Module == module
    
    # 2. Check if each value is greater than the corresponding LOQ value.
    Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                       FUN = function(x) {
                         x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
  }
  
  
  # 3. Ensure ordering since this is stored outside of the geomxSet
  LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]
  
  
  #####################################
  ## Segment gene detection rate plot
  #####################################
  
  # 1. Save detection rate information to pheno data
  pData(target_demoData)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
  pData(target_demoData)$GeneDetectionRate <-
    pData(target_demoData)$GenesDetected / nrow(target_demoData)
  
  
  # 2. Classify each sample based on the percentage of genes above the LOQ. 
  pData(target_demoData)$DetectionThreshold <- 
    cut(pData(target_demoData)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  
  
  # 3. Create a barplot of the different classes. (1%, 5%, 10%, 15%)
  
  geneplot <- ggplot(pData(target_demoData),
                     aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = get(ROI))) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type", title = "Number of Samples in Gene Detection Rate Samples")
  
  print(geneplot)
  
  #######################################
  ### Filter segments based on 
  ### gene detection rate
  #######################################
  
  target_demoData <- target_demoData[, pData(target_demoData)$GeneDetectionRate >= .01]
  
  
  #######################################
  ## Gene gene detection rate plot
  #######################################
  
  # Filter LOQ matrix
  LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
  
  
  # Calculate the number of samples detected for each sample.
  fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
  
  
  # Save the dection rate for each gene.
  fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))
  
  
  
  # Create a data frame with for plotting with a frequency column.
  plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
  
  
  # Calculate number of genes above each threshold.
  plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
  
  
  # Calculate the percentage of genes above this threshold. 
  plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
  rownames(plot_detect) <- plot_detect$Freq
  
  
  # Plot gene detection rate
  print( ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
           geom_bar(stat = "identity") +
           geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
                     vjust = 1.6, color = "black", size = 4) +
           scale_fill_gradient2(low = "orange2", mid = "lightblue",
                                high = "dodgerblue3", midpoint = 0.65,
                                limits = c(0,1),
                                labels = scales::percent) +
           theme_bw() +
           scale_y_continuous(labels = scales::percent, limits = c(0,1),
                              expand = expansion(mult = c(0, 0))) +
           labs(x = "% of Segments",
                y = "Genes Detected, % of Panel > LOQ")
         
  )
  # Save plot
  #ggsave("../Outputs/Plots/LOQGeneDetectionRatePlot.png")
  
  
  #########################################
  ### Filter genes based on gene detection
  ### Rate
  ##########################################
  #  Subset to target genes detected in at least 10% of the samples.
  negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)
  
  
  # Filter out genes detected in less than 10% of genes. 
  target_demoData <- 
    target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                      fData(target_demoData)$TargetName %in% neg_probes, ]
  
}    


################################################################################
### Step 3: Normalization
###   There are two types of normalization possible
###     - Q3 normalization
###     - Background Normalization
###     - Transcripts per million (default)
################################################################################

##################################
### Q3 Normalization
##################################

if(norm_method == "q_norm"){
  
  # Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
  target_demoData <- normalize(target_demoData ,
                               norm_method = "quant", 
                               desiredQuantile = .75,
                               toElt = norm_method)
  
}

###################################
## Background normalization
###################################

if(norm_method == "b_norm"){
  
  # Background normalization for WTA/CTA without custom spike-in
  target_demoData <- normalize(target_demoData ,
                               norm_method = "neg", 
                               fromElt = "exprs",
                               toElt = norm_method)
  
  
}


#############################################
## Transcripts per-million standardization
#############################################

if(norm_method== "TPM"){
  
  
  TPM_norm = apply(target_demoData@assayData$exprs, 2,
                   FUN = function(x){ (x/sum(x))*1000000 })
  
  
  # Add assay data
  assayDataElement(target_demoData,"TPM") = TPM_norm
  
  
}



##############################################
### Log transformation
##############################################

if(log_transformation ==T) {
  
  
  # log transform assay data
  assayLogT = log(target_demoData@assayData[[norm_method]], base = base)
  
  
  # Add assay data to the 
  assayDataElement(target_demoData,paste0("logT-base",base)) = assayLogT
  
  
}





################################################################################
## Batch Effects (here batch info is captured in "Expiriment variable)
##. There are two cases for this because if we used log transformation then
## we want to use the log transformation assay slot. 
################################################################################
if(batchEff==T){
  
  # If we didnt log transform then we want to use the norm_method slot. 
  if(log_transformation==F){
    
    # Do batch correction
    adjusted_expression <- ComBat(dat = target_demoData@assayData[[norm_method]],
                                  batch = target_demoData@phenoData@data$Experiment,mod = NULL, prior.plots = F)
    
    # Save batch corrected results in nanostring object
    assayDataElement(target_demoData,"batch_corrected") = adjusted_expression
  }
  
  
  # if log transformed, we will want to you to log transformation assay slot
  if(log_transformation==T){
    
    # Do batch correction
    adjusted_expression <- ComBat(dat = target_demoData@assayData[[paste0("logT-base",base)]],
                                  batch = target_demoData@phenoData@data$Experiment,mod = NULL, prior.plots = F)
    
    # Save batch corrected results in nanostring object
    assayDataElement(target_demoData,"batch_corrected") = adjusted_expression
  }
  
  
}



##################################
## Save Data
##################################
save(target_demoData,file = "../Data/Processed/NormalizedData.Rdata")
