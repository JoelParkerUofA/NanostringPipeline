###############################################################################
## Eploratory Analysis
##    This consists of two parts
##      - Dimintionality reduction
##      - Heatmaps
################################################################################

UMAP = TRUE

Heatmap = TRUE



######################################
## Load Processed Data
######################################

# You may need to edit this code to update path
load("Data/Processed/NormalizedData.Rdata") 


# Which assay to use? to see options use: target_demoData@assayData$
elt = "batch_corrected"


######################################
## Sample Distribution
######################################

knitr::kable(target_demoData@phenoData@data %>%
               select(SlideName,SegmentLabel,region) %>%
               group_by(SlideName,region,SegmentLabel) %>%
               summarise(N=n()) %>%
               dcast(SlideName~region + SegmentLabel),
             caption = "Distribution of Samples for each slide acrross ROI and (AOI)") 



######################################
## UMAP 
######################################

if(Heatmap == TRUE){
    # update defaults for umap to contain a stable random_state (seed)
    custom_umap <- umap::umap.defaults 
    custom_umap$random_state <- 42 
    
    
    # run UMAP
    umap_out <-
      umap(t(assayDataElement(target_demoData , elt = elt)), 
           config = custom_umap)  # Uses default umap settings
    
    
    # Save the results from the UMAP in the data object
    pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)] 
    
    # Plot using ggplot
    ggplot(pData(target_demoData), 
           aes(x = UMAP1, y = UMAP2, color = region, shape = Experiment)) + 
      geom_point(size = 3) + 
      theme_bw() + ggtitle("Dim Reduction using UMAP.")

}


######################################
## Heatmap
######################################

if(Heatmap == TRUE){

    # Log transform gene expression data. 
    log2Trans <- target_demoData@assayData[[elt]]
    
    
    # Summarize genes by variance and get the top 100
    var_genes <- apply(log2Trans,1,var)
    
    top50Genes <- var_genes[order(var_genes, decreasing = T)[1:50]]
    
    
    annot <- target_demoData@phenoData@data %>%
      select(region,SegmentLabel, SlideName) %>%
      arrange(SlideName, SegmentLabel, region) 
    
    
    
    
    pheatmap(log2Trans[names(top50Genes),rownames(annot)], cluster_rows = F, cluster_cols = F,
             annotation_col = annot, show_rownames = F, show_colnames = F, 
             main = "Top 50 variable genes (log2-transformed) gene expression", border_color = "black" )
}
