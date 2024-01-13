###############################################################################
## Eploratory Analysis
##    This consists of two parts
##      - Dimintionality reduction
##      - Heatmaps
################################################################################

## Which exploratory analysis section would you like to include?
UMAP = TRUE

Heatmap = TRUE

# What is the path to the processed data?
data_path = "Data/Processed/NormalizedData.Rdata"

# Which assay to use? to see options use: target_demoData@assayData$
elt = "batch_corrected"

# Provide the variable names for ROI, AOI, and slide names
ROI = "region"

AOI = "SegmentLabel"

batch = "Experiment"

slideName = "SlideName"

######################################
## Load Processed Data
######################################

# You may need to edit this code to update path
load(data_path) 



######################################
## Sample Distribution
######################################

knitr::kable(target_demoData@phenoData@data %>%
               select(as.name(slideName),as.name(AOI),as.name(ROI)) %>%
               group_by(across(all_of(c(slideName, ROI, AOI))))%>%
               summarise(N=n()) %>%
               dcast(formula(paste0(slideName,"~",ROI,"+",AOI))),
             caption = "Distribution of Samples for each slide acrross ROI and (AOI)") 



######################################
## UMAP 
######################################

if(UMAP == TRUE){
    # update defaults for umap to contain a stable random_state (seed)
    custom_umap <- umap::umap.defaults 
    custom_umap$random_state <- 42 
    
    
    # run UMAP
    umap_out <-
      umap(t(assayDataElement(target_demoData , elt = elt)), 
           config = custom_umap)  # Uses default umap settings
    
    
    # Save the results from the UMAP in the data object
    pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)] 
    
    if(is.null(batch)){
      
    # Plot using ggplot
    ggplot(pData(target_demoData), 
           aes(x = UMAP1, y = UMAP2, color = !!as.name(ROI))) + 
      geom_point(size = 3) + 
      theme_bw() + ggtitle("Dim Reduction Using UMAP.")
      
    } else {
      # Plot using ggplot
      ggplot(pData(target_demoData), 
             aes(x = UMAP1, y = UMAP2, color = !!as.name(ROI), shape = !!as.name(batch))) + 
        geom_point(size = 3) + 
        theme_bw() + ggtitle("Dim Reduction Using UMAP.")
      
      
    }
}


######################################
## Heatmap
######################################

if(Heatmap == TRUE){

    # Log transform gene expression data. 
    mat <- target_demoData@assayData[[elt]]
    
    
    # Summarize genes by variance and get the top 100
    var_genes <- apply(mat,1,var)
    
    top50Genes <- var_genes[order(var_genes, decreasing = T)[1:50]]
    
    
    annot <- target_demoData@phenoData@data %>%
      select(all_of(c(ROI,AOI,slideName))) %>%
      arrange(!!as.name(slideName), !!as.name(AOI), !!as.name(ROI)) 
    
    
    
    
    pheatmap(mat[names(top50Genes),rownames(annot)], cluster_rows = F, cluster_cols = F,
             annotation_col = annot, show_rownames = F, show_colnames = F, 
             main = "Top 50 variable genes (log2-transformed) gene expression", border_color = "black" )
}
