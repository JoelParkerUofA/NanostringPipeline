#############################################################################
## Trajectory analysis. Here we are testing for a trajectory from no disease to 
## Disease. 
#############################################################################


# What is the path to the processed data?
data_path = "Data/Processed/NormalizedData.Rdata"

# Which assay to use? to see options use: target_demoData@assayData$ 
# here we are using the batch corrected log-transformed
elt = "batch_corrected"

# Provide the variable names for ROI, AOI, and slide names
ROI = "region"

AOI = "SegmentLabel"

subject_id = "SlideName"


######################################
## Load Processed Data
######################################

# You may need to edit this code to update path
load(data_path) 



######################################
## Sampling distributions
######################################
target_demoData@phenoData@data %>%
               select(as.name(slideName),as.name(AOI),as.name(ROI)) %>%
               group_by(across(all_of(c(slideName, ROI, AOI))))%>%
               summarise(N=n()) %>%
               dcast(formula(paste0(slideName,"~",ROI,"+",AOI)))





##############################################
## Create model data so gene expression
## and metadata are together in one dataset
##############################################

# Get the assay data for the modeling
# CPM
assay <- data.frame(t(assayDataElement(target_demoData, elt = elt)))


# Here we merge the metadata and the assay data into one dataframe and add the 
# trajectory to the regions. 
model_dat <- merge(target_demoData@phenoData@data, assay, by="row.names") %>%
  mutate(
    regionT =  case_when(region == "Sun Protected" ~1,
                         region=="AK Edge"~2,
                         region == "AK Center"~ 3), 
    
    subject_ID = as.factor(!!as.name(subject_id))
  )


# Split the data into CKPlus and CKMinus 
model_datPlus <- model_dat %>%
  filter(SegmentLabel=="panckPos")


model_datMinus <- model_dat %>%
  filter(SegmentLabel=="panckNeg")



# Run trajectory analysis on one gene
model_datPlus = model_datPlus %>% 
  arrange(subject_ID)


geeM <- gee(A2M~regionT, data = dat,id=subject_ID, 
            corstr = "independence")

geeMpack <- geepack::geeglm(A2M~  regionT, family=gaussian, id=subject_ID, data = model_datPlus,
                            corstr = "exchangeable")
    

summary(geeM)

summary(geeMpack)



# Create gee function to loop through all genes
trajectories <- function(dat,genes){
  # Create function for modeling each gene
      geeF <- function(X){
        geeM <- gee(X~regionT, data = dat,id=subject_ID, 
                    corstr = "independence")
        
        # Output beta coeficient and pvalue
        c("Trajectory"= geeM$coefficients[2],
          pvalue= 2*pnorm(abs(summary(geeM)$coefficients[2,5]), lower.tail = F))
      }
      
      # Orderdata based on subject ID
      dat = dat %>%
        arrange(subject_ID)
      
      # Run gene trajectoris on the data
      traj <- apply(dat[,genes], MARGIN = 2, FUN = geeF)
      
      # Return the result
      return(t(traj))    
        
}


#################################################################
## Run trajectory analysis
################################################################

# Identify the genes we want to run it on. 
genes <-colnames(assay)


# Plus trajectories
trajectories_plus <- trajectories(model_datPlus,genes = genes[1:10])

trajectories_plus <- t(trajectories_plus)

# Minus trajectories
trajectories_minus <-trajectories(model_datMinus,genes = genes[1:10])




#############################################################################
ggplot(model_datMinus, aes(x=regionT,y=CST3))+
  geom_point() +
  geom_smooth(method = "lm")




