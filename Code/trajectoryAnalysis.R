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

geeCKPlus <- geeglm(A2M~  regionT, family=gaussian, id=subject_ID, data = model_datPlus,
                            corstr = "exchangeable")

geeCKMinus <- geeglm(A2M~  regionT, family=gaussian, id=subject_ID, data = model_datMinus,
                              corstr = "exchangeable") 



summary(geeCKPlus)

summary(geeMpack)



# Create gee function to loop through all genes
trajectories <- function(dat,genes){
  # Create function for modeling each gene
      geeF <- function(X){
        
        geeM <-  geeglm(X~  regionT, family=gaussian, id=subject_ID, data =dat,
                        corstr = "exchangeable")
        
        # Get summary of results
        sum_res <- summary(geeM)
        
        # Output beta coeficient and pvalue
        c("Trajectory"= geeM$coefficients[2],
          se = sum_res$coefficients[2,2],
          "pvalue"= sum_res$coefficients[2,4],
        "intercept" = sum_res$coefficients[1,1]) 
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
trajectories_plus <- trajectories(model_datPlus,genes = genes[1:1000])


# Minus trajectories
trajectories_minus <-trajectories(model_datMinus,genes = genes[1:1000])




############################################################################
## Tables to summarise results
###########################################################################


as.data.frame(trajectories_plus) %>%
  mutate(
    # How many had a significant trajectory?
    sig = case_when(pvalue<0.05~"sig",
                    pvalue >= 0.05 ~ "NotSig"), 
    
    # UPward or downward trajectory?
    traject = case_when(
      Trajectory.regionT < 0  & sig == "sig" ~ "Downward",
      Trajectory.regionT > 0  & sig == "sig" ~ "Upward")
  ) %>%
  summarise(SigTrajecory = sum(sig=="sig"), sigProp = SigTrajecory /n(),
            Upward = sum(traject =="Upward", na.rm = T),
            Downward = sum(traject =="Downward", na.rm = T))



as.data.frame(trajectories_minus) %>%
  mutate(
    # How many had a significant trajectory?
    sig = case_when(pvalue<0.05~"sig",
                    pvalue >= 0.05 ~ "NotSig"), 
    
    # UPward or downward trajectory?
    traject = case_when(
      Trajectory.regionT < 0  & sig == "sig" ~ "Downward",
      Trajectory.regionT > 0  & sig == "sig" ~ "Upward")
  ) %>%
  summarise(SigTrajecory = sum(sig=="sig"), sigProp = SigTrajecory /n(),
            Upward = sum(traject =="Upward", na.rm = T),
            Downward = sum(traject =="Downward", na.rm = T))


#############################################################################
# Plot genes of interest
############################################################################

slopeG = trajectories_minus["ACADS",1]
interceptG = trajectories_minus["ACADS",4]

seG = trajectories_minus["ACADS",2]

upperSlopeG = slopeG + 1.96 * seG 
lowerSlopeG = slopeG - 1.96 * seG 

ggplot(model_datMinus, aes(x=regionT,y=ACADS))+
  geom_point() +
  geom_abline(slope = slopeG, intercept = interceptG) +
  geom_abline(slope = upperSlopeG, intercept = interceptG) + # Upper bound
  geom_abline(slope = lowerSlopeG, intercept = interceptG) # Lower Bound
  




