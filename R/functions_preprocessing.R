# Preprocessing functions

#' Generate synthetic additive groups
#'
#' Five synthetic additive groups were manually constructed based on the advice
#' of SSA analysts, using reference_group numbers:
#' 11) Literacy required
#' 89) Pre-employment training: Certification
#' 90) Pre-employment training: License
#' 91) Pre-employment training: Educational Certificate
#' 78) Sitting or standing/walking
#'
#' @param blsRaw Original data (read in from CSV file)
#' @return Data edited for new additive groups, and containing only relevant
#' records (mean and standard error estimates) and columns
#' @export
syntheticAddGroups <- function(blsRaw) {
  
  blsRaw[which(blsRaw$data_element_text=="Literacy required" &
                 (blsRaw$data_type_text=="YES" | blsRaw$data_type_text=="NO") &
                 blsRaw$estimate_type=="Percent"),"additive_group"] <- 11
  blsRaw[which(blsRaw$data_element_text=="Pre-employment training: Certification" &
                 (blsRaw$data_type_text=="YES" | blsRaw$data_type_text=="NO") &
                 blsRaw$estimate_type=="Percent"),"additive_group"] <- 89
  blsRaw[which(blsRaw$data_element_text=="Pre-employment training: License" &
                 (blsRaw$data_type_text=="YES" | blsRaw$data_type_text=="NO") &
                 blsRaw$estimate_type=="Percent"),"additive_group"] <- 90
  blsRaw[which(blsRaw$data_element_text=="Pre-employment training: Educational Certificate" &
                 (blsRaw$data_type_text=="YES" | blsRaw$data_type_text=="NO") &
                 blsRaw$estimate_type=="Percent"),"additive_group"] <- 91
  
  # In original analysis we used '100' for the new additive group
  # '100' is actually already taken by 'Pace: Pause Control' (a COG requirement)
  # Instead for the package we use '78', i.e. the reference group of 'Sitting'
  blsRaw[(blsRaw$data_element_text=="Sitting" | blsRaw$data_element_text=="Standing/walking") &
           blsRaw$unit_of_measure=="Percent of day" & blsRaw$estimate_type=="Mean","additive_group"] <- 78
  
  # Select only required observations and columns
  bls <- blsRaw[!is.na(blsRaw$additive_group),c("occupation_text","upper_soc_code","data_element_text",
                                                "data_type_text","estimate_type_text","additive_group","value")]
  
  return(bls)
}


#' Generate missing observations
#'
#' Observations that were missing in the data (i.e., no associated record or 
#' estimate) are generated here.
#' This is done by "completing" the data using tidyr::complete.
#'
#' @param bls Data frame of relevant observations (output of syntheticAddGroups())
#' @return Data with newly generated observations added (completed data)
#' @export
fillMissingObservations <- function(bls) {
  blsC = bls%>% tidyr::complete(occupation_text,nesting(data_type_text,additive_group,
                                                 estimate_type_text,data_element_text))
  
  for (occ in unique(blsC$occupation_text)) {    # fill in upper SOC code for new observations
    usoc <- as.numeric(names(table(blsC[blsC$occupation_text==occ,"upper_soc_code"])))
    if (length(usoc)==1) {
      blsC[blsC$occupation_text==occ,"upper_soc_code"] <- usoc
    }
  }
  
  return(blsC)
}


#' Generate additional (optional) missing observations
#'
#' SSA analysts requested the addition of a few other missing observations,
#' based on requirement/level, listed below:
#' 
#' Humidity: CONSTANTLY
#' Heavy vibrations: CONSTANTLY, FREQUENTLY, SELDOM
#' Extreme cold: CONSTANTLY, FREQUENTLY
#' Climbing ladders, ropes, or scaffolds: CONSTANTLY, FREQUENTLY
#' Climbing ramps or stairs (work-related): CONSTANTLY, FREQUENTLY
#' Stooping is missing: CONSTANTLY
#' Kneeling is missing: CONSTANTLY, FREQUENTLY
#' Crawling is missing: CONSTANTLY, FREQUENTLY
#'
#' @param blsC Completed data (output of fillMissingObservations())
#' @return Completed data, with requested requirement/levels added
#' @export
otherMissingObservations <- function(blsC) {
  # Humidity
  hu <- droplevels(blsC[which(blsC$data_element_text=="Humidity" &
                                blsC$data_type_text=="NOT PRESENT"),])
  hu$data_type_text <- "CONSTANTLY"
  hu$value <- NA
  new.obs <- hu
  rm(hu)
  
  # Heavy vibrations
  hv <- droplevels(blsC[which(blsC$data_element_text=="Heavy vibrations" &
                                blsC$data_type_text=="NOT PRESENT"),])
  hv$data_type_text <- "CONSTANTLY"
  hv$value <- NA
  new.obs <- rbind(new.obs,hv)
  hv$data_type_text <- "FREQUENTLY"
  hv$value <- NA
  new.obs <- rbind(new.obs,hv)
  hv$data_type_text <- "SELDOM"
  hv$value <- NA
  new.obs <- rbind(new.obs,hv)
  rm(hv)
  
  # Extreme cold
  ec <- droplevels(blsC[which(blsC$data_element_text=="Extreme cold" &
                                blsC$data_type_text=="NOT PRESENT"),])
  ec$data_type_text <- "CONSTANTLY"
  ec$value <- NA
  new.obs <- rbind(new.obs,ec)
  ec$data_type_text <- "FREQUENTLY"
  ec$value <- NA
  new.obs <- rbind(new.obs,ec)
  rm(ec)
  
  # Climbing ladders, ropes, or scaffolds
  cl <- droplevels(blsC[which(blsC$data_element_text=="Climbing ladders, ropes, or scaffolds" &
                                blsC$data_type_text=="NOT PRESENT"),])
  cl$data_type_text <- "CONSTANTLY"
  cl$value <- NA
  new.obs <- rbind(new.obs,cl)
  cl$data_type_text <- "FREQUENTLY"
  cl$value <- NA
  new.obs <- rbind(new.obs,cl)
  rm(cl)
  
  # Climbing ramps or stairs (work-related)
  cr <- droplevels(blsC[which(blsC$data_element_text=="Climbing ramps or stairs (work-related)" &
                                blsC$data_type_text=="NOT PRESENT"),])
  cr$data_type_text <- "CONSTANTLY"
  cr$value <- NA
  new.obs <- rbind(new.obs,cr)
  cr$data_type_text <- "FREQUENTLY"
  cr$value <- NA
  new.obs <- rbind(new.obs,cr)
  rm(cr)
  
  # Stooping
  stooping <- droplevels(blsC[which(blsC$data_element_text=="Stooping" &
                                      blsC$data_type_text=="NOT PRESENT"),])
  stooping$data_type_text <- "CONSTANTLY"
  stooping$value <- NA
  new.obs <- rbind(new.obs,stooping)
  rm(stooping)
  
  # Kneeling
  kneeling <- droplevels(blsC[which(blsC$data_element_text=="Kneeling" &
                                      blsC$data_type_text=="NOT PRESENT"),])
  kneeling$data_type_text <- "CONSTANTLY"
  kneeling$value <- NA
  new.obs <- rbind(new.obs,kneeling)
  kneeling$data_type_text <- "FREQUENTLY"
  kneeling$value <- NA
  new.obs <- rbind(new.obs,kneeling)
  rm(kneeling)
  
  # Crawling
  crawling <- droplevels(blsC[which(blsC$data_element_text=="Crawling" &
                                      blsC$data_type_text=="NOT PRESENT"),])
  crawling$data_type_text <- "CONSTANTLY"
  crawling$value <- NA
  new.obs <- rbind(new.obs,crawling)
  crawling$data_type_text <- "FREQUENTLY"
  crawling$value <- NA
  new.obs <- rbind(new.obs,crawling)
  rm(crawling)
  
  
  # Append new observations to data
  blsC <- rbind(blsC,new.obs)
  blsC <- distinct(blsC)
  rm(new.obs)
  return(blsC)
}


#' Make corrections to data
#'
#' During our analysis process we noticed some corrections/adjustments that 
#' needed to be made:
#' 
#' (1) Requirement 'Hearing requirements: Pass a hearing test' was split between
#' 'Hearing requirements: Pass a hearing test' and
#' 'Hearing requirements: Pass a hearing Test'; these are standardized to the 
#' former.
#' 
#' (2) Remove observations with 'FULLY MITIGATED' in the data_type_text field.
#' 
#' (3) Levels for the synthetic additive group (requirement) 78, ('Sitting or 
#' standing/walking') are set as 'SIT' and 'STAND' under data_type_text.
#' 
#' (4) Ensure that data is formatted as a data frame.
#'
#' @param blsC Completed data (output of fillMissingObservations(), or
#' otherMissingObservations())
#' @return Data, updated with corrections
#' @export
dataCorrections <- function(blsC) {
  # Data corrections
  blsC$data_element_text <- gsub("hearing Test","hearing test",blsC$data_element_text)
  blsC <- blsC[-which(blsC$data_type_text=="FULLY MITIGATED"),]
  blsC[blsC$data_element_text=="Sitting","data_type_text"] <- "SIT"
  blsC[blsC$data_element_text=="Standing/walking","data_type_text"] <- "STAND"
  
  # Formatting
  blsC <- as.data.frame(blsC)
}


#' Get mean estimates
#'
#' Select only the mean estimates from the data.
#'
#' @param blsC Completed, corrected data (output of dataCorrections())
#' @return Mean estimates only
#' @export
getMeanEstimates <- function(blsC) {
  blsMnC <- blsC[blsC$estimate_type_text=="Estimate",c("occupation_text","upper_soc_code","data_element_text",
                                                       "data_type_text","additive_group","value")]
  blsMnC <- blsMnC[order(blsMnC$additive_group),][order(blsMnC$occupation_text),] # reorder data
  
  return(blsMnC)
}


#' Get standard error estimates
#'
#' Select only the standard error estimates from the data.
#'
#' @param blsC Completed, corrected data (output of dataCorrections())
#' @return Standard error estimates only
#' @export
getErrors <- function(blsC) {
  blsE <- blsC[blsC$estimate_type_text=="Standard error",c("occupation_text","upper_soc_code",
                                                           "data_element_text","data_type_text",
                                                           "additive_group","value")]
  blsE <- blsE[order(blsE$additive_group),][order(blsE$occupation_text),] # reorder data
  
  return(blsE)
}


#' Calculate mean estimates for N-1 groups
#'
#' Some of the occupational groups are only missing a single mean estimate.
#' E.g., binary occupational groups (with requirement levels 'YES'/'NO') that 
#' have one estimate available can be filled in based on the known data. 
#' Similarly, if an occupational group has 5 levels, and 4 of them have known
#' estimates, then the final one can be calculated. This is because the total
#' value of any occupational group must be 100%.
#'
#' @param blsMnC Mean estimates (output of getMeanEstimates())
#' @return N-1 group completed data
#' @export
fillNminusOneGroups <- function(blsMnC) {
  yes.no <- droplevels(blsMnC[blsMnC$data_type_text=="YES" | blsMnC$data_type_text=="NO",])
  
  present.not <- droplevels(blsMnC[((blsMnC$data_type_text=="PRESENT" | 
                                       blsMnC$data_type_text=="SELDOM" | 
                                       blsMnC$data_type_text=="OCCASIONALLY" | 
                                       blsMnC$data_type_text=="FREQUENTLY" | 
                                       blsMnC$data_type_text=="CONSTANTLY" | 
                                       blsMnC$data_type_text=="NOT PRESENT") & 
                                      !is.na(blsMnC$additive_group)),])
  
  yn.addGroups <- unique(yes.no$additive_group)
  yn.occupations <- unique(yes.no$occupation_text)
  
  pn.addGroups <- unique(present.not$additive_group)
  pn.occupations <- unique(present.not$occupation_text)
  
  rm(yes.no)
  rm(present.not)
  
  for (ag in yn.addGroups) {
    for (oc in yn.occupations) {
      est.sum <- sum(blsMnC[blsMnC$additive_group==ag & blsMnC$occupation_text==oc,"value"],na.rm=TRUE)
      
      if (est.sum>100) {   # takes care of those whose values sum to > 100 (usually bc of rounding error)
        next
      }
      
      if (est.sum==0) {   # takes care of those with NA's or 0's for value (for both YES/NO)
        next
      }
      
      if (est.sum==100) {
        to.replace <- which(blsMnC$additive_group==ag & blsMnC$occupation_text==oc & is.na(blsMnC$value))
        if (length(to.replace)>=1) {
          blsMnC[to.replace,"value"] <- 100 - est.sum   # this will be 0 
        }
      }
      
      if (est.sum<100) {
        to.replace <- which(blsMnC$additive_group==ag & blsMnC$occupation_text==oc & is.na(blsMnC$value))
        if(length(to.replace)==1) {
          blsMnC[to.replace,"value"] <- 100 - est.sum
        }
      }
    }
  }
  
  for (ag in pn.addGroups) {
    for (oc in pn.occupations) {
      est.sum <- sum(blsMnC[blsMnC$additive_group==ag & blsMnC$occupation_text==oc,"value"],na.rm=TRUE)
      
      if (est.sum>100) {   # takes care of those whose values sum to > 100 (usually bc of rounding error)
        next
      }
      
      if (est.sum==0) {   # takes care of those with NA's or 0's for value (for all levels)
        next
      }
      
      if (est.sum==100) {
        to.replace <- which(blsMnC$additive_group==ag & blsMnC$occupation_text==oc & is.na(blsMnC$value))
        if (length(to.replace)>=1) {
          blsMnC[to.replace,"value"] <- 100 - est.sum   # this will be 0
        }
      }
      
      if (est.sum<100) {
        to.replace <- which(blsMnC$additive_group==ag & blsMnC$occupation_text==oc & is.na(blsMnC$value))
        if(length(to.replace)==1) {
          blsMnC[to.replace,"value"] <- 100 - est.sum
        }
      }
    }
  }
  
  blsMnGC <- blsMnC
  rm(blsMnC)
  
  return(blsMnGC)
}


#' Assign errors to known mean estimates, and bounds for all estimates (known
#' and unknown)
#'
#' The errors are collected alongside the mean estimates, and can allow for
#' empirical bounding of the data during modeling. For records lacking a mean
#' estimate/error, we can still bound the missing estimate using the remaining
#' value in an occupational group. For example, if there are three records in an
#' occupational group, of which two estimates are known and sum to X%, the upper
#' bound of the missing estimate is 100%-X%, and the lower bound is 0%.
#'
#' @param blsE Errors associated with mean estimates (output of getErrors())
#' @param blsMnGC N-1 group completed data (output of fillNminusOneGroups())
#' @return Data with relevant errors and bounds assigned to each record
#' @export
errorsAndBounding <- function(blsE,blsMnGC) {
  # Recall from the footnotes that:
  # stderr == -5 ==> stderr < 0.05
  # stderr == -6 ==> stderr < 0.5
  # estimate == -1 ==> estimate < 0.5
  blsE[which(blsE$value==-5),"value"] <- 0.05
  blsE[which(blsE$value==-6),"value"] <- 0.50
  blsMnGC[which(blsMnGC$value==-1),"value"] <- 0
  
  # Format data for easier handling
  blsE.names <- gsub(" ","",paste(blsE$occupation_text,blsE$additive_group,blsE$data_type_text))
  rownames(blsE) <- blsE.names
  blsMnGC.names <- gsub(" ","",paste(blsMnGC$occupation_text,blsMnGC$additive_group,blsMnGC$data_type_text))
  rownames(blsMnGC) <- blsMnGC.names
  
  # Assign blsE
  blsMnGC$std.error <- NA
  blsMnGC[blsE.names[which(blsE.names %in% blsMnGC.names)],"std.error"] <- blsE[blsE.names[which(blsE.names %in% blsMnGC.names)],"value"]
  blsMnGC[!is.na(blsMnGC$value) & is.na(blsMnGC$std.error),"std.error"] <- 0  # set missing blsE to 0???
  
  # Assign residual percentages
  blsMnGC$pct.remaining <- NA
  occ <- as.character(unique(blsMnGC$occupation_text))
  adg <- unique(blsMnGC$additive_group)
  for (i in occ) {
    #print(i)
    for (j in adg) {
      #print(j)
      current.group <- blsMnGC[blsMnGC$occupation_text==i & blsMnGC$additive_group==j,]
      resid.pct <- 100 - sum(current.group$value,na.rm=TRUE)
      current.group.resid <- current.group[which(is.na(current.group$value)),]
      blsMnGC[rownames(current.group.resid),"pct.remaining"] <- resid.pct
    }
  }
  
  # Reformat values as proportions of 1
  blsMnGC$value <- blsMnGC$value/100
  blsMnGC$std.error <- blsMnGC$std.error/100
  blsMnGC$pct.remaining <- blsMnGC$pct.remaining/100
  
  # Upper and lower bounds
  blsMnGC$upper_bound <- NA
  blsMnGC$lower_bound <- NA
  ub1 <- apply(cbind(blsMnGC$value + blsMnGC$std.error,rep(1,nrow(blsMnGC))),1,min)
  lb1 <- apply(cbind(blsMnGC$value - blsMnGC$std.error,rep(0,nrow(blsMnGC))),1,max)
  ub2 <- blsMnGC$pct.remaining
  lb2 <- blsMnGC$pct.remaining - blsMnGC$pct.remaining    # keeps NAs for known values for coalescing step (below)
  blsMnGC$upper_bound <- coalesce(ub1,ub2)
  blsMnGC$lower_bound <- coalesce(lb1,lb2)
  rm(list=c("ub1","ub2","lb1","lb2"))
  
  return(blsMnGC)
}


#' Populate predictors based on data mapping
#'
#' (1) The data_element_text field was restructured into a new field, 
#' Requirement, and added to the data.
#' 
#' (2) The data_element_text and data_type_text fields have somewhat overlapping
#' roles in the original data.  Moreover their structure makes for a large 
#' increase in the dimensions of the data. We devised a way to map these
#' two fields to numeric fields (Frequency and Intensity), leading to a
#' significant reduction in dimension. This mapping is applied here, and the new
#' fields are added to the data.
#' 
#' (3) Also added to the data are fields for SOC2, SOC3, and SOC4 codes, as well
#' as a field for the relevant requirement category.
#' 
#' (4) Generate an indicator column that differentiates between known estimates
#' and missing estimates.
#'
#' @param blsMnGCe Data with errors and bounds (output of errorsAndBounding())
#' @param predictors.data A data frame with data/predictor mappings
#' @return Data augmented with relevant predictors
#' @export
getPredictors <- function(blsMnGCe,predictors.data) {
  # Format existing columns as necessary
  blsMnGCe$occupation_text <- as.factor(blsMnGCe$occupation_text)
  blsMnGCe$additive_group <- as.numeric(as.character(blsMnGCe$additive_group))
  blsMnGCe$data_type_text <- as.character(blsMnGCe$data_type_text)
  
  # New predictor columns
  blsMnGCe$upSOC2 <- as.factor(substr(blsMnGCe[,"upper_soc_code"],1,2))
  blsMnGCe$upSOC3 <- as.factor(substr(blsMnGCe[,"upper_soc_code"],1,3))
  blsMnGCe$upSOC4 <- as.factor(substr(blsMnGCe[,"upper_soc_code"],1,4))
  blsMnGCe$characteristic <- NA
  blsMnGCe$frequency <- NA
  blsMnGCe$intensity <- NA
  blsMnGCe$req.cat <- NA
  
  # Indicator column
  blsMnGCe$known.val <- as.numeric(!is.na(blsMnGCe$value))
  # Overwrite complete synthetic additive groups 11, 89, 90, 91 that don't sum to 1 (this is a known issue)
  for (adg in c("11","89","90","91")) {
    for (occ in levels(blsMnGCe$occupation_text)) {
      current.occ.group <- blsMnGCe[as.character(blsMnGCe$occupation_text)==occ
                                    & as.character(blsMnGCe$additive_group)==adg,]
      if(sum(current.occ.group$value,na.rm=TRUE)!=1 &
         sum(current.occ.group$known.val,na.rm=TRUE)==2) # these are all binary additive groups!
      {
        blsMnGCe[rownames(current.occ.group),"known.val"] <- 2
      }
    }
  }
  
  # Assign predictors  
  bls.addGroups <- unique(blsMnGCe$additive_group)
  for (i in c(1:length(bls.addGroups))) {
    dtt <- droplevels(predictors.data[predictors.data$additive_group==bls.addGroups[i],])$data_type_text
    for (j in c(1:length(dtt))) {
      if(nrow(blsMnGCe[which(blsMnGCe$additive_group==bls.addGroups[i] & blsMnGCe$data_type_text==dtt[j]),])==0){
        next
      }
      
      if(nrow(predictors.data[predictors.data$additive_group==bls.addGroups[i] & predictors.data$data_type_text==dtt[j],])==0){
        next
      }
      
      blsMnGCe[which(blsMnGCe$additive_group==bls.addGroups[i] & blsMnGCe$data_type_text==dtt[j]),"characteristic"] <-
        predictors.data[predictors.data$additive_group==bls.addGroups[i] & predictors.data$data_type_text==dtt[j],"data_element_text_new"]
      blsMnGCe[which(blsMnGCe$additive_group==bls.addGroups[i] & blsMnGCe$data_type_text==dtt[j]),"frequency"] <-
        predictors.data[predictors.data$additive_group==bls.addGroups[i] & predictors.data$data_type_text==dtt[j],"Frequency"]
      blsMnGCe[which(blsMnGCe$additive_group==bls.addGroups[i] & blsMnGCe$data_type_text==dtt[j]),"intensity"] <-
        predictors.data[predictors.data$additive_group==bls.addGroups[i] & predictors.data$data_type_text==dtt[j],"Intensity"]
      blsMnGCe[which(blsMnGCe$additive_group==bls.addGroups[i] & blsMnGCe$data_type_text==dtt[j]),"req.cat"] <-
        predictors.data[predictors.data$additive_group==bls.addGroups[i] & predictors.data$data_type_text==dtt[j],"Requirement_category"]
    }
  }
  
  # Format new predictor columns as necessary
  blsMnGCe$additive_group <- as.factor(blsMnGCe$additive_group)
  blsMnGCe$characteristic <- as.factor(blsMnGCe$characteristic)
  blsMnGCe$req.cat <- as.factor(blsMnGCe$req.cat)
  
  # Final processing
  bls.data <- blsMnGCe  
  rm(blsMnGCe)
  return(bls.data)
}


#' Assign default modeling weights
#' 
#' In the modeling stage, data points can be weighted differently. This function
#' assigns default weights to records based on the presence/absence of a mean
#' estimate.
#'
#' @param bls.data Data augmented with relevant predictors (output of 
#' getPredictors())
#' @return Data augmented with relevant predictors and default modeling weights
#' @export
setDefaultModelingWeights <- function(bls.data) {
  # Assign default modeling weights based on indicator column
  bls.data$weight <- bls.data$known.val
  bls.data[bls.data$weight==0,"weight"] <- 0.5
  bls.data[bls.data$weight==2,"weight"] <- 1
  
  return(bls.data)
}
