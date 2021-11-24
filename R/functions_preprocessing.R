# Preprocessing functions

#' Get raw ORS data from BLS
#'
#' Retrieve raw Occupational Requirements Survey (ORS) data from the Bureau of 
#' Labor Statistics' (BLS) website.
#'
#' @param orsLink pointing to excel data file at BLS;
#' default is https://www.bls.gov/ors/xlsx/2018_excel_output.xlsx
#' @return Raw ORS data from BLS website
#' @export
getORS <- function(orsLink="https://www.bls.gov/ors/xlsx/2018_excel_output.xlsx") {
  url1 <- orsLink
  p1f <- tempfile()
  download.file(url1,p1f,mode="wb")
  orsRaw <- readxl::read_excel(path=p1f,sheet=2)
  return(orsRaw)
}


#' Transform raw ORS data
#'
#' The data retrieved from BLS website needs to be transformed into the format 
#' expected by all subsequent preprocessing (and other) functions in this 
#' package. This function accomplishes this task.
#'
#' @param 
#' @return 
#' @export
transformORS <- function() {

}


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
#' @param orsRaw Original data (read in from CSV file)
#' @return Data edited for new additive groups, and containing only relevant
#' records (mean and standard error estimates) and columns
#' @export
syntheticAddGroups <- function(orsRaw) {
  
  orsRaw[which(orsRaw$data_element_text=="Literacy required" &
                 (orsRaw$data_type_text=="YES" | orsRaw$data_type_text=="NO") &
                 orsRaw$estimate_type=="Percent"),"additive_group"] <- 11
  orsRaw[which(orsRaw$data_element_text=="Pre-employment training: Certification" &
                 (orsRaw$data_type_text=="YES" | orsRaw$data_type_text=="NO") &
                 orsRaw$estimate_type=="Percent"),"additive_group"] <- 89
  orsRaw[which(orsRaw$data_element_text=="Pre-employment training: License" &
                 (orsRaw$data_type_text=="YES" | orsRaw$data_type_text=="NO") &
                 orsRaw$estimate_type=="Percent"),"additive_group"] <- 90
  orsRaw[which(orsRaw$data_element_text=="Pre-employment training: Educational Certificate" &
                 (orsRaw$data_type_text=="YES" | orsRaw$data_type_text=="NO") &
                 orsRaw$estimate_type=="Percent"),"additive_group"] <- 91
  
  # In original analysis we used '100' for the new additive group
  # '100' is actually already taken by 'Pace: Pause Control' (a COG requirement)
  # Instead for the package we use '78', i.e. the reference group of 'Sitting'
  orsRaw[(orsRaw$data_element_text=="Sitting" | orsRaw$data_element_text=="Standing/walking") &
           orsRaw$unit_of_measure=="Percent of day" & orsRaw$estimate_type=="Mean","additive_group"] <- 78
  
  # Select only required observations and columns
  ors <- orsRaw[!is.na(orsRaw$additive_group),c("occupation_text","upper_soc_code","data_element_text",
                                                "data_type_text","estimate_type_text","additive_group","value")]
  
  return(ors)
}


#' Generate missing observations
#'
#' Observations that were missing in the data (i.e., no associated record or 
#' estimate) are generated here.
#' This is done by "completing" the data using tidyr::complete.
#'
#' @param ors Data frame of relevant observations (output of syntheticAddGroups())
#' @return Data with newly generated observations added (completed data)
#' @export
fillMissingObservations <- function(ors) {
  orsC = ors%>% tidyr::complete(occupation_text,nesting(data_type_text,additive_group,
                                                        estimate_type_text,data_element_text))
  
  for (occ in unique(orsC$occupation_text)) {    # fill in upper SOC code for new observations
    usoc <- as.numeric(names(table(orsC[orsC$occupation_text==occ,"upper_soc_code"])))
    if (length(usoc)==1) {
      orsC[orsC$occupation_text==occ,"upper_soc_code"] <- usoc
    }
  }
  
  return(orsC)
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
#' @param orsC Completed data (output of fillMissingObservations())
#' @return Completed data, with requested requirement/levels added
#' @export
otherMissingObservations <- function(orsC) {
  # Humidity
  hu <- droplevels(orsC[which(orsC$data_element_text=="Humidity" &
                                orsC$data_type_text=="NOT PRESENT"),])
  hu$data_type_text <- "CONSTANTLY"
  hu$value <- NA
  new.obs <- hu
  rm(hu)
  
  # Heavy vibrations
  hv <- droplevels(orsC[which(orsC$data_element_text=="Heavy vibrations" &
                                orsC$data_type_text=="NOT PRESENT"),])
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
  ec <- droplevels(orsC[which(orsC$data_element_text=="Extreme cold" &
                                orsC$data_type_text=="NOT PRESENT"),])
  ec$data_type_text <- "CONSTANTLY"
  ec$value <- NA
  new.obs <- rbind(new.obs,ec)
  ec$data_type_text <- "FREQUENTLY"
  ec$value <- NA
  new.obs <- rbind(new.obs,ec)
  rm(ec)
  
  # Climbing ladders, ropes, or scaffolds
  cl <- droplevels(orsC[which(orsC$data_element_text=="Climbing ladders, ropes, or scaffolds" &
                                orsC$data_type_text=="NOT PRESENT"),])
  cl$data_type_text <- "CONSTANTLY"
  cl$value <- NA
  new.obs <- rbind(new.obs,cl)
  cl$data_type_text <- "FREQUENTLY"
  cl$value <- NA
  new.obs <- rbind(new.obs,cl)
  rm(cl)
  
  # Climbing ramps or stairs (work-related)
  cr <- droplevels(orsC[which(orsC$data_element_text=="Climbing ramps or stairs (work-related)" &
                                orsC$data_type_text=="NOT PRESENT"),])
  cr$data_type_text <- "CONSTANTLY"
  cr$value <- NA
  new.obs <- rbind(new.obs,cr)
  cr$data_type_text <- "FREQUENTLY"
  cr$value <- NA
  new.obs <- rbind(new.obs,cr)
  rm(cr)
  
  # Stooping
  stooping <- droplevels(orsC[which(orsC$data_element_text=="Stooping" &
                                      orsC$data_type_text=="NOT PRESENT"),])
  stooping$data_type_text <- "CONSTANTLY"
  stooping$value <- NA
  new.obs <- rbind(new.obs,stooping)
  rm(stooping)
  
  # Kneeling
  kneeling <- droplevels(orsC[which(orsC$data_element_text=="Kneeling" &
                                      orsC$data_type_text=="NOT PRESENT"),])
  kneeling$data_type_text <- "CONSTANTLY"
  kneeling$value <- NA
  new.obs <- rbind(new.obs,kneeling)
  kneeling$data_type_text <- "FREQUENTLY"
  kneeling$value <- NA
  new.obs <- rbind(new.obs,kneeling)
  rm(kneeling)
  
  # Crawling
  crawling <- droplevels(orsC[which(orsC$data_element_text=="Crawling" &
                                      orsC$data_type_text=="NOT PRESENT"),])
  crawling$data_type_text <- "CONSTANTLY"
  crawling$value <- NA
  new.obs <- rbind(new.obs,crawling)
  crawling$data_type_text <- "FREQUENTLY"
  crawling$value <- NA
  new.obs <- rbind(new.obs,crawling)
  rm(crawling)
  
  
  # Append new observations to data
  orsC <- rbind(orsC,new.obs)
  orsC <- dplyr::distinct(orsC)
  rm(new.obs)
  return(orsC)
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
#' @param orsC Completed data (output of fillMissingObservations(), or
#' otherMissingObservations())
#' @return Data, updated with corrections
#' @export
dataCorrections <- function(orsC) {
  # Data corrections
  orsC$data_element_text <- gsub("hearing Test","hearing test",orsC$data_element_text)
  orsC <- orsC[-which(orsC$data_type_text=="FULLY MITIGATED"),]
  orsC[orsC$data_element_text=="Sitting","data_type_text"] <- "SIT"
  orsC[orsC$data_element_text=="Standing/walking","data_type_text"] <- "STAND"
  
  # Formatting
  orsC <- as.data.frame(orsC)
}


#' Get mean estimates
#'
#' Select only the mean estimates from the data.
#'
#' @param orsC Completed, corrected data (output of dataCorrections())
#' @return Mean estimates only
#' @export
getMeanEstimates <- function(orsC) {
  orsMnC <- orsC[orsC$estimate_type_text=="Estimate",c("occupation_text","upper_soc_code","data_element_text",
                                                       "data_type_text","additive_group","value")]
  orsMnC <- orsMnC[order(orsMnC$additive_group),][order(orsMnC$occupation_text),] # reorder data
  
  return(orsMnC)
}


#' Get standard error estimates
#'
#' Select only the standard error estimates from the data.
#'
#' @param orsC Completed, corrected data (output of dataCorrections())
#' @return Standard error estimates only
#' @export
getErrors <- function(orsC) {
  orsE <- orsC[orsC$estimate_type_text=="Standard error",c("occupation_text","upper_soc_code",
                                                           "data_element_text","data_type_text",
                                                           "additive_group","value")]
  orsE <- orsE[order(orsE$additive_group),][order(orsE$occupation_text),] # reorder data
  
  return(orsE)
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
#' @param orsMnC Mean estimates (output of getMeanEstimates())
#' @return N-1 group completed data
#' @export
fillNminusOneGroups <- function(orsMnC) {
  yes.no <- droplevels(orsMnC[orsMnC$data_type_text=="YES" | orsMnC$data_type_text=="NO",])
  
  present.not <- droplevels(orsMnC[((orsMnC$data_type_text=="PRESENT" | 
                                       orsMnC$data_type_text=="SELDOM" | 
                                       orsMnC$data_type_text=="OCCASIONALLY" | 
                                       orsMnC$data_type_text=="FREQUENTLY" | 
                                       orsMnC$data_type_text=="CONSTANTLY" | 
                                       orsMnC$data_type_text=="NOT PRESENT") & 
                                      !is.na(orsMnC$additive_group)),])
  
  yn.addGroups <- unique(yes.no$additive_group)
  yn.occupations <- unique(yes.no$occupation_text)
  
  pn.addGroups <- unique(present.not$additive_group)
  pn.occupations <- unique(present.not$occupation_text)
  
  rm(yes.no)
  rm(present.not)
  
  for (ag in yn.addGroups) {
    for (oc in yn.occupations) {
      est.sum <- sum(orsMnC[orsMnC$additive_group==ag & orsMnC$occupation_text==oc,"value"],na.rm=TRUE)
      
      if (est.sum>100) {   # takes care of those whose values sum to > 100 (usually bc of rounding error)
        next
      }
      
      if (est.sum==0) {   # takes care of those with NA's or 0's for value (for both YES/NO)
        next
      }
      
      if (est.sum==100) {
        to.replace <- which(orsMnC$additive_group==ag & orsMnC$occupation_text==oc & is.na(orsMnC$value))
        if (length(to.replace)>=1) {
          orsMnC[to.replace,"value"] <- 100 - est.sum   # this will be 0 
        }
      }
      
      if (est.sum<100) {
        to.replace <- which(orsMnC$additive_group==ag & orsMnC$occupation_text==oc & is.na(orsMnC$value))
        if(length(to.replace)==1) {
          orsMnC[to.replace,"value"] <- 100 - est.sum
        }
      }
    }
  }
  
  for (ag in pn.addGroups) {
    for (oc in pn.occupations) {
      est.sum <- sum(orsMnC[orsMnC$additive_group==ag & orsMnC$occupation_text==oc,"value"],na.rm=TRUE)
      
      if (est.sum>100) {   # takes care of those whose values sum to > 100 (usually bc of rounding error)
        next
      }
      
      if (est.sum==0) {   # takes care of those with NA's or 0's for value (for all levels)
        next
      }
      
      if (est.sum==100) {
        to.replace <- which(orsMnC$additive_group==ag & orsMnC$occupation_text==oc & is.na(orsMnC$value))
        if (length(to.replace)>=1) {
          orsMnC[to.replace,"value"] <- 100 - est.sum   # this will be 0
        }
      }
      
      if (est.sum<100) {
        to.replace <- which(orsMnC$additive_group==ag & orsMnC$occupation_text==oc & is.na(orsMnC$value))
        if(length(to.replace)==1) {
          orsMnC[to.replace,"value"] <- 100 - est.sum
        }
      }
    }
  }
  
  orsMnGC <- orsMnC
  rm(orsMnC)
  
  return(orsMnGC)
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
#' @param orsE Errors associated with mean estimates (output of getErrors())
#' @param orsMnGC N-1 group completed data (output of fillNminusOneGroups())
#' @return Data with relevant errors and bounds assigned to each record
#' @export
errorsAndBounding <- function(orsE,orsMnGC) {
  # Recall from the footnotes that:
  # stderr == -5 ==> stderr < 0.05
  # stderr == -6 ==> stderr < 0.5
  # estimate == -1 ==> estimate < 0.5
  orsE[which(orsE$value==-5),"value"] <- 0.05
  orsE[which(orsE$value==-6),"value"] <- 0.50
  orsMnGC[which(orsMnGC$value==-1),"value"] <- 0
  
  # Format data for easier handling
  orsE.names <- gsub(" ","",paste(orsE$occupation_text,orsE$additive_group,orsE$data_type_text))
  rownames(orsE) <- orsE.names
  orsMnGC.names <- gsub(" ","",paste(orsMnGC$occupation_text,orsMnGC$additive_group,orsMnGC$data_type_text))
  rownames(orsMnGC) <- orsMnGC.names
  
  # Assign orsE
  orsMnGC$std.error <- NA
  orsMnGC[orsE.names[which(orsE.names %in% orsMnGC.names)],"std.error"] <- orsE[orsE.names[which(orsE.names %in% orsMnGC.names)],"value"]
  orsMnGC[!is.na(orsMnGC$value) & is.na(orsMnGC$std.error),"std.error"] <- 0  # set missing orsE to 0???
  
  # Assign residual percentages
  orsMnGC$pct.remaining <- NA
  occ <- as.character(unique(orsMnGC$occupation_text))
  adg <- unique(orsMnGC$additive_group)
  for (i in occ) {
    #print(i)
    for (j in adg) {
      #print(j)
      current.group <- orsMnGC[orsMnGC$occupation_text==i & orsMnGC$additive_group==j,]
      resid.pct <- 100 - sum(current.group$value,na.rm=TRUE)
      current.group.resid <- current.group[which(is.na(current.group$value)),]
      orsMnGC[rownames(current.group.resid),"pct.remaining"] <- resid.pct
    }
  }
  
  # Reformat values as proportions of 1
  orsMnGC$value <- orsMnGC$value/100
  orsMnGC$std.error <- orsMnGC$std.error/100
  orsMnGC$pct.remaining <- orsMnGC$pct.remaining/100
  
  # Upper and lower bounds
  orsMnGC$upper_bound <- NA
  orsMnGC$lower_bound <- NA
  ub1 <- apply(cbind(orsMnGC$value + orsMnGC$std.error,rep(1,nrow(orsMnGC))),1,min)
  lb1 <- apply(cbind(orsMnGC$value - orsMnGC$std.error,rep(0,nrow(orsMnGC))),1,max)
  ub2 <- orsMnGC$pct.remaining
  lb2 <- orsMnGC$pct.remaining - orsMnGC$pct.remaining    # keeps NAs for known values for coalescing step (below)
  orsMnGC$upper_bound <- dplyr::coalesce(ub1,ub2)
  orsMnGC$lower_bound <- dplyr::coalesce(lb1,lb2)
  rm(list=c("ub1","ub2","lb1","lb2"))
  
  return(orsMnGC)
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
#' as a field for the relevant Requirement Category.
#' 
#' (4) Generate an indicator column that differentiates between known estimates
#' and missing estimates.
#'
#' @param orsMnGCe Data with errors and bounds (output of errorsAndBounding())
#' @return Data augmented with relevant predictors
#' @export
getPredictors <- function(orsMnGCe) {
  # Get mapping for data/predictors
  predictors.data <- imputeORS:::predictors.data
  
  # Format existing columns as necessary
  orsMnGCe$occupation_text <- as.factor(orsMnGCe$occupation_text)
  orsMnGCe$additive_group <- as.numeric(as.character(orsMnGCe$additive_group))
  orsMnGCe$data_type_text <- as.character(orsMnGCe$data_type_text)
  
  # New predictor columns
  orsMnGCe$upSOC2 <- as.factor(substr(orsMnGCe[,"upper_soc_code"],1,2))
  orsMnGCe$upSOC3 <- as.factor(substr(orsMnGCe[,"upper_soc_code"],1,3))
  orsMnGCe$upSOC4 <- as.factor(substr(orsMnGCe[,"upper_soc_code"],1,4))
  orsMnGCe$requirement <- NA
  orsMnGCe$frequency <- NA
  orsMnGCe$intensity <- NA
  orsMnGCe$req.cat <- NA
  
  # Indicator column
  orsMnGCe$known.val <- as.numeric(!is.na(orsMnGCe$value))
  # Overwrite complete synthetic additive groups 11, 89, 90, 91 that don't sum to 1 (this is a known issue)
  for (adg in c("11","89","90","91")) {
    for (occ in levels(orsMnGCe$occupation_text)) {
      current.occ.group <- orsMnGCe[as.character(orsMnGCe$occupation_text)==occ
                                    & as.character(orsMnGCe$additive_group)==adg,]
      if(sum(current.occ.group$value,na.rm=TRUE)!=1 &
         sum(current.occ.group$known.val,na.rm=TRUE)==2) # these are all binary additive groups!
      {
        orsMnGCe[rownames(current.occ.group),"known.val"] <- 2
      }
    }
  }
  
  # Assign predictors  
  ors.addGroups <- unique(orsMnGCe$additive_group)
  for (i in c(1:length(ors.addGroups))) {
    dtt <- droplevels(predictors.data[predictors.data$additive_group==ors.addGroups[i],])$data_type_text
    for (j in c(1:length(dtt))) {
      if(nrow(orsMnGCe[which(orsMnGCe$additive_group==ors.addGroups[i] & orsMnGCe$data_type_text==dtt[j]),])==0){
        next
      }
      
      if(nrow(predictors.data[predictors.data$additive_group==ors.addGroups[i] & predictors.data$data_type_text==dtt[j],])==0){
        next
      }
      
      orsMnGCe[which(orsMnGCe$additive_group==ors.addGroups[i] & orsMnGCe$data_type_text==dtt[j]),"requirement"] <-
        predictors.data[predictors.data$additive_group==ors.addGroups[i] & predictors.data$data_type_text==dtt[j],"Requirement"]
      orsMnGCe[which(orsMnGCe$additive_group==ors.addGroups[i] & orsMnGCe$data_type_text==dtt[j]),"frequency"] <-
        predictors.data[predictors.data$additive_group==ors.addGroups[i] & predictors.data$data_type_text==dtt[j],"Frequency"]
      orsMnGCe[which(orsMnGCe$additive_group==ors.addGroups[i] & orsMnGCe$data_type_text==dtt[j]),"intensity"] <-
        predictors.data[predictors.data$additive_group==ors.addGroups[i] & predictors.data$data_type_text==dtt[j],"Intensity"]
      orsMnGCe[which(orsMnGCe$additive_group==ors.addGroups[i] & orsMnGCe$data_type_text==dtt[j]),"req.cat"] <-
        predictors.data[predictors.data$additive_group==ors.addGroups[i] & predictors.data$data_type_text==dtt[j],"Requirement_category"]
    }
  }
  
  # Format new predictor columns as necessary
  orsMnGCe$additive_group <- as.factor(orsMnGCe$additive_group)
  orsMnGCe$requirement <- as.factor(orsMnGCe$requirement)
  orsMnGCe$req.cat <- as.factor(orsMnGCe$req.cat)
  
  # Final processing
  ors.data <- orsMnGCe  
  rm(orsMnGCe)
  return(ors.data)
}


#' Assign default modeling weights
#' 
#' In the modeling stage, data points can be weighted differently. This function
#' assigns default weights to records based on the presence/absence of a mean
#' estimate.
#'
#' @param ors.data Data augmented with relevant predictors (output of 
#' getPredictors())
#' @param unknown.weight Modeling weight to assign to unknown estimates; default
#' is 0.5
#' @return Data augmented with relevant predictors and default modeling weights
#' @export
setDefaultModelingWeights <- function(ors.data,unknown.weight=0.5) {
  # Assign default modeling weights based on indicator column
  ors.data$weight <- ors.data$known.val
  ors.data[ors.data$weight==0,"weight"] <- unknown.weight
  ors.data[ors.data$weight==2,"weight"] <- 1
  
  return(ors.data)
}

  