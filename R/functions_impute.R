# Functions for imputing missing estimates


#' Smart Guessing for missing estimates
#'
#' Function to perform "smart guessing" procedure to generate initial guess for 
#' imputation of missing estimates. This approach leverages some of the 
#' structure built into the data. Namely, occupations within a given 2- or 
#' 3-digit SOC group are similar in nature. Thus, it is reasonable to assume 
#' that the requirements for the constituent occupations also follow similar 
#' distributions. This function uses information from the other members of a 
#' given SOC group to produce an initial guess for a particular occupation, and 
#' is briefly described below.
#'
#' The following procedure is followed within each SOC group. Each requirement 
#' is searched for an occupation with the "best" distribution, i.e. the job with
#' the maximum number of known estimates. In cases where there are multiple such
#' jobs, their requirement distributions are averaged to arrive at a single best
#' distribution. Then, each job (within a given requirement) is compared to this
#'  best distribution, and falls into one of three cases:
#' 
#' (1) Overlap between current job and the best distribution
#' (2) No overlap between current job and the best distribution
#' (3) Current job has no associated estimates (subset of case 2, above)
#' 
#' In the first case, missing estimates are populated as follows. A scaling 
#' factor was first computed based on the overlapping observations in the 
#' current job and the best distribution. This scaling factor is then multiplied
#' by the sum total of the estimates associated with observations in the best 
#' distribution that did not have counterparts in the current job, yielding some
#'  value x. The value of x is then evenly distributed across all the 
#' observations that were missing estimates in the current job, but had 
#' estimates in the best distribution. Finally, the sum of all the values in the
#' current job (both known, and guessed) is subtracted from 1, and this 
#' remaining value is evenly distributed across any outstanding observations in
#' the current job.
#' 
#' In the second case, the missing estimates are simply populated with the naive
#' guess for their value. For example, if the known estimates in the current job
#' sum to 0.8, and there are two observations with missing estimates, each one 
#' is given a value of 0.2 / 2 = 0.1. 
#' 
#' In the third case, the observations in the current job whose counterparts 
#' have known estimates in the best distribution simply receive the value of the
#' counterpart's estimate. The remaining estimates are populated using the naive
#' approach described in case 2.
#' 
#' The above procedure is completed per requirement, per SOC group. All guesses
#' are then adjusted to adhere to boundary conditions on the data (all estimates
#' must be in the range [0,1], and the sum of all estimates within an 
#' occupational group must be <=1). Note that the modeling weights associated 
#' with guessed values are altered based on which of the three cases they fall 
#' into, with those falling in cases 1 and 3 receiving higher weights, and those
#' falling in case 2 receiving lower weights. These weights are used in the 
#' iterative modeling step.
#'
#' @param ors.data.sims Original data augmented with relevant predictors, i.e. 
#' all records, including both known and missing estimates, as well as simulated
#' data (output of getSimulations())
#' @param sim.no Which simulation to run smart guessing on
#' @param wt.low Model weight to assign to low-confidence smart guesses; default
#' is 0
#' @param wt.mid Model weight to assign to mid-confidence smart guesses; default
#' is 0.5
#' @param wt.high Model weight to assign to high-confidence smart guesses; 
#' default is 0.5
#' @param verbose Should messages be printed; default is FALSE (mute messages) 
#' @return Input data frame, with missing values filled in by smart guesses
#' @export
smartGuess <- function(ors.data.sims,sim.no,
                       wt.low=0,wt.mid=0.5,wt.high=0.5,
                       verbose=FALSE) {
  
  if (is.numeric(sim.no)) {
    sim.col <- paste("valsim",sim.no,sep="")
    ors.data.sims$value <- ors.data.sims[,sim.col]
  }
  ors.data.sims <- ors.data.sims[,-grep("valsim",colnames(ors.data.sims))]
  
  # Pre-populate known/missing values into prediction column
  ors.data.sims$prediction <- ors.data.sims$value
  
  # Set missing value weights to minimum
  ors.data.sims[which(ors.data.sims$known.val!=1),"weight"] <- wt.low
  
  # Column to flag predictions that are forced scaled/bounded on [0,1]
  ors.data.sims$flag <- 0
  
  # Run smart guessing procedure
  for (adg in levels(ors.data.sims$additive_group)) {
    best.pct <- -1
    best.occ <- ""
    best.dist <- NA
    
    # Identify best distribution(s)
    if (verbose) {
      cat("\n",paste("Identifying most complete distribution for additive group: ",adg,sep=""),"\n")
    }
    for (occ in levels(ors.data.sims$occupation_text)) {
      current.occ.group <- ors.data.sims[as.character(ors.data.sims$occupation_text)==occ
                                         & as.character(ors.data.sims$additive_group)==adg,]
      known.pct <- sum(!is.na(current.occ.group$value)) / nrow(current.occ.group)
      
      if (nrow(current.occ.group)==0) {
        next
      }
      
      if (known.pct > best.pct[1]) {
        best.pct <- known.pct
        best.occ <- occ
        best.dist <- current.occ.group$value
        names(best.dist) <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
      }
      else if (known.pct==best.pct[1]) {
        best.pct <- c(best.pct,known.pct)
        best.occ <- c(best.occ,occ)
        best.dist <- NA
      }
    }
    
    # Average distribution (in case of equivalent best)
    if (length(best.occ) > 1) {
      if (verbose) {
        print("*** Averaging over multiple best distributions...")
      }
      bocc.tab <- data.frame(matrix(nrow=nrow(current.occ.group),ncol=length(best.occ)),
                             row.names=paste(current.occ.group$frequency,current.occ.group$intensity,sep="-"))
      colnames(bocc.tab) <- best.occ
      
      for (bocc in best.occ) {
        current.occ.group <- ors.data.sims[as.character(ors.data.sims$occupation_text)==bocc
                                           & as.character(ors.data.sims$additive_group)==adg,]
        current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
        bocc.tab[current.occ.group$to.sort,bocc] <- current.occ.group$value
      }
      best.dist <- rowMeans(bocc.tab,na.rm=TRUE)
      names(best.dist) <- rownames(bocc.tab)
    }
    best.dist[is.nan(best.dist)] <- NA
    if(sum(best.dist,na.rm=TRUE) > 1) {
      best.dist <- best.dist / sum(best.dist,na.rm=TRUE)
    }
    if (verbose) {
      print(best.dist)
    }
    
    # Assign upper bounds and guesses by occupational group (OG)
    for (occ in levels(ors.data.sims$occupation_text)) {
      current.occ.group <- ors.data.sims[as.character(ors.data.sims$occupation_text)==occ
                                         & as.character(ors.data.sims$additive_group)==adg,]
      current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
      current.missing <- rownames(current.occ.group)[which(is.na(current.occ.group$value))]
      
      # Calculate remaining percentage
      rem.pct <- (1 - sum(current.occ.group$value,na.rm=TRUE))
      if(rem.pct < 0) {
        rem.pct <- 0
      }
      
      # Reassign bounds (and pct.remaining) - not necessary, but retained for book-keeping
      ors.data.sims[current.missing,"upper_bound"] <- rem.pct
      ors.data.sims[current.missing,"lower_bound"] <- 0
      ors.data.sims[current.missing,"pct.remaining"] <- rem.pct
      
      # Check that the current OG has missing values
      if (length(current.missing)==0) {
        next
      }
      
      # Check that there is remaining percentage to assign
      if (rem.pct==0) {
        ors.data.sims[current.missing,"prediction"] <- 0
        next
      }
      
      # If there are values in the best distribution...
      if (sum(!is.na(best.dist)) >= 1) {
        
        # Determine which (if any) levels overlap between best distribution and current OG
        best.known <- names(best.dist)[which(!is.na(best.dist))]
        current.known <- current.occ.group$to.sort[which(!is.na(current.occ.group$value))]
        common.known <- intersect(best.known,current.known)
        
        # If there is overlap between best distribution and current OG..
        if (length(common.known) > 0) {  # make sure there is overlap
          if (verbose) {
            print(paste("Smart guessing for occupation: ",occ,sep=""))
          }
          
          # Common values
          c1 <- sum(best.dist[common.known],na.rm=TRUE)
          nc1 <- sum(best.dist,na.rm=TRUE) - c1
          c2 <- sum(current.occ.group[which(current.occ.group$to.sort %in% common.known),"value"],na.rm=TRUE)
          nc2 <- sum(current.occ.group[,"value"],na.rm=TRUE) - c2
          scaling.factor <- (1-c2) / (1-c1)
          
          # Levels of relevant missing values (those that violate boundary conditions)
          nj2.levels <- rownames(current.occ.group)[which(current.occ.group$to.sort %in% best.known & !(current.occ.group$to.sort %in% current.known))]
          jm2.levels <- rownames(current.occ.group)[which(!(current.occ.group$to.sort %in% best.known | current.occ.group$to.sort %in% current.known))]
          
          # Missing values
          if (is.na(scaling.factor) | scaling.factor==Inf) {
            nj2 <- 1 - sum(c2,nc2)
            ors.data.sims[c(nj2.levels,jm2.levels),"flag"] <- 1
            ors.data.sims[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else if ((nc1 * scaling.factor) > (1 - sum(c2,nc2))) {
            nj2 <- 1 - sum(c2,nc2)
            ors.data.sims[c(nj2.levels,jm2.levels),"flag"] <- 1
            ors.data.sims[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else {
            nj2 <- nc1 * scaling.factor
            ors.data.sims[c(nj2.levels,jm2.levels),"weight"] <- wt.high
          }
          jm2 <- 1 - sum(c2,nc2,nj2)
          
          # Populate smart guesses for relevant observations
          if (length(nj2.levels) > 0 & length(jm2.levels) > 0) {
            ors.data.sims[nj2.levels,"prediction"] <- nj2 / length(nj2.levels)
            ors.data.sims[jm2.levels,"prediction"] <- jm2 / length(jm2.levels)
          }
          else if ((length(nj2.levels) > 0 & length(jm2.levels)==0)) {
            ors.data.sims[nj2.levels,"prediction"] <- rem.pct / length(nj2.levels)
          }
          else if ((length(nj2.levels)==0 & length(jm2.levels) > 0)) {
            ors.data.sims[jm2.levels,"prediction"] <- rem.pct / length(jm2.levels)
          }
          else {
            next
          }
          
          # Scale guesses to [0,1]
          if (sum(ors.data.sims[rownames(current.occ.group),"prediction"])!=1) {
            ors.data.sims[current.missing,"prediction"] <- rem.pct * (ors.data.sims[current.missing,"prediction"] / sum(ors.data.sims[current.missing,"prediction"]))
          }
        }
        
        # If the current OG is an empty distribution (all NAs, i.e. special case of no overlap)...
        else if (sum(is.na(current.occ.group$value))==nrow(current.occ.group)) {  # make sure the OG is empty
          if (verbose) {
            print(paste("Smart guessing for occupation: ",occ,sep=""))
          }
          
          # Get observations that can take smart guesses/bounds (all)
          smart.guess.labels <- current.occ.group[which(is.na(current.occ.group$value)),"to.sort"]
          # Populate values from best distribution as smart guesses for current OG
          ors.data.sims[current.missing,"prediction"] <- best.dist[smart.guess.labels]
          ors.data.sims[current.missing,"weight"] <- wt.high
          
          # Generate naive guesses for the others, if necessary (i.e. best distribution has some NAs)
          if (sum(is.na(best.dist)) >= 1) {
            naive.guess <- rownames(current.occ.group)[which(current.occ.group$to.sort %in% names(best.dist)[which(is.na(best.dist))])]
            ors.data.sims[naive.guess,"prediction"] <- 
              (1 - sum(ors.data.sims[rownames(current.occ.group),"prediction"],na.rm=TRUE)) / sum(is.na(ors.data.sims[rownames(current.occ.group),"prediction"]))
          }
        }
        
        # Otherwise (non-empty OG, but no overlap)...
        else {
          # Generate naive prediction
          ors.data.sims[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group$value))
          ors.data.sims[current.missing,"weight"] <- wt.low
        }
      }
      # If the best distribution is empty...
      else {
        # Generate naive prediction
        ors.data.sims[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group$value))
        ors.data.sims[current.missing,"weight"] <- wt.low
      }
    }
  }
  
  return(ors.data.sims)
}



#' Impute missing estimates, with Smart Guessing
#'
#' Function to impute missing estimates in an iterative fashion. The first 
#' prediction (Iteration 0) is the result of the smart guessing procedure. All
#' subsequent iterations rely on XGBoost to produce preliminary predictions. 
#' These predictions are then adjusted to adhere to boundary conditions on the 
#' data (all estimates must be in the range [0,1], and the sum of all estimates
#' within an occupational group must be <=1).
#'
#' @param ors.data.sims Original data augmented with relevant predictors, i.e. 
#' all records, including both known and missing estimates, as well as simulated
#' data (output of getSimulations())
#' @param n.iter Number of times to iterate/adjust the model
#' @param weight.step Increment by which to increase modeling weight of test 
#' fold data with each iteration
#' @param mdl.d Tree model maximum depth; default is 14
#' @param mdl.n Tree model rounds; default is 200
#' @param mdl.e Tree model eta; default is 0.6
#' @param sg.soc.code SOC code to use for smart guessing, either "upSOC2", 
#' "upSOC3", or "upSOC4"
#' @return A list containing the results of iterative modeling, for each 
#' simulation
#' @export
iterateModel <- function(ors.data.sims,n.iter,weight.step,
                         mdl.d=14,mdl.n=200,mdl.e=0.6,
                         sg.soc.code) {
  
  # Separate observations by SOC2/3/4 code
  ors.data.sims.split <- list()
  for (i in levels(ors.data.sims[,sg.soc.code])) {
    ors.data.sims.split[[i]] <- droplevels(ors.data.sims[as.character(ors.data.sims[,sg.soc.code])==i,])
  }
  n.sims <- length(grep("valsim",colnames(ors.data.sims)))
  
  # Smart guess missing estimates (for each simulation)
  sim.list <- list()
  for (i in c(1:n.sims)) {
    print(paste("Smart guessing for simulation ",i," (Iteration 0)...",sep=""))
    sim.list[[i]] <- do.call(rbind,lapply(ors.data.sims.split,smartGuess,sim.no=i))
    rownames(sim.list[[i]]) <- gsub("^[0-9][0-9]*\\.","",rownames(sim.list[[i]]))
  }
  
  # Run iterative modeling process for each simulation
  print("-------ITERATING MODELS-------")
  doParallel::registerDoParallel(parallel::makeCluster(n.sims)) # use n.sims threads for parallel processing
  model.results <- foreach::foreach(current.sim=sim.list,sim.num=icount(),
                                    .packages=c('Matrix','xgboost','ModelMetrics')) %dopar% {
    
    print(paste("Running simulation ",sim.num,":",sep=""))
    
    # Initialize list
    model.iterations <- list()
    model.iterations[[1]] <- list(data=current.sim)
    
    # Get column names of predictors/response
    select.cols <- c("occupation_text","requirement","frequency","intensity","req.cat","upSOC2","upSOC3","prediction")
    
    # Iterate
    for (iter in c(1:n.iter)) {
      print(paste("Iteration ",iter,"...",sep=""))
      
      # Get data/predictions from previous fit
      data.train <- model.iterations[[iter]]$data
      
      # Update weights of guesses to reflect new guess (assigned below)
      if (iter > 1) {
        data.train[data.train$weight<1,"weight"] <- data.train[data.train$weight<1,"weight"] + weight.step
        data.train[which(data.train$weight<1 & data.train$weight>=0.75),"weight"] <- 0.75   # cap at 0.75
        data.train[data.train$weight==weight.step,"weight"] <- 0    # naive guesses get 0
      }
      
      # Format data and assign bounds
      sparse.train <- xgboost::xgb.DMatrix(data=sparse.model.matrix(prediction ~ .,data=data.train[,select.cols])[,-1],label=data.train$prediction)
      xgboost::setinfo(sparse.train,'label_lower_bound',data.train$lower_bound)
      xgboost::setinfo(sparse.train,'label_upper_bound',data.train$upper_bound)
      xgboost::setinfo(sparse.train,'weight',data.train$weight)
      
      # Fit the model and make predictions
      # Note: https://github.com/uber/causalml/issues/96
      print(paste("Evaluating model... "))
      mdl <- xgboost::xgboost(data=sparse.train,booster="gbtree",verbose=0,max_depth=mdl.d,nrounds=mdl.n,
                              eta=mdl.e,nthread=20,objective="reg:squarederror")
      data.train$current.prediction <- predict(mdl,newdata=sparse.train)
      data.train$prediction.raw.bounded <- data.train$current.prediction
      data.train$prediction.raw <- data.train$current.prediction
      
      # Correct/scale new predictions based on known information (bounds, known values, percentage rules, etc.)
      print("Adjusting predictions...")
      # Bounds
      data.train[which(data.train$current.prediction < 0),"current.prediction"] <- 0
      data.train[which(data.train$current.prediction > 1),"current.prediction"] <- 1
      data.train$prediction.raw.bounded <- data.train$current.prediction
      # Known values
      data.train[which(data.train$known.val==1),"current.prediction"] <- data.train[which(data.train$known.val==1),"value"]
      # Percentage rules
      to.adjust <- rownames(data.train)[which(data.train$known.val!=1)]
      for (occ in levels(data.train$occupation_text)) {
        for (adg in levels(data.train$additive_group)) {
          current.full <- data.train[which(as.character(data.train$occupation_text)==occ & 
                                             as.character(data.train$additive_group)==adg),]
          current.to.adjust <- to.adjust[to.adjust %in% rownames(current.full)]
          
          # Additive groups not summing to 1
          if(nrow(current.full) >= 1 &
             length(current.to.adjust) > 0 & 
             sum(current.full$current.prediction) != 1 &
             sum(current.full[current.to.adjust,"current.prediction"]) != 0) {
            
            data.train[current.to.adjust,"current.prediction"] <- 
              (data.train[current.to.adjust,"current.prediction"] / sum(data.train[current.to.adjust,"current.prediction"])) * mean(data.train[current.to.adjust,"upper_bound"],na.rm=TRUE)
          }
          
          # Check for additive groups that sum to <1 and have all missing estimates with (current) prediction of 0
          # This will happen due to bounding adjustments (and above <if> statement will change nothing)
          # For these observations, assign original smart guess
          if (nrow(current.full) >= 1 &
              length(current.to.adjust) > 0 & 
              sum(current.full$current.prediction) < 1 & 
              # sum(current.full[current.to.adjust,"current.prediction"]==0,na.rm=TRUE)>0 &
              sum(current.full[current.to.adjust,"current.prediction"]) == 0) {
            data.train[current.to.adjust,"current.prediction"] <- model.iterations[[1]]$data[current.to.adjust,"prediction"]
          }
        }
      }
      
      # Format results
      res <- data.frame(previous.prediction=data.train$prediction,
                        current.prediction=data.train$current.prediction,
                        row.names=rownames(data.train))
      data.train$prediction <- data.train$current.prediction
      data.train$current.prediction <- NULL
      
      # Get predictions, model, data, MAE, and ME
      iter.results <- list(results=res,model=mdl,data=data.train,
                           missing.mae=ModelMetrics::mae(res[to.adjust,"previous.prediction"],res[to.adjust,"current.prediction"]),
                           missing.me=mean(res[to.adjust,"previous.prediction"]-res[to.adjust,"current.prediction"]),
                           known.mae=ModelMetrics::mae(data.train[which(data.train$known.val==1),"value"],
                                                       data.train[which(data.train$known.val==1),"prediction.raw.bounded"]),
                           known.me=mean(data.train[which(data.train$known.val==1),"value"]-
                                           data.train[which(data.train$known.val==1),"prediction.raw.bounded"]))
      
      # Store results
      model.iterations[[iter+1]] <- iter.results
    }
    
    names(model.iterations) <- paste("iteration",c(0:n.iter),sep="")
    list(model.iterations)[[1]]
  }
  
  names(model.results) <- paste("simulation",c(1:length(model.results)),sep="")
  return(model.results)
  
}


#' Get predictions (one per simulation) and overall mean prediction for a given
#' iteration
#' 
#' DETAILED DESCRIPTION
#'
#' @param simulations 
#' @param iteration 
#' @param ors.data.sims 
#' @return 
#' @export
computeMeanPredictions <- function(simulations,iteration,ors.data.sims) {
  # Compute convergence iterations and get relevant predictions
  # Uses population means - mean of simulation predictions for missing observations, actual values for known
  n.sim <- length(simulations)
  data <- data.frame(matrix(nrow=nrow(simulations[[1]]$iteration0$data),ncol=n.sim),
                     row.names=rownames(simulations[[1]]$iteration0$data))
  colnames(data) <- c(paste("sim",c(1:n.sim),".prediction",sep=""))
  data <- cbind(simulations[[1]]$iteration0$data,data)
  data <- data[,-grep("^prediction",colnames(data))]
  data <- data[,-grep("value",colnames(data))]
  data <- data[,-grep("weight",colnames(data))]
  data <- data[,-grep("flag",colnames(data))]
  
  if (is.numeric(iteration)) {
    print(paste("Calculating mean prediction at iteration",iteration))
    for (i in c(1:n.sim)) {
      data[,paste("sim",i,".prediction",sep="")] <- simulations[[i]][[iteration+1]]$data$prediction
    }
  }
  else {
    print("Unexpected input for 'iteration'")
    return()
  }
  
  data$mean.prediction <- rowMeans(data[,paste("sim",c(1:n.sim),".prediction",sep="")])
  data[data$known.val==1,"mean.prediction"] <- ors.data.sims[rownames(data)[data$known.val==1],"value"]
  
  return(data) 
}


#' Blend model results according to calculated proportions
#'
#' DETAILED DESCRIPTION
#'
#' @param model.results.soc2 
#' @param model.results.soc3 
#' @param conv.iter.soc2 
#' @param conv.iter.soc3 
#' @param soc2.prop 
#' @param soc3.prop 
#' @param ors.data.sims 
#' @param write.files 
#' @return 
#' @export
blendImputations <- function(model.results.soc2,model.results.soc3,
                             conv.iter.soc2,conv.iter.soc3,soc2.prop,soc3.prop,
                             write.files=FALSE) {
  
  # Format results from model 1
  model.summary.soc2 <- imputeORS::computeMeanPredictions(model.results.soc2,conv.iter.soc2,ors.data.sims)
  model.summary.soc2$data_element_text <- gsub("Lifting/carrying","Lifting/carrying:",model.summary.soc2$data_element_text)
  model.summary.soc2$data_element_text <- gsub("Communicating verbally","Speaking",model.summary.soc2$data_element_text)
  model.summary.soc2$data_element_text <- as.factor(model.summary.soc2$data_element_text)

  # Format results from model 2
  model.summary.soc3 <- imputeORS::computeMeanPredictions(model.results.soc3,conv.iter.soc3,ors.data.sims)
  model.summary.soc3$data_element_text <- gsub("Lifting/carrying","Lifting/carrying:",model.summary.soc3$data_element_text)
  model.summary.soc3$data_element_text <- gsub("Communicating verbally","Speaking",model.summary.soc3$data_element_text)
  model.summary.soc3$data_element_text <- as.factor(model.summary.soc3$data_element_text)
  model.summary.soc3 <- model.summary.soc3[rownames(model.summary.soc2),]
  
  # Blend
  select.cols <- c(paste("sim",c(1:length(model.results.soc2)),".prediction",sep=""),"mean.prediction")
  blended.results <- model.summary.soc2
  blended.results[,select.cols] <- (soc2.prop * model.summary.soc2[,select.cols]) + (soc3.prop * model.summary.soc3[,select.cols])
  
  # Write results
  # blended.results$actual <- ors.data.sims[rownames(blended.results),"value"]
  if (write.files) {
    write.csv(blended.results,file=paste(soc2.prop,"soc2-",soc3.prop,"soc3",".csv",sep=""))
    save(blended.results,file=paste(soc2.prop,"soc2-",soc3.prop,"soc3",".Rdata",sep="")) 
  }
  
  return(blended.results)
}


#' Approximate uncertainty due to models
#' 
#' DETAILED DESCRIPTION
#' 
#' @param model.results.soc2 
#' @param model.results.soc3 
#' @param conv.iter.soc2 
#' @param conv.iter.soc3 
#' @param soc2.prop 
#' @param soc3.prop 
#' @return 
#' @export
computeModelUncertainty <- function(model.results.soc2,model.results.soc3,
                                    conv.iter.soc2,conv.iter.soc3,soc2.prop,soc3.prop) {
  mdl.uncert.mae <-vector()
  mdl.uncert.me <-vector()
  for (i in c(1:length(model.results.soc2))) {
    # print(i) # simulation number
    
    # Get data for known values
    a <- model.results.soc2[[i]][[conv.iter.soc2+1]]$data
    a <- a[a$known.val==1,]
    b <- model.results.soc3[[i]][[conv.iter.soc3+1]]$data
    b <- b[b$known.val==1,]
    
    # Get bounded predictions only
    c <- a$value
    b <- b[rownames(a),"prediction.raw.bounded"]
    a <- a$prediction.raw.bounded
    
    # Calculations
    ab <- (soc2.prop * a) + (soc3.prop * b)
    mdl.uncert.mae <- c(mdl.uncert.mae,ModelMetrics::mae(c,ab))
    mdl.uncert.me <- c(mdl.uncert.me,mean(c-ab))
  }
  
  return(list(mae.by.sim=mdl.uncert.mae,me.by.sim=mdl.uncert.me,
              avg.mae=mean(mdl.uncert.mae),avg.me=mean(mdl.uncert.me)))
}
