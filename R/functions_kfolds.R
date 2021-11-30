# Functions for k-folds cross validation

#' Smart Guessing for k-folds CV
#'
#' Function to perform "smart guessing" procedure to generate initial guess for 
#' iterative k-folds cross validation. This approach leverages some of the 
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
#' must be in the range \[0,1\], and the sum of all estimates within an 
#' occupational group must be <=1). Note that the modeling weights associated 
#' with guessed values are altered based on which of the three cases they fall 
#' into, with those falling in cases 1 and 3 receiving higher weights, and those
#' falling in case 2 receiving lower weights. These weights are used in the 
#' iterative modeling step.
#'
#' @param ors.data.for.sg Data with "missing" estimates to populate with smart 
#' guesses
#' @param ors.data Original data augmented with relevant predictors, i.e. all 
#' records, including both known and missing estimates (output of 
#' getPredictors(), or computeSimulations())
#' @param wt.low Model weight to assign to low-confidence smart guesses; default
#' is 0
#' @param wt.mid Model weight to assign to mid-confidence smart guesses; default
#' is 0.5
#' @param wt.high Model weight to assign to high-confidence smart guesses; 
#' default is 0.5
#' @param verbose Should messages be printed; default is FALSE (mute messages) 
#' @return Input data frame, with missing values filled in by smart guesses
#' @export
smartGuessKFCV <- function(ors.data.for.sg,ors.data,
                           wt.low=0,wt.mid=0.5,wt.high=0.5,
                           verbose=FALSE) {
  
  # Pre-populate known/missing values into prediction column
  ors.data.for.sg$prediction <- ors.data.for.sg[,"value"]
  
  # Set missing value weights to minimum
  ors.data.for.sg[which(ors.data.for.sg$known.val!=1),"weight"] <- wt.low
  
  # Column to flag predictions that are forced scaled/bounded on [0,1]
  ors.data.for.sg$sg.flag <- 0
  
  # Run smart guessing procedure
  for (adg in levels(ors.data.for.sg$additive_group)) {
    best.pct <- -1
    best.occ <- ""
    best.dist <- NA
    
    # Identify best distribution(s)
    if (verbose) {
      cat("\n",paste("Identifying most complete distribution for additive group: ",adg,sep=""),"\n")
    }
    for (occ in levels(ors.data.for.sg$occupation_text)) {
      current.occ.group <- ors.data.for.sg[as.character(ors.data.for.sg$occupation_text)==occ
                                           & as.character(ors.data.for.sg$additive_group)==adg,]
      known.pct <- sum(!is.na(current.occ.group[,"value"])) / nrow(current.occ.group)
      
      if (nrow(current.occ.group)==0) {
        next
      }
      else if (known.pct > best.pct[1]) {
        best.pct <- known.pct
        best.occ <- occ
        best.dist <- current.occ.group[,"value"]
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
        current.occ.group <- ors.data.for.sg[as.character(ors.data.for.sg$occupation_text)==bocc
                                             & as.character(ors.data.for.sg$additive_group)==adg,]
        current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
        bocc.tab[current.occ.group$to.sort,bocc] <- current.occ.group[,"value"]
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
    for (occ in levels(ors.data.for.sg$occupation_text)) {
      current.occ.group <- ors.data.for.sg[as.character(ors.data.for.sg$occupation_text)==occ
                                           & as.character(ors.data.for.sg$additive_group)==adg,]
      orig.occ.group <- ors.data[as.character(ors.data$occupation_text)==occ
                                 & as.character(ors.data$additive_group)==adg,]
      orig.levels <- nrow(orig.occ.group)
      
      if (nrow(current.occ.group)==0) {
        next
      }
      
      current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
      current.missing <- rownames(current.occ.group)[which(is.na(current.occ.group[,"value"]))]
      
      # Check that the current OG has missing values
      if (length(current.missing)==0) {
        next
      }
      
      # Calculate remaining percentage - take into account distribution of original data
      rem.pct <- length(current.missing) * ((1 - sum(current.occ.group[,"value"],na.rm=TRUE)) / (orig.levels-sum(!(is.na(current.occ.group[,"value"])))))
      if(rem.pct < 0) {
        rem.pct <- 0
      }
      
      # Reassign bounds (and pct.remaining)
      ors.data.for.sg[current.missing,"upper_bound"] <- rem.pct
      ors.data.for.sg[current.missing,"lower_bound"] <- 0
      ors.data.for.sg[current.missing,"pct.remaining"] <- rem.pct
      
      # Check that there is remaining percentage to assign
      if (rem.pct==0) {
        ors.data.for.sg[current.missing,"prediction"] <- 0
        ors.data.for.sg[current.missing,"weight"] <- wt.high
        next
      }
      
      # If there are values in the best distribution...
      if (sum(!is.na(best.dist)) >= 1) {
        
        # Determine which (if any) levels overlap between best distribution and current OG
        best.known <- names(best.dist)[which(!is.na(best.dist))]
        current.known <- current.occ.group$to.sort[which(!is.na(current.occ.group[,"value"]))]
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
            ors.data.for.sg[c(nj2.levels,jm2.levels),"sg.flag"] <- 1
            ors.data.for.sg[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else if ((nc1 * scaling.factor) > (1 - sum(c2,nc2))) {
            nj2 <- 1 - sum(c2,nc2)
            ors.data.for.sg[c(nj2.levels,jm2.levels),"sg.flag"] <- 1
            ors.data.for.sg[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else {
            nj2 <- nc1 * scaling.factor
            ors.data.for.sg[c(nj2.levels,jm2.levels),"weight"] <- wt.high
          }
          jm2 <- 1 - sum(c2,nc2,nj2)
          
          # Populate smart guesses for relevant observations
          if (length(nj2.levels) > 0 & length(jm2.levels) > 0) {
            ors.data.for.sg[nj2.levels,"prediction"] <- nj2 / length(nj2.levels)
            ors.data.for.sg[jm2.levels,"prediction"] <- jm2 / length(jm2.levels)
          }
          else if ((length(nj2.levels) > 0 & length(jm2.levels)==0)) {
            ors.data.for.sg[nj2.levels,"prediction"] <- rem.pct / length(nj2.levels)
          }
          else if ((length(nj2.levels)==0 & length(jm2.levels) > 0)) {
            ors.data.for.sg[jm2.levels,"prediction"] <- rem.pct / length(jm2.levels)
          }
          else {
            next
          }
          
          # Scale guesses to [0,1]
          # if (sum(ors.data.for.sg[rownames(current.occ.group),"prediction"])!=1) {
          if (sum(ors.data.for.sg[rownames(current.occ.group),"prediction"])>1) {
            ors.data.for.sg[current.missing,"prediction"] <- rem.pct * (ors.data.for.sg[current.missing,"prediction"] / sum(ors.data.for.sg[current.missing,"prediction"]))
          }
        }
        
        # If the current OG is an empty distribution (all NAs, i.e. special case of no overlap)...
        else if (sum(is.na(current.occ.group[,"value"]))==nrow(current.occ.group)) {  # make sure the OG is empty
          if (verbose) {
            print(paste("Smart guessing for occupation: ",occ,sep=""))
          }
          
          # Get observations that can take smart guesses/bounds (all)
          smart.guess.labels <- current.occ.group[which(is.na(current.occ.group[,"value"])),"to.sort"]
          
          # Populate values from best distribution as smart guesses for current OG
          ors.data.for.sg[current.missing,"prediction"] <- best.dist[smart.guess.labels]
          ors.data.for.sg[current.missing,"weight"] <- wt.high
          
          # Generate naive guesses for the others, if necessary
          # I.e. best distribution has some NAs, or current OG has observations not in best distribution
          # Latter is the result of using only known values (and simulating missing)
          if (sum(is.na(best.dist)) >= 1 | 
              sum(is.na(ors.data.for.sg[current.missing,"prediction"])) >= 1 ) {
            naive.guess <- current.missing[which(is.na(ors.data.for.sg[current.missing,"prediction"]))]
            ors.data.for.sg[naive.guess,"prediction"] <- (rem.pct/length(current.missing)) / length(naive.guess)
          }
        }
        
        # Otherwise (non-empty OG, but no overlap)...
        else {
          # Generate naive prediction
          ors.data.for.sg[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group[,"value"]))
          ors.data.for.sg[current.missing,"weight"] <- wt.low
        }
      }
      
      # If the best distribution is empty...
      else {
        # Generate naive prediction
        ors.data.for.sg[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group[,"value"]))
        ors.data.for.sg[current.missing,"weight"] <- wt.low
      }
    }
  }
  
  return(ors.data.for.sg)
}


#' K-folds CV, with Smart Guessing
#'
#' Function to do k-folds cross validation in an iterative fashion. The first 
#' prediction (Iteration 0) is the result of the smart guessing procedure. All
#' subsequent iterations rely on XGBoost to produce preliminary predictions. 
#' These predictions are then adjusted to adhere to boundary conditions on the 
#' data (all estimates must be in the range \[0,1\], and the sum of all estimates
#' within an occupational group must be <=1).
#'
#' @param ors.data Original data augmented with relevant predictors, i.e. all 
#' records, including both known and missing estimates (output of 
#' getPredictors(), or computeSimulations())
#' @param n.iter Number of times to iterate/adjust the model
#' @param weight.step Increment by which to increase modeling weight of test 
#' fold data with each iteration
#' @param mdl.d Tree model maximum depth; default is 14
#' @param mdl.n Tree model rounds; default is 200
#' @param mdl.e Tree model eta; default is 0.6
#' @param fold.list A list of preset folds, if these have to be fixed across 
#' runs; default is to generate a new set of folds
#' @param sg.soc.code SOC code to use for smart guessing, either "upSOC2", 
#' "upSOC3", or "upSOC4"
#' @return A list of length two, containing a list of hold out indices for each
#' test fold, and the results of iterative modeling
#' @export
iterateModelKFCV <- function(ors.data,n.iter,weight.step,
                             mdl.d=14,mdl.n=200,mdl.e=0.6,
                             fold.list=NULL,sg.soc.code) {
  
  # Get known values only
  ors.data.known <- ors.data[ors.data$known.val==1,]
  
  # Initialize list
  model.iterations <- list()
  
  # Get column names of predictors/response
  select.cols <- c("occupation_text","requirement","frequency","intensity","req.cat","upSOC2","upSOC3","prediction")
  
  doParallel::registerDoParallel(parallel::makeCluster(10)) # use 10 threads (1 per fold) for parallel cross-validation (CV)
  for (iter in c(0:n.iter)) {
    print(paste("Iteration ",iter,"...",sep=""))
    if(iter==0) {
      
      if (is.null(fold.list)) {
        print("Generating folds...")
        cv.folds <- caret::createFolds(ors.data.known$upSOC3,k=10,list=TRUE) # Create folds, balanced on SOC3 code
      }
      else {
        print("Using provided folds...")
        cv.folds <- fold.list
      }
      
      # 'dopar' will run this on multiple threads (change to just 'do' for synchronous runs)
      cv.results <- foreach::foreach(fold=cv.folds,fold.num=icount(),
                                     .packages=c('Matrix','xgboost','ModelMetrics','imputeORS')) %dopar% {
        
        ors.data.known$is.test <- 0
        ors.data.known$orig.value <- ors.data.known$value
        
        # Split data in to test/training folds
        data.train <- droplevels(ors.data.known[-fold,]) # Get the opposite of the test observations to train on
        data.test <- droplevels(ors.data.known[fold,])
        data.test$value <- NA
        
        # Binary indicator of test fold
        data.test$is.test <- 1
        
        # Merge training and test sets back together
        data.train <- rbind(data.train,data.test)
        
        # Separate observations by SOC2/3/4 code, and smart guess based on this distribution
        data.train.split <- list()
        for (i in levels(data.train[,sg.soc.code])) {
          data.train.split[[i]] <- droplevels(data.train[as.character(data.train[,sg.soc.code])==i,])
        }
        
        data.train.smart <- lapply(data.train.split,imputeORS::smartGuessKFCV,ors.data,wt.low=0,wt.mid=0.5,wt.high=0.5)
        data.train.smart <- do.call(rbind,data.train.smart)
        rownames(data.train.smart) <- gsub("^[0-9][0-9]*\\.","",rownames(data.train.smart))
        data.train <- data.train.smart[rownames(data.train),]
        rm(data.train.smart)
        
        # Make sure only true guesses have weight != 1
        data.train[which(data.train$prediction==data.train$value),"weight"] <- 1
        
        # Reset values
        data.train$value <- data.train$orig.value
        data.train$orig.value <- NULL
        
        # Return predictions vs. actual values, data, and fold number
        # Initial iteration (0) is just smart guessing
        list(results=data.frame(actual=data.train$value,current.prediction=data.train$prediction,
                                row.names=rownames(data.train)),
             data=data.train,which.test.fold=fold.num)
      }
      
      names(cv.results) <- paste("fold",c(1:length(cv.folds)),sep="")
      model.iterations[[iter+1]] <- cv.results
      
    }
    else {
      cv.results <- foreach::foreach(res=model.iterations[[iter]],fold.num=icount(),
                                     .packages=c('Matrix','xgboost','ModelMetrics','imputeORS')) %dopar% {
        
        # Get data from previous iteration
        data.train <- res$data
        
        # Update weights of predictions for new fit
        if (iter>1) {
          data.train[data.train$weight<1,"weight"] <- data.train[data.train$weight<1,"weight"] + weight.step
          data.train[which(data.train$weight<1 & data.train$weight>=0.75),"weight"] <- 0.75   # cap at 0.75
          data.train[data.train$weight==weight.step,"weight"] <- 0    # naive guesses get 0
        }
        
        # Format data and assign bounds
        sparse.train <- xgboost::xgb.DMatrix(data=Matrix::sparse.model.matrix(prediction ~ .,data=data.train[,select.cols])[,-1],label=data.train$prediction)
        xgboost::setinfo(sparse.train,'label_lower_bound',data.train$lower_bound)
        xgboost::setinfo(sparse.train,'label_upper_bound',data.train$upper_bound)
        xgboost::setinfo(sparse.train,'weight',data.train$weight)
        
        # Fit the model and make predictions
        # Note: https://github.com/uber/causalml/issues/96
        print(paste("Evaluating model ",fold.num,"...",sep=""))
        mdl <- xgboost::xgboost(data=sparse.train,booster="gbtree",max_depth=mdl.d,eta=mdl.e,
                                nthread=20,nrounds=mdl.n,objective="reg:squarederror")
        data.train$prediction <- predict(mdl,newdata=sparse.train)
        
        # Replace only estimates in test fold with their predicted values
        to.adjust <- rownames(data.train)[which(data.train$prediction!=data.train$value & data.train$is.test==1,)]
        data.train$temp.col <- data.train$value
        data.train[to.adjust,"temp.col"] <- data.train[to.adjust,"prediction"]
        data.train$prediction <- data.train$temp.col
        data.train$temp.col <- NULL
        
        # Correct/scale predictions based on known information
        data.train[data.train$prediction < 0,"prediction"] <- 0
        data.train[data.train$prediction > 1,"prediction"] <- 1
        for (occ in levels(data.train$occupation_text)) {
          for (adg in unique(data.train$additive_group)) {
            current.full <- data.train[data.train$occupation_text==occ & data.train$additive_group==adg,]
            current.to.adjust <- to.adjust[to.adjust %in% rownames(current.full)]
            
            if(length(current.to.adjust) > 0 & 
               sum(current.full$prediction) > 1 &
               length(unique(!is.na(data.train[current.to.adjust,"upper_bound"])))==1) {
              
              data.train[current.to.adjust,"prediction"] <- 
                (data.train[current.to.adjust,"prediction"] / sum(data.train[current.to.adjust,"prediction"])) * mean(data.train[current.to.adjust,"upper_bound"],na.rm=TRUE)
            }
          }
        }
        
        # Return predictions vs. actual values, model, data, and fold number
        list(results=data.frame(actual=data.train$value,current.prediction=data.train$prediction,
                                previous.prediction=res$results$current.prediction,row.names=rownames(data.train)),
             model=mdl,data=data.train,which.test.fold=fold.num)
      }
      
      names(cv.results) <- paste("fold",c(1:length(cv.folds)),sep="")
      model.iterations[[iter+1]] <- cv.results
      
    }
  }
  
  names(model.iterations) <- paste("iteration",c(0:n.iter),sep="")
  return(list(folds=cv.folds,model.iterations=model.iterations))
}


#' Collate predictions for all 10 test folds, for each iteration
#' 
#' Takes the results of a k-folds CV run and collates all the test fold data. 
#' In the KFCV analysis, each of 10 folds is treated as the test fold, and a 
#' model is generated and iterated over each of these test folds. The
#' predictions resulting from each test fold are combined for each iteration,
#' such that there is a complete predicted dataset (all records) obtained from
#' each iteration. Recall that Iteration 0 is simply the output of smart
#' guessing (no modeling involved).
#'
#' @param model.results Modeling results (output of iterateModelKFCVwSG())
#' @return Collated predictions of each test fold, by iteration
#' @export
getTestFoldData <- function(model.results) {
  model.iterations <- model.results$model.iterations
  
  # Set up final data frame
  test.fold.data <- data.frame(fold=rep(NA,nrow(model.iterations[[1]][[1]]$data)),
                               req.cat=model.iterations[[1]][[1]]$data$req.cat,
                               actual=model.iterations[[1]][[1]]$data$value,
                               row.names=rownames(model.iterations[[1]][[1]]$data))
  
  # By iteration...
  for(i in c(1:length(model.iterations))) {
    mi <- model.iterations[[i]]
    current.iter <- as.data.frame(matrix(nrow=0,ncol=3))
    
    # By fold...
    for(j in c(1:length(mi))){
      ids <- rownames(mi[[j]]$data[mi[[j]]$data$is.test==1,])
      
      # NOTE: "value" contains actual value
      # NOTE: "prediction" contains initial smart guess (iter0), along with subsequent adjusted predictions (iter1+)
      p <- mi[[j]]$data[ids,"prediction"]
      f <- rep(j,length(ids))
      
      current.iter <- rbind(current.iter,data.frame(prediction=p,fold=f,row.names=ids))
    }
    
    colnames(current.iter) <- paste(c("Prediction","fold"),(i-1),sep="")
    test.fold.data <- cbind(test.fold.data,current.iter[rownames(test.fold.data),])
    test.fold.data$fold <- test.fold.data[,ncol(test.fold.data)]
    test.fold.data[,ncol(test.fold.data)] <- NULL
  }
  
  return(test.fold.data)
}


#' Compute convergence iteration
#'
#' For a set of iterative k-folds modeling results, it is necessary to determine
#' the convergence point of the model. This function computes this using the 
#' RMSE of iterated predictions vs. the actual known values (this can be 
#' accomplished because the k-folds procedure relies on known estimates only,
#' and simulates missing estimates from the known set).
#' 
#' We defined convergence as when the difference between the RMSE of consecutive
#' iterations was < 0.001. If this convergence criteria is not met, the function
#' simply returns the final iteration (and a message stating that convergence 
#' was not reached). Recall that Iteration 0 is simply the output of smart 
#' guessing (no modeling involved).
#'
#' @param test.fold.data Test fold predictions (output of getTestFoldData())
#' @param verbose Should messages be printed; default is TRUE (print messages)
#' @return Convergence iteration of k-folds cross validation modeling
#' @export
computeConvergence <- function(test.fold.data,verbose=TRUE) {
  
  rmse.by.iter <- vector()
  for (i in c(1:(ncol(test.fold.data)-3))) {
    rmse.by.iter <- c(rmse.by.iter,ModelMetrics::rmse(test.fold.data[,"actual"],test.fold.data[,i+3]))
  }
  names(rmse.by.iter) <- paste("Prediction",c(0:(length(rmse.by.iter)-1)),sep="")
  
  i <- 0  # start at 0 to maintain consistency with iteration numbering (0-n)
  d <- 1
  while (d >= 0.001 & i < (length(rmse.by.iter)-1)) {
    d <- rmse.by.iter[i+1] - rmse.by.iter[i+2]
    i <- i+1
  }
  
  if (verbose) {
    cat("RMSE by iteration:\n")
    print(rmse.by.iter)
    cat("\nDifferences in RMSE:\n")
    cat(rmse.by.iter[1:(length(rmse.by.iter)-1)]-rmse.by.iter[2:length(rmse.by.iter)])
    cat("\n\n")
  }
  plot(0:(length(rmse.by.iter)-1),rmse.by.iter,xlab="Iteration",ylab="RMSE (prediction vs. actual)")
  
  if (d >= 0.001) {
    cat("Convergence criteria unmet, returning final iteration\n")
    return(i)
  }
  else {
    return(i-1)
  }
}


#' Plot predicted vs. actual values
#'
#' For each iteration, plot the predicted values vs. actual values of 
#' observations in test folds. Recall that Iteration 0 is simply the output of 
#' smart guessing (no modeling involved). Note that n.iter=10 will plot 
#' Iterations 0-10, for a total of 11 iterations.
#'
#' @param test.fold.data Test fold predictions (output of getTestFoldData())
#' @param plot.alpha Alpha (transparency) to use for points; default is 0.1
#' @param n.iter How many iterations are to be plotted; default is to plot all
#' iterations, including Iteration 0 (smart guessing)
#' @param reqcat.palette A character vector specifying colors to use for each 
#' requirement category; defaults are pink (COG), black (ENV), yellow (PHY), and
#' cyan (EDU)
#' @param line.color A string specifying color to use for the line y=x; default 
#' is red
#' @param titles A character vector specifying the titles of each subplot; 
#' default is the column names of test.fold.data
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @return Returns a list of plot objects, with one object for each iteration; 
#' optionally produces a plot (.png file) of all iterations
#' @export
plotTestFolds <- function(test.fold.data,plot.alpha=0.1,n.iter=(ncol(test.fold.data)-4),
                          reqcat.palette=c("palevioletred1","black","yellow","cyan"),
                          line.color=c("red"),
                          titles=gsub("([a-z])([0-9])","\\1 \\2",colnames(test.fold.data)),
                          print.plot=TRUE) {

  test.fold.data <- test.fold.data[,c(1:(n.iter+4))]
  test.fold.data$req.cat <- as.character(test.fold.data$req.cat)
  test.fold.data[test.fold.data$req.cat=="COG","req.cat"] <- "Cognitive and mental demands"
  test.fold.data[test.fold.data$req.cat=="ENV","req.cat"] <- "Environmental conditions"
  test.fold.data[test.fold.data$req.cat=="PHY","req.cat"] <- "Physical demands"
  test.fold.data[test.fold.data$req.cat=="SVP","req.cat"] <- "Education, training, and experience"
  test.fold.data$req.cat <- factor(test.fold.data$req.cat)
  # c("Cognitive and mental","Environmental conditions","Physical demands","Education, training, and experience")
  
  plts <- list()
  
  for(i in c(4:ncol(test.fold.data))) {
    selected.cols <- test.fold.data[,c("req.cat","actual",colnames(test.fold.data)[i])]
    colnames(selected.cols) <- c("req.cat","actual","calculated")
    
    plts[[(i-3)]] <- ggplot2::ggplot(selected.cols) + 
      ggplot2::geom_point(aes(x=calculated,y=actual,color=req.cat),alpha=plot.alpha) +
      ggplot2::scale_colour_manual(values=c("Cognitive and mental demands"=reqcat.palette[1],
                                            "Environmental conditions"=reqcat.palette[2],
                                            "Physical demands"=reqcat.palette[3],
                                            "Education, training, and experience"=reqcat.palette[4]),name="") +
      ggplot2::labs(title=paste(titles[i],"vs. Actual")) + ggplot2::xlab(titles[i]) + ggplot2::ylab("Actual") +
      ggplot2::geom_abline(intercept=0,slope=1,color=line.color) + 
      ggplot2::theme_bw() + ggplot2::theme(legend.position="bottom")
  }
  
  if (print.plot==TRUE) {
    plots.grid <- ggpubr::ggarrange(plotlist=plts,nrow=ceiling(length(plts)/2),ncol=2,labels="auto",
                                    # legend.grob=get_legend(plts[[1]]),
                                    # common.legend=TRUE,
                                    legend="none")
    ggplot2::ggsave("kfcv-residPlot.png",plots.grid,width=8.5,height=(ceiling(length(plts)/2)*4 + 0.36))
    # return(plots.grid)
  }
  
  return(plts)
}


#' Compute model blending ratio
#'
#' Our procedure involves running the model twice over. In the first run, the
#' initial smart guesses were computed based on separating the data based on
#' their SOC2 codes. In the second run, the same was done using the SOC3 codes
#' instead. This function computes a blending ratio, i.e. the contribution of
#' each model to our final predictions. This is done by determining the 
#' convergence iteration of each model, then computing a weighted average of the
#' relevant predictions in steps of 0.01, e.g. 0.01(run1) + 0.99*(run2). The 
#' RMSE of these prediction weighted averages were then calculated. The weighted
#' average that resulted in the lowest RMSE was selected to determine the
#' blending ratio. For example, if the prediction weighted average of 0.45(run1)
#' + 0.55(run2) yielded the lowest RMSE, then the blending ratio would be 45:55.
#' 
#' Note that we specify models based on SOC2 and SOC3 smart guesses to provide a
#' concrete example. In reality, the two models are not limited to this specific
#' combination.
#'
#' @param model.results.soc2 Modeling results of SOC2 smart guessed data (output
#' of iterateModelKFCVwSG())
#' @param model.results.soc3 Modeling results of SOC3 smart guessed data (output
#' of iterateModelKFCVwSG())
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @return A list of length five, containing the convergence iterations of
#' model.results.soc2 and model.results.soc3, the contributions of 
#' model.results.soc2 and model.results.soc3 (as proportions), and the blended
#' predictions of each iteration, which can be fed to plotTestFolds(); 
#' optionally produces a plot (.png file) of RMSE of blended predictions at 
#' convergence
#' @export
computeBlendingRatio <- function(model.results.soc2,model.results.soc3,print.plot=TRUE) {
  # Get aggregated test fold data
  test.folds.soc2 <- getTestFoldData(model.results.soc2)
  test.folds.soc3 <- getTestFoldData(model.results.soc3)
  
  # Compute convergence iterations
  cat("Computing convergence for model 1...\n\n")
  conv.iter.soc2 <- computeConvergence(test.folds.soc2)
  cat("Computing convergence for model 2...\n\n")
  conv.iter.soc3 <- computeConvergence(test.folds.soc3)
  cat(paste("Convergence iteration for model 1: ",conv.iter.soc2,"\n",sep=""))
  cat(paste("Convergence iteration for model 2: ",conv.iter.soc3,"\n",sep=""))
  conv.iter <- c(conv.iter.soc2,conv.iter.soc3)
  conv.iter.soc2 <- paste("Prediction",conv.iter.soc2,sep="")
  conv.iter.soc3 <- paste("Prediction",conv.iter.soc3,sep="")
  
  # Blending ratio
  blending.ratios <- vector()
  for (i in seq(0,1,0.01)) {
    lin.com <- (i*test.folds.soc2[rownames(test.folds.soc2),conv.iter.soc2]) + 
      ((1-i)*test.folds.soc3[rownames(test.folds.soc2),conv.iter.soc3])
    blending.ratios <- c(blending.ratios,ModelMetrics::rmse(test.folds.soc2$actual,lin.com))
  }
  rm(lin.com)
  blending.ratios <- data.frame(RMSE=blending.ratios,
                                SOC2.contribution=seq(0,1,0.01),
                                SOC3.contribution=rev(seq(0,1,0.01)),
                                color=0)
  blending.ratios[which.min(blending.ratios$RMSE),"color"] <- 1
  blending.ratios$color <- as.factor(blending.ratios$color)
  
  # Plot
  if (print.plot) {
    p.blend <- ggplot2::ggplot(blending.ratios,aes(x=SOC2.contribution,y=RMSE)) + 
      ggplot2::geom_point(aes(color=color),alpha=0.5) + 
      ggplot2::scale_color_manual(values=c("grey50","red")) + 
      ggplot2::labs(title="RMSE of blended predictions (at convergence)") + 
      ggplot2::theme_bw() + ggplot2::theme(legend.position="none") + 
      ggplot2::scale_x_continuous("Model 2 contribution",
                                  sec.axis=sec_axis(~ . *-1 + 1,name="Model 3 contribution"))
    ggplot2::ggsave(paste("kfcv-blended-rmse.png",sep=""),p.blend,width=6,height=3.75)
  }
  
  # Get proportions
  soc2.prop <- blending.ratios[blending.ratios$color==1,"SOC2.contribution"]
  soc3.prop <- blending.ratios[blending.ratios$color==1,"SOC3.contribution"]
  
  # Blended data
  test.folds.blend <- test.folds.soc2[,1:3]
  test.folds.blend <- cbind(test.folds.blend,
                            (soc2.prop*test.folds.soc2[,4:ncol(test.folds.soc2)] + 
                               soc3.prop*test.folds.soc3[rownames(test.folds.soc2),4:ncol(test.folds.soc3)]))
  
  # Return data
  return(list(soc2.conv.iter=conv.iter[1],
              soc3.conv.iter=conv.iter[2],
              soc2.proportion=soc2.prop,
              soc3.proportion=soc3.prop,
              test.folds.blended=test.folds.blend))
}



# Determine neighbor counts
# neighbors <- vector()
# for (i in levels(iter4$upSOC3)) {
#   for (j in levels(iter4$characteristic)) {
#     x.temp <- nrow(iter4[iter4$upSOC3==i & iter4$characteristic==j & 
#                            (iter4$data_type_text=="CONSTANTLY" |
#                               iter4$data_type_text=="FREQUENTLY" |
#                               iter4$data_type_text=="OCCASIONALLY" |
#                               iter4$data_type_text=="SELDOM" |
#                               iter4$data_type_text=="NOT PRESENT"),])
#     if (x.temp!=0) {
#       neighbors <- c(neighbors,
#                      nrow(iter4[iter4$upSOC3==i & iter4$characteristic==j & iter4$data_type_text=="CONSTANTLY",]))
#     }
#   }
# }
# rm(x.temp)
