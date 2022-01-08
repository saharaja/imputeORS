# Functions for train-test portion of model tuning

#' Train-test approach for model tuning
#'
#' Function to do train-test modeling in an iterative fashion, for the purposes 
#' of model tuning. The first prediction (Iteration 0) is just the naive guess. 
#' All subsequent iterations rely on XGBoost to produce preliminary predictions. 
#' These predictions are then adjusted to adhere to boundary conditions on the 
#' data (all estimates must be in the range \[0,1\], and the sum of all estimates
#' within an occupational group must be <=1).
#'
#' @param ors.data Original data augmented with relevant predictors, i.e. all 
#' records, including both known and missing estimates (output of 
#' [setDefaultModelingWeights()], or [computeSimulations()])
#' @param n.iter Number of times to iterate/adjust the model
#' @param weight.step Increment by which to increase modeling weight of test 
#' set data with each iteration
#' @param mdl.d Tree model maximum depth; default is 14
#' @param mdl.n Tree model rounds; default is 200
#' @param mdl.e Tree model eta; default is 0.6
#' @param fold.list A list of preset folds, if these have to be fixed across 
#' runs; default is to generate a new set of folds
#' @return A list of length two, containing a list of hold out indices for the
#' test set, and the results of iterative modeling
#' @seealso [setDefaultModelingWeights()]
#' @seealso [computeSimulations()]
#' @export
iterateModelTT <- function(ors.data,n.iter,weight.step,
                           mdl.d=14,mdl.n=200,mdl.e=0.6,
                           fold.list=NULL) {
  
  # Get known values only
  ors.data.known <- ors.data[ors.data$known.val==1,]
  
  # Initialize list
  model.iterations <- list()
  
  # Get column names of predictors/response
  select.cols <- c("occupation_text","requirement","frequency","intensity","req.cat","upSOC2","upSOC3","prediction")
  
  for (iter in c(0:n.iter)) {
    print(paste("Iteration ",iter,"...",sep=""))
    if(iter==0) {
      
      if (is.null(fold.list)) {
        print("Generating test/train/validate sets...")
        tt.folds <- caret::createFolds(ors.data.known$upSOC3,k=10,list=TRUE) # Create folds, balanced on SOC3 code
        tt.folds$GroupA <- unlist(tt.folds[1:8],recursive=FALSE,use.names=FALSE)
        tt.folds$GroupB <- tt.folds[[9]]
        tt.folds$GroupC <- tt.folds[[10]]
      }
      else {
        print("Using provided sets...")
        tt.folds <- fold.list
      }
      
      ors.data.known$is.test <- 0
      ors.data.known$prediction <- ors.data.known$value
      
      # Split data in to test/training folds
      data.train <- droplevels(ors.data.known[tt.folds$GroupA,]) # Get the opposite of the test observations to train on
      data.test <- droplevels(ors.data.known[tt.folds$GroupB,])
      data.test$prediction <- NA
      
      # Binary indicator of test fold
      data.test$is.test <- 1
      
      # Naive Guess:
      # For observations in test set, use averaged percent remaining for value
      # This is simply the naive guess (initial prediction)
      # Have to recalculate percent remaining based on data in test vs. training sets
      print("Generating guesses in test set...")
      for (occ in levels(data.test$occupation_text)) {
        for (adg in unique(data.test$additive_group)) {
          current.full <- ors.data[as.character(ors.data$occupation_text)==occ & 
                                     ors.data$additive_group==adg,]
          current.test <- data.test[as.character(data.test$occupation_text)==occ & data.test$additive_group==adg,]
          
          if(nrow(current.test) > 0) {
            new.pct.remaining <- 1 - (sum(current.full$value,na.rm=TRUE) - sum(current.test$value,na.rm=TRUE))
            
            data.test[rownames(current.test),"prediction"] <- new.pct.remaining/(sum(is.na(current.full$value)) + nrow(current.test))
            data.test[rownames(current.test),"lower_bound"] <- 0
            data.test[rownames(current.test),"upper_bound"] <- new.pct.remaining
            data.test[rownames(current.test),"weight"] <-  0.25
          }
        }
      }
      
      # Merge training and test sets back together
      data.train <- rbind(data.train,data.test)
      
      # Make sure only true guesses have weight != 1
      data.train[which(data.train$prediction==data.train$value),"weight"] <- 1
      
      # Return predictions vs. actual values, data, and fold number
      # Initial iteration (0) is just naive guessing
      tt.results <- list(results=data.frame(actual=data.train$value,
                                            current.prediction=data.train$prediction,
                                            row.names=rownames(data.train)),
                         data=data.train)
      
      model.iterations[[iter+1]] <- tt.results
      
    }
    else {  
      # Get data from previous iteration
      data.train <- model.iterations[[iter]]$data
      
      # Update weights of predictions for new fit
      if (iter>1) {
        data.train[data.train$weight<1,"weight"] <- data.train[data.train$weight<1,"weight"] + weight.step
        data.train[which(data.train$weight<1 & data.train$weight>=0.75),"weight"] <- 0.75   # cap at 0.75
        # data.train[data.train$weight==weight.step,"weight"] <- 0    # naive guesses get 0 - for when smartGuess is used
      }
      
      # Format data and assign bounds
      sparse.train <- xgboost::xgb.DMatrix(data=Matrix::sparse.model.matrix(prediction ~ .,data=data.train[,select.cols])[,-1],label=data.train$prediction)
      xgboost::setinfo(sparse.train,'label_lower_bound',data.train$lower_bound)
      xgboost::setinfo(sparse.train,'label_upper_bound',data.train$upper_bound)
      xgboost::setinfo(sparse.train,'weight',data.train$weight)
      
      # Fit the model and make predictions
      # Note: https://github.com/uber/causalml/issues/96
      print(paste("Evaluating model...",sep=""))
      mdl <- xgboost::xgboost(data=sparse.train,booster="gbtree",max_depth=mdl.d,eta=mdl.e,
                              nthread=20,nrounds=mdl.n,objective="reg:squarederror",verbose=0)
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
      tt.results <- list(results=data.frame(actual=data.train$value,
                                            current.prediction=data.train$prediction,
                                            previous.prediction=model.iterations[[iter]]$results$current.prediction,
                                            row.names=rownames(data.train)),
                         model=mdl,data=data.train)
      
      model.iterations[[iter+1]] <- tt.results
      
    }
  }
  
  names(model.iterations) <- paste("iteration",c(0:n.iter),sep="")
  return(list(folds=tt.folds,model.iterations=model.iterations))
}


#' Collate predictions for each iteration for the train-test approach
#' 
#' Takes the results of a train-test run and collates all the test set data. 
#' The model is generated and iterated over the test set, producing predictions
#' for about 10% of the data. These predictions are extracted and organized by
#' iteration. Recall that Iteration 0 is simply the naive guess (no modeling 
#' involved).
#'
#' @param model.results Modeling results (output of [iterateModelTT()])
#' @return Predictions from the test set, by iteration
#' @seealso [iterateModelTT()]
#' @export
getTestSetData <- function(model.results) {
  model.iterations <- model.results$model.iterations
  
  # Set up final data frame
  test.obs <- rownames(model.iterations[[1]]$data)[which(model.iterations[[1]]$data$is.test==1)]
  test.set.data <- data.frame(req.cat=model.iterations[[1]]$data[test.obs,"req.cat"],
                              actual=model.iterations[[1]]$data[test.obs,"value"],
                              row.names=test.obs)
  
  # Get predictions by iteration
  for(i in c(1:length(model.iterations))) {
    mi <- model.iterations[[i]]
    p <- mi$data[test.obs,"prediction"]
    test.set.data <- cbind(test.set.data,p)
  }
  
  colnames(test.set.data) <- c("req.cat","actual",paste("Prediction",c(0:(i-1)),sep=""))
  return(test.set.data)
}


