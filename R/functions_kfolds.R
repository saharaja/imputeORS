# Functions for k-folds cross validation

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
#' setDefaultModelingWeights(), or computeSimulations())
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
  
  # Initialize necessary columns
  ors.data$is.test <- 0
  ors.data$orig.value <- ors.data$value
  
  # Separate known and missing values
  ors.data.known <- droplevels(ors.data[which(ors.data$known.val==1),])
  ors.data.missing <- droplevels(ors.data[which(ors.data$known.val!=1),])
  
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
                                      
        # Split known data in to test/training folds
        data.train <- droplevels(ors.data.known[-fold,]) # Get the opposite of the test observations to train on
        data.test <- droplevels(ors.data.known[fold,])
        data.test$value <- NA
        
        # Binary indicator of test fold
        data.test$is.test <- 1
        
        # Merge training and test sets back together
        data.train <- rbind(data.train,data.test)

        # Merge known and missing data back together
        data.train <- rbind(data.train,ors.data.missing)
        
        # Separate observations by SOC2/3/4 code, and smart guess based on this distribution
        data.train.split <- list()
        for (i in levels(data.train[,sg.soc.code])) {
          data.train.split[[i]] <- droplevels(data.train[as.character(data.train[,sg.soc.code])==i,])
        }
        
        data.train.smart <- lapply(data.train.split,imputeORS::smartGuess,ors.data,wt.low=0,wt.mid=0.5,wt.high=0.5)
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
                                is.test=data.train$is.test,known.val=data.train$known.val,
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
        mdl <- xgboost::xgboost(data=sparse.train,booster="gbtree",verbose=0,max_depth=mdl.d,eta=mdl.e,
                                nthread=20,nrounds=mdl.n,objective="reg:squarederror")
        data.train$prediction <- predict(mdl,newdata=sparse.train)
        
        # Correct/scale predictions based on known information
        print("Adjusting predictions...")
        # Bounds
        data.train[data.train$prediction < 0,"prediction"] <- 0
        data.train[data.train$prediction > 1,"prediction"] <- 1
        # Known values
        to.reset <- rownames(data.train)[which(data.train$known.val==1 & data.train$is.test==0)]
        data.train[to.reset,"prediction"] <- data.train[to.reset,"value"]
        # Percentage rules
        to.adjust <- rownames(data.train)[which((data.train$known.val!=1) | 
                                                  (data.train$known.val==1 & data.train$is.test==1))]
        for (occ in levels(data.train$occupation_text)) {
          for (adg in unique(data.train$additive_group)) {
            current.full <- data.train[data.train$occupation_text==occ & data.train$additive_group==adg,]
            current.to.adjust <- to.adjust[to.adjust %in% rownames(current.full)]
            
            # Additive groups not summing to 1
            if(nrow(current.full) >= 1 &
               length(current.to.adjust) > 0 & 
               sum(current.full$prediction) != 1 &
               sum(current.full[current.to.adjust,"prediction"]) != 0) {
              
              x <- data.train[current.to.adjust,"prediction"] / sum(data.train[current.to.adjust,"prediction"])
              y <- 1 - sum(current.full[current.full$is.test==0,"value"],na.rm=TRUE)
              if (y < 0) {  # an unlikely edge case
                data.train[current.to.adjust,"prediction"] <- 0
              }
              else {
                data.train[current.to.adjust,"prediction"] <- x * y
              }
            }
            
            # Check for additive groups that sum to <1 and have all missing estimates with (current) prediction of 0
            # This will happen due to bounding adjustments (and above <if> statement will change nothing)
            # For these observations, assign original smart guess
            if (nrow(current.full) >= 1 &
                length(current.to.adjust) > 0 &
                sum(current.full$prediction) < 1 &
                sum(current.full[current.to.adjust,"prediction"]) == 0) {
              
              data.train[current.to.adjust,"prediction"] <- model.iterations[[1]][[fold.num]]$data[current.to.adjust,"prediction"]
            }
          }
        }

        # Return predictions vs. actual values, model, data, and fold number
        list(results=data.frame(actual=data.train$value,current.prediction=data.train$prediction,
                                previous.prediction=res$results$current.prediction,
                                is.test=data.train$is.test,known.val=data.train$known.val,
                                row.names=rownames(data.train)),
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
  template.data <- model.iterations[[1]][[1]]$data[which(model.iterations[[1]][[1]]$data$known.val==1),]
  test.fold.data <- data.frame(fold=rep(NA,nrow(template.data)),
                               req.cat=template.data$req.cat,
                               actual=template.data$value,
                               row.names=rownames(template.data))
  rm(template.data)
  
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
#' @param show.plot Should plot of RMSE by iteration be displayed; default is 
#' TRUE (display plot)
#' @return Convergence iteration of k-folds cross validation modeling
#' @export
computeConvergence <- function(test.fold.data,verbose=TRUE,show.plot=TRUE) {
  
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
  if (show.plot) {
    plot(0:(length(rmse.by.iter)-1),rmse.by.iter,xlab="Iteration",ylab="RMSE (prediction vs. actual)")
  }
  
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
#' blending ratio. For example, if the prediction weighted average of 0.45(run1) + 
#' 0.55(run2) yielded the lowest RMSE, then the blending ratio would be 45:55.
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
#' @return A list of length six, containing the convergence iterations of
#' model.results.soc2 and model.results.soc3, the contributions of 
#' model.results.soc2 and model.results.soc3 (as proportions), the blended
#' predictions of each iteration (which can be fed to plotTestFolds()), and a 
#' plot object of the RMSE of blended predictions at convergence; optionally 
#' produces a plot (.png file) of RMSE of blended predictions at convergence
#' @export
computeBlendingRatio <- function(model.results.soc2,model.results.soc3,print.plot=TRUE) {
  # Get aggregated test fold data
  test.folds.soc2 <- getTestFoldData(model.results.soc2)
  test.folds.soc3 <- getTestFoldData(model.results.soc3)
  
  # Compute convergence iterations
  cat("Computing convergence for model 1...\n")
  conv.iter.soc2 <- computeConvergence(test.folds.soc2,verbose=FALSE,show.plot=FALSE)
  cat("\nComputing convergence for model 2...\n")
  conv.iter.soc3 <- computeConvergence(test.folds.soc3,verbose=FALSE,show.plot=FALSE)
  cat(paste("\nConvergence iteration for model 1: ",conv.iter.soc2,"\n",sep=""))
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
  p.blend <- ggplot2::ggplot(blending.ratios,aes(x=SOC2.contribution,y=RMSE)) + 
    ggplot2::geom_point(aes(color=color),alpha=0.5) + 
    ggplot2::scale_color_manual(values=c("grey50","red")) + 
    ggplot2::labs(title="RMSE of blended predictions (at convergence)") + 
    ggplot2::theme_bw() + ggplot2::theme(legend.position="none") + 
    ggplot2::scale_x_continuous("Model 2 contribution",
                                sec.axis=sec_axis(~ . *-1 + 1,name="Model 3 contribution"))
  if (print.plot) {
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
              test.folds.blended=test.folds.blend,
              plot.blend=p.blend))
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
