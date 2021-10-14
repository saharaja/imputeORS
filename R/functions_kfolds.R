# Functions for k-folds cross validation


#' Smart Guessing for k-folds CV
#'
#' DETAILED DESCRIPTION lord have mercy on my soul
#'
#' @param ors.data Data with "missing" estimates to populate with smart guesses
#' @param original.ors.full Original data (all records, including both known and 
#' missing estimates)
#' @param sim.no Is the smart guessing being run on a simulation, or the 
#' original value; default is original value
#' @param wt.low Model weight to assign to low-confidence smart guesses; default
#' is 0
#' @param wt.mid Model weight to assign to mid-confidence smart guesses; default
#' is 0.5
#' @param wt.high Model weight to assign to high-confidence smart guesses; 
#' default is 0.5
#' @return Original data frame, with missing values filled in by smart guesses
#' @export
smartGuessKFCV <- function(ors.data,original.ors.full,sim.no="value",
                           wt.low=0,wt.mid=0.5,wt.high=0.5) {
  if(sim.no=="value") {
    val.col <- sim.no
  }
  else if(is.numeric(sim.no)) {
    val.col <- paste("valsim",sim.no,sep="")
  }
  
  # Pre-populate known/missing values into prediction column
  ors.data$prediction <- ors.data[,val.col]
  
  # Set missing value weights to minimum
  ors.data[which(ors.data$known.val!=1),"weight"] <- wt.low
  
  # Column to flag predictions that are forced scaled/bounded on [0,1]
  ors.data$sg.flag <- 0
  
  # Run smart guessing procedure
  for (adg in levels(ors.data$additive_group)) {
    best.pct <- -1
    best.occ <- ""
    best.dist <- NA
    
    # Identify best distribution(s)
    cat("\n",paste("Identifying most complete distribution for additive group: ",adg,sep=""),"\n")
    for (occ in levels(ors.data$occupation_text)) {
      current.occ.group <- ors.data[as.character(ors.data$occupation_text)==occ
                                    & as.character(ors.data$additive_group)==adg,]
      known.pct <- sum(!is.na(current.occ.group[,val.col])) / nrow(current.occ.group)
      
      if (nrow(current.occ.group)==0) {
        next
      }
      else if (known.pct > best.pct[1]) {
        best.pct <- known.pct
        best.occ <- occ
        best.dist <- current.occ.group[,val.col]
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
      print("*** Averaging over multiple best distributions...")
      bocc.tab <- data.frame(matrix(nrow=nrow(current.occ.group),ncol=length(best.occ)),
                             row.names=paste(current.occ.group$frequency,current.occ.group$intensity,sep="-"))
      colnames(bocc.tab) <- best.occ
      
      for (bocc in best.occ) {
        current.occ.group <- ors.data[as.character(ors.data$occupation_text)==bocc
                                      & as.character(ors.data$additive_group)==adg,]
        current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
        bocc.tab[current.occ.group$to.sort,bocc] <- current.occ.group[,val.col]
      }
      best.dist <- rowMeans(bocc.tab,na.rm=TRUE)
      names(best.dist) <- rownames(bocc.tab)
    }
    best.dist[is.nan(best.dist)] <- NA
    if(sum(best.dist,na.rm=TRUE) > 1) {
      best.dist <- best.dist / sum(best.dist,na.rm=TRUE)
    }
    # print(best.dist)
    
    # Assign upper bounds and guesses by occupational group (OG)
    for (occ in levels(ors.data$occupation_text)) {
      current.occ.group <- ors.data[as.character(ors.data$occupation_text)==occ
                                    & as.character(ors.data$additive_group)==adg,]
      orig.occ.group <- original.ors.full[as.character(original.ors.full$occupation_text)==occ
                                          & as.character(original.ors.full$additive_group)==adg,]
      orig.levels <- nrow(orig.occ.group)
      
      if (nrow(current.occ.group)==0) {
        next
      }
      
      current.occ.group$to.sort <- paste(current.occ.group$frequency,current.occ.group$intensity,sep="-")
      current.missing <- rownames(current.occ.group)[which(is.na(current.occ.group[,val.col]))]
      
      # Check that the current OG has missing values
      if (length(current.missing)==0) {
        next
      }
      
      # Calculate remaining percentage - take into account distribution of original data
      rem.pct <- length(current.missing) * ((1 - sum(current.occ.group[,val.col],na.rm=TRUE)) / (orig.levels-sum(!(is.na(current.occ.group[,val.col])))))
      if(rem.pct < 0) {
        rem.pct <- 0
      }
      
      # Reassign bounds (and pct.remaining)
      ors.data[current.missing,"upper_bound"] <- rem.pct
      ors.data[current.missing,"lower_bound"] <- 0
      ors.data[current.missing,"pct.remaining"] <- rem.pct
      
      # Check that there is remaining percentage to assign
      if (rem.pct==0) {
        ors.data[current.missing,"prediction"] <- 0
        ors.data[current.missing,"weight"] <- wt.high
        next
      }
      
      # If there are values in the best distribution...
      if (sum(!is.na(best.dist)) >= 1) {
        
        # Determine which (if any) levels overlap between best distribution and current OG
        best.known <- names(best.dist)[which(!is.na(best.dist))]
        current.known <- current.occ.group$to.sort[which(!is.na(current.occ.group[,val.col]))]
        common.known <- intersect(best.known,current.known)
        
        # If there is overlap between best distribution and current OG..
        if (length(common.known) > 0) {  # make sure there is overlap
          # print(paste("Smart guessing for occupation: ",occ,sep=""))
          
          # Common values
          c1 <- sum(best.dist[common.known],na.rm=TRUE)
          nc1 <- sum(best.dist,na.rm=TRUE) - c1
          c2 <- sum(current.occ.group[which(current.occ.group$to.sort %in% common.known),val.col],na.rm=TRUE)
          nc2 <- sum(current.occ.group[,val.col],na.rm=TRUE) - c2
          scaling.factor <- (1-c2) / (1-c1)
          
          # Levels of relevant missing values (those that violate boundary conditions)
          nj2.levels <- rownames(current.occ.group)[which(current.occ.group$to.sort %in% best.known & !(current.occ.group$to.sort %in% current.known))]
          jm2.levels <- rownames(current.occ.group)[which(!(current.occ.group$to.sort %in% best.known | current.occ.group$to.sort %in% current.known))]
          
          # Missing values
          if (is.na(scaling.factor) | scaling.factor==Inf) {
            nj2 <- 1 - sum(c2,nc2)
            ors.data[c(nj2.levels,jm2.levels),"sg.flag"] <- 1
            ors.data[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else if ((nc1 * scaling.factor) > (1 - sum(c2,nc2))) {
            nj2 <- 1 - sum(c2,nc2)
            ors.data[c(nj2.levels,jm2.levels),"sg.flag"] <- 1
            ors.data[c(nj2.levels,jm2.levels),"weight"] <- wt.mid
          }
          else {
            nj2 <- nc1 * scaling.factor
            ors.data[c(nj2.levels,jm2.levels),"weight"] <- wt.high
          }
          jm2 <- 1 - sum(c2,nc2,nj2)
          
          # Populate smart guesses for relevant observations
          if (length(nj2.levels) > 0 & length(jm2.levels) > 0) {
            ors.data[nj2.levels,"prediction"] <- nj2 / length(nj2.levels)
            ors.data[jm2.levels,"prediction"] <- jm2 / length(jm2.levels)
          }
          else if ((length(nj2.levels) > 0 & length(jm2.levels)==0)) {
            ors.data[nj2.levels,"prediction"] <- rem.pct / length(nj2.levels)
          }
          else if ((length(nj2.levels)==0 & length(jm2.levels) > 0)) {
            ors.data[jm2.levels,"prediction"] <- rem.pct / length(jm2.levels)
          }
          else {
            next
          }
          
          # Scale guesses to [0,1]
          # if (sum(ors.data[rownames(current.occ.group),"prediction"])!=1) {
          if (sum(ors.data[rownames(current.occ.group),"prediction"])>1) {
            ors.data[current.missing,"prediction"] <- rem.pct * (ors.data[current.missing,"prediction"] / sum(ors.data[current.missing,"prediction"]))
          }
        }
        
        # If the current OG is an empty distribution (all NAs, i.e. special case of no overlap)...
        else if (sum(is.na(current.occ.group[,val.col]))==nrow(current.occ.group)) {  # make sure the OG is empty
          print(paste("Smart guessing for occupation: ",occ,sep=""))
          
          # Get observations that can take smart guesses/bounds (all)
          smart.guess.labels <- current.occ.group[which(is.na(current.occ.group[,val.col])),"to.sort"]
          
          # Populate values from best distribution as smart guesses for current OG
          ors.data[current.missing,"prediction"] <- best.dist[smart.guess.labels]
          ors.data[current.missing,"weight"] <- wt.high
          
          # Generate naive guesses for the others, if necessary
          # I.e. best distribution has some NAs, or current OG has observations not in best distribution
          # Latter is the result of using only known values (and simulating missing)
          if (sum(is.na(best.dist)) >= 1 | 
              sum(is.na(ors.data[current.missing,"prediction"])) >= 1 ) {
            naive.guess <- current.missing[which(is.na(ors.data[current.missing,"prediction"]))]
            ors.data[naive.guess,"prediction"] <- (rem.pct/length(current.missing)) / length(naive.guess)
          }
        }
        
        # Otherwise (non-empty OG, but no overlap)...
        else {
          # Generate naive prediction
          ors.data[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group[,val.col]))
          ors.data[current.missing,"weight"] <- wt.low
        }
      }
      
      # If the best distribution is empty...
      else {
        # Generate naive prediction
        ors.data[current.missing,"prediction"] <- rem.pct / sum(is.na(current.occ.group[,val.col]))
        ors.data[current.missing,"weight"] <- wt.low
      }
    }
  }
  
  return(ors.data)
}


#' K-folds CV with Smart Guessing
#'
#' DETAILED DESCRIPTION lord have mercy on my soul
#'
#' @param original.ors.full Original data (all records, including both known and 
#' missing estimates)
#' @param n.iter Number of times to adjust/iterate the model
#' @param initial.guess.weight
#' @param prediction.contribution
#' @param guess.contribution
#' @param weight.step
#' @param mdl.d
#' @param mdl.n
#' @param mdl.e
#' @param fold.list
#' @param sg.soc.code SOC code to use for smart guessing, either "upSOC2", 
#' "upSOC3", or "upSOC4"
#' @return 
#' @export
iterateModelKFCVwSG <- function(original.ors.full,n.iter,initial.guess.weight,
                                prediction.contribution,guess.contribution,weight.step,
                                mdl.d,mdl.n,mdl.e,fold.list=NULL,sg.soc.code) {
  
  # Get known values only
  original.ors.known <- original.ors.full[original.ors.full$known.val==1,]
  
  # Initialize list
  model.iterations <- list()
  
  # Get column names of predictors/response
  select.cols <- c("occupation_text","characteristic","frequency","intensity","req.cat","upSOC2","upSOC3","value")
  
  for (iter in c(0:n.iter)) {
    print(paste("Iteration ",iter,"...",sep=""))
    if(iter==0) {
      
      if (is.null(fold.list)) {
        cv.folds <- caret::createFolds(original.ors.known$upSOC3,k=10,list=TRUE) # Create folds, balanced on SOC3 code
        print("Generating folds...")
      }
      else {
        cv.folds <- fold.list
        print("Using provided folds...")
      }
      
      
      # 'dopar' will run this on multiple threads (change to just 'do' for synchronous runs)
      cv.results <- foreach(fold=cv.folds,fold.num=icount(),.packages=c('Matrix','xgboost','ModelMetrics','imputeORS')) %dopar% {
        
        original.ors.known$orig.value <- original.ors.known$value
        original.ors.known$is.test <- 0
        
        # Split data in to test/training folds
        data.train <- droplevels(original.ors.known[-fold,]) # Get the opposite of the test observations to train on
        data.test <- droplevels(original.ors.known[fold,])
        data.test$value <- NA
        if (is.numeric(initial.guess.weight)) {
          data.test$weight <- initial.guess.weight
        }
        
        # Binary indicator of test fold
        data.test$is.test <- 1
        
        # Merge training and test sets back together
        data.train <- rbind(data.train,data.test)
        
        # Separate observations by SOC2/3/4 code, and smart guess based on this distribution
        data.train.split <- list()
        for (i in levels(data.train[,sg.soc.code])) {
          data.train.split[[i]] <- droplevels(data.train[as.character(data.train[,sg.soc.code])==i,])
        }
        data.train.smart <- lapply(data.train.split,smartGuessKFCV,original.ors.full,"value",wt.low=0,wt.mid=0.5,wt.high=0.5)
        data.train.smart <- do.call(rbind,data.train.smart)
        
        rownames(data.train.smart) <- gsub("^[0-9][0-9]*\\.","",rownames(data.train.smart))
        data.train.smart$orig.value <- data.train[rownames(data.train.smart),"orig.value"]
        data.train.smart$is.test <- data.train[rownames(data.train.smart),"is.test"]
        
        data.train <- data.train.smart[rownames(data.train),]
        data.train$value <- data.train$prediction
        data.train$prediction <- NULL
        
        # Make sure only true guesses have weight != 1
        data.train[which(data.train$value==data.train$orig.value),"weight"] <- 1
        
        # Format data and assign bounds
        sparse.train <- xgb.DMatrix(data=sparse.model.matrix(value ~ .,data=data.train[,select.cols])[,-1],label=data.train$value)
        setinfo(sparse.train,'label_lower_bound',data.train$lower_bound)
        setinfo(sparse.train,'label_upper_bound',data.train$upper_bound)
        setinfo(sparse.train,'weight',data.train$weight)
        
        # Fit the model and make predictions
        print(paste("Evaluating model ",fold.num,"...",sep=""))
        mdl <- xgboost(data=sparse.train,booster="gbtree",max_depth=mdl.d,eta=mdl.e,
                       nthread=20,nrounds=mdl.n,objective="reg:squarederror") # Note: https://github.com/uber/causalml/issues/96
        train.est.pred <- predict(mdl,newdata=sparse.train)
        train.est.val <- data.train$value
        
        # Return predictions vs. actual values, RMSE, model, fold number, hold out fold indices, and data
        list(results=data.frame(val=train.est.val,pred=train.est.pred,row.names=rownames(data.train)),
             results.rmse=rmse(train.est.val,train.est.pred),
             model=mdl,fold=fold.num,data=data.train)
      }
      
      model.iterations[[iter]] <- cv.results
      
    }
    else {
      cv.results <- foreach(res=model.iterations[[iter-1]],fold.num=icount(),.packages=c('Matrix','xgboost','ModelMetrics','imputeORS')) %dopar% {
        
        # Get data from previous fit
        data.train <- res$data
        
        # Update weights of guesses to reflect new guess (assigned below)
        if (is.numeric(initial.guess.weight)) {
          data.train[data.train$weight<1,"weight"] <- min(initial.guess.weight+((iter-1)*weight.step),0.75)
        }
        else {
          data.train[data.train$weight<1,"weight"] <- data.train[data.train$weight<1,"weight"] + weight.step
          data.train[which(data.train$weight<1 & data.train$weight>=0.75),"weight"] <- 0.75   # cap at 0.75
          data.train[data.train$weight==weight.step,"weight"] <- 0    # naive guesses get 0
        }
        
        # Get previous predictions
        data.train$prev.pred <- res$results$pred
        
        # For observations in test fold, use previous guess and previous prediction to generate new guess
        print("Recalculating values in test fold...")
        to.adjust <- rownames(data.train)[which(data.train$value!=data.train$orig.value & data.train$is.test==1,)]
        data.train[to.adjust,"value"] <- ((prediction.contribution*data.train[to.adjust,"prev.pred"]) + 
                                            (guess.contribution*data.train[to.adjust,"value"]))
        data.train$prev.pred <- NULL
        
        # Correct/scale new guesses based on known information
        data.train[data.train$value < 0,"value"] <- 0
        data.train[data.train$value > 1,"value"] <- 1
        for (occ in levels(data.train$occupation_text)) {
          for (adg in unique(data.train$additive_group)) {
            current.full <- data.train[data.train$occupation_text==occ & data.train$additive_group==adg,]
            current.to.adjust <- to.adjust[to.adjust %in% rownames(current.full)]
            
            if(length(current.to.adjust) > 0 & 
               sum(current.full$value) > 1 &
               length(unique(!is.na(data.train[current.to.adjust,"upper_bound"])))==1) {
              
              data.train[current.to.adjust,"value"] <- 
                (data.train[current.to.adjust,"value"] / sum(data.train[current.to.adjust,"value"])) * mean(data.train[current.to.adjust,"upper_bound"],na.rm=TRUE)
            }
          }
        }
        
        # Format data and assign bounds
        sparse.train <- xgb.DMatrix(data=sparse.model.matrix(value ~ .,data=data.train[,select.cols])[,-1],label=data.train$value)
        setinfo(sparse.train,'label_lower_bound',data.train$lower_bound)
        setinfo(sparse.train,'label_upper_bound',data.train$upper_bound)
        setinfo(sparse.train,'weight',data.train$weight)
        
        # Fit the model and make predictions
        print(paste("Evaluating model ",fold.num,"...",sep=""))
        mdl <- xgboost(data=sparse.train,booster="gbtree",max_depth=mdl.d,eta=mdl.e,
                       nthread=20,nrounds=mdl.n,objective="reg:squarederror") # Note: https://github.com/uber/causalml/issues/96
        train.est.pred <- predict(mdl,newdata=sparse.train)
        train.est.val <- data.train$value
        
        # Return predictions vs. actual values, RMSE, model, fold number, hold out fold indices, and data
        list(results=data.frame(val=train.est.val,pred=train.est.pred,row.names=rownames(data.train)),
             results.rmse=rmse(train.est.val,train.est.pred),
             model=mdl,fold=fold.num,data=data.train)
      }
      
      model.iterations[[iter]] <- cv.results
      
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
#' each iteration.
#'
#' @param model.iterations Modeling results (output of iterateModelKFCVwSG(), 
#' specifically output$model.iterations)
#' @return Collated predictions of each test fold, by iteration
#' @export
getTestFoldData <- function(model.iterations) {
  
  test.fold.data <- data.frame(fold=rep(NA,nrow(model.iterations[[1]][[1]]$data)),
                               req.cat=model.iterations[[1]][[1]]$data$req.cat,
                               actual=model.iterations[[1]][[1]]$data$orig.value,
                               row.names=rownames(model.iterations[[1]][[1]]$data))
  
  for(i in c(1:length(model.iterations))) {
    mi <- model.iterations[[i]]
    
    current.iter <- as.data.frame(matrix(nrow=0,ncol=3))
    for(j in c(1:length(mi))){
      ids <- rownames(mi[[j]]$data[mi[[j]]$data$is.test==1,])
      
      # NOTE: "value" contains initial smart guess (iter0), along with subsequent adjusted predictions (iter1-n)
      # NOTE: "prediction" contains raw predictions (pre-adjustment)
      g <- mi[[j]]$data[ids,"value"]
      # p <- mi[[j]]$results[ids,"pred"]
      f <- rep(j,length(ids))
      
      # current.iter <- rbind(current.iter,data.frame(guess=g,predict=p,fold=f,row.names=ids))
      current.iter <- rbind(current.iter,data.frame(guess=g,fold=f,row.names=ids))
      
    }
    colnames(current.iter) <- paste(c("Prediction","fold"),(i-1),sep="")
    test.fold.data <- cbind(test.fold.data,current.iter[rownames(test.fold.data),])
    test.fold.data$fold <- test.fold.data[,ncol(test.fold.data)]
    test.fold.data[,ncol(test.fold.data)] <- NULL
  }
  
  return(test.fold.data)
}


#' Plot predicted vs. actual values
#'
#' DETAILED DESCRIPTION
#'
#' @param test.fold.data Test fold predictions (output of getTestFoldData())
#' @param plot.alpha Alpha (transparency) value to use for points
#' @param n.iter How many iterations are to be plotted
#' @param reqcat.palette A character vector specifying colors to use for each 
#' requirement category; defaults are black (ENV), yellow (PHY), and cyan (EDU)
#' @param line.color A string specifying color to use for the line y=x; default 
#' is red
#' @param titles A character vector specifying the titles of each subplot; 
#' default is the column names of test.fold.data
#' @return Produces a plot (.png file)
#' @export
plotTestFolds <- function(test.fold.data,plot.alpha,n.iter,
                          reqcat.palette=c("black","yellow","cyan"),
                          line.color=c("red"),
                          titles=gsub("([a-z])([0-9])","\\1 \\2",colnames(test.fold.data))) {

  test.fold.data <- test.fold.data[,c(1:(n.iter+3))]
  reqcat.labels <- c("Environmental conditions","Physical demands","Education, training, and experience")

  plts <- list()
  
  for(i in c(4:ncol(test.fold.data))) {
    selected.cols <- test.fold.data[,c("req.cat","actual",colnames(test.fold.data)[i])]
    colnames(selected.cols) <- c("req.cat","actual","calculated")
    
    plts[[(i-3)]] <- ggplot(selected.cols) + geom_point(aes(x=calculated,y=actual,color=req.cat),alpha=plot.alpha) +
      theme_bw() + scale_colour_manual(values=reqcat.palette,labels=reqcat.labels,name="") +
      labs(title=paste(titles[i],"vs. Actual")) + xlab(titles[i]) + ylab("Actual") +
      geom_abline(intercept=0,slope=1,color=line.color) + theme(legend.position="bottom")
  }
  # return(plts)
  
  plots.grid <- ggpubr::ggarrange(plotlist=plts,nrow=ceiling(length(plts)/2),ncol=2,labels="auto",
                                  legend.grob=get_legend(plts[[1]]),legend="bottom")
  # return(plots.grid)
  
  ggsave("residPlot.png",plots.grid,width=8.5,height=(ceiling(length(plts)/2)*4 + 0.36))
}


#' Compute model blending ratio
#'
#' DETAILED DESCRIPTION
#'
#' @param model.results.soc2
#' @param model.results.soc3
#' @param convergence.iteration
#' @return 
#' @export
blendModels <- function(model.results.soc2,model.results.soc3,convergence.iteration) {
  test.folds.soc2 <- getTestFoldData(model.results.soc2$model.iterations)
  test.folds.soc3 <- getTestFoldData(model.results.soc3$model.iterations)
  conv.iter.label <- paste("Prediction",convergence.iteration,sep="")

  # Blending ratio
  blending.ratios <- vector()
  for (i in seq(0,1,0.01)) {
    lin.com <- (i*test.folds.soc2[rownames(test.folds.soc2),conv.iter.label]) + 
      ((1-i)*test.folds.soc3[rownames(test.folds.soc2),conv.iter.label])
    blending.ratios <- c(blending.ratios,rmse(test.folds.soc2$actual,lin.com))
  }
  rm(lin.com)
  blending.ratios <- data.frame(RMSE=blending.ratios,
                                SOC2.contribution=seq(0,1,0.01),
                                SOC3.contribution=rev(seq(0,1,0.01)),
                                color=0)
  blending.ratios[which.min(blending.ratios$RMSE),"color"] <- 1
  blending.ratios$color <- as.factor(blending.ratios$color)
  p.blend <- ggplot(blending.ratios,aes(x=SOC2.contribution,y=RMSE)) + 
    geom_point(aes(color=color),alpha=0.5) + scale_color_manual(values=c("grey50","red")) + 
    labs(title="RMSE of blended predictions (at convergence)") + theme_bw() + theme(legend.position="none") + 
    scale_x_continuous("SOC2 model contribution",sec.axis=sec_axis(~ . *-1 + 1,name="SOC3 model contribution"))
  ggsave(paste("kfcv-blended-rmse-iter",convergence.iteration,".png",sep=""),p.blend,width=6,height=3.75)
  
  soc2.prop <- blending.ratios[blending.ratios$color==1,"SOC2.contribution"]
  soc3.prop <- blending.ratios[blending.ratios$color==1,"SOC3.contribution"]
  
  # Blended data
  test.folds.blend <- test.folds.soc2[,1:3]
  test.folds.blend <- cbind(test.folds.blend,
                            (soc2.prop*test.folds.soc2[,4:ncol(test.folds.soc2)] + 
                               soc3.prop*test.folds.soc3[rownames(test.folds.soc2),4:ncol(test.folds.soc3)]))
  
  # Return data
  return(list(soc2.proportion=soc2.prop,
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
