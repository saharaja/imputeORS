

#' Compute sums of occupational groups
#' 
#' DETAILED DESCRIPTION
#'
#' @param ors.data
#' @param column.name
#' @return 
#' @export
checkSums <- function(ors.data,column.name="prediction") {
  pct.sums <- vector()
  for (occ in levels(ors.data$occupation_text)) {
    for (adg in levels(ors.data$additive_group)) {
      current.occ.adg <- ors.data[as.character(ors.data$occupation_text)==occ & 
                                    as.character(ors.data$additive_group)==adg,]
      
      if (nrow(current.occ.adg)==0) {
        next
      }
      names.temp <- c(names(pct.sums),paste(occ,adg,sep="-"))
      pct.sums <- c(pct.sums,sum(current.occ.adg[,column.name],na.rm=TRUE))
      names(pct.sums) <- names.temp
      rm(names.temp)
    }
  }
  return(pct.sums)
}


#' Get convergence predictions (one per simulation) and overall mean prediction
#' 
#' DETAILED DESCRIPTION
#'
#' @param simulations
#' @param iteration
#' @param ors.data.sims
#' @return 
#' @export
getMeanPredictions <- function(simulations,iteration,ors.data.sims) {
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


#' Calculate confidence intervals for predictions (missing observations only)
#' 
#' DETAILED DESCRIPTION
#'
#' @param model.summary.average
#' @param confidence.level
#' @return 
#' @export
predictCIs <- function(model.summary.average,confidence.level=0.95) {
  # Get predictions from relevant iteration of each simulation to create distributions
  # Use MISSING OBSERVATIONS only
  p <- model.summary.average[model.summary.average$known.val==0,
                             grep(paste("sim[0-9]*\\.prediction",sep=""),colnames(model.summary.average))]
  
  # Get CIs
  returnConfInt <- function(distrib,con.lvl=confidence.level) {
    tryCatch(
      {
        return(c(t.test(distrib,conf.level=con.lvl)$conf.int[1:2],mean(distrib),0))
        # return(c(confint(lm(unlist(distrib)~1),level=confidence.level)[1:2],mean(distrib),0))   # equivalent method
      } , error=function(e) {
        print(e)
        return(c(NA,NA,mean(distrib),1))
        # return(c(mean(distrib),mean(distrib),1))
      }
    )
  }
  # p <- p[1:1000,]  # for debugging
  ci <- apply(p,1,returnConfInt)
  ci <- as.data.frame(t(ci))
  colnames(ci) <- c("lower","upper","mean.prediction","error.flag")
  
  return(list(distributions=p,CIs=ci))
}


## FUNCTION to compute expected values
#' SUMMARY
#' 
#' DETAILED DESCRIPTION
#'
#' @param simulations
#' @return 
#' @export
computeEVs <- function(simulations) {
  # Get mean predictions
  data <- getMeanPredictions(simulations)
  
  # Compute score
  data$score <- data$mean.prediction * data$frequency * data$intensity
  
  # Reformat and get aggregate scores to get expected values
  expectVal.matrix <- reshape2::dcast(data,additive_group~occupation_text,value.var="score",fun.aggregate=sum)
  rownames(expectVal.matrix) <- expectVal.matrix$additive_group
  expectVal.matrix$additive_group <- NULL
  
  return(expectVal.matrix)
}


## FUNCTION to compute correlation of expected values and produce relevant heatmap (organized by SOC2 code)
#' SUMMARY
#' 
#' DETAILED DESCRIPTION
#'
#' @param simulations
#' @return 
#' @export
correlationPlot <- function(simulations) {
  # Compute expected values
  expectVal.matrix <- computeEVs(simulations)
  
  # Correlation
  cor.mat <- cor(expectVal.matrix)
  
  # Get SOC2 codes
  get.soc2.codes <- data[,c("upSOC2","occupation_text")]
  get.soc2.codes <- distinct(get.soc2.codes)
  rownames(get.soc2.codes) <- get.soc2.codes$occupation_text
  rownames(cor.mat) <- paste(get.soc2.codes[rownames(cor.mat),"upSOC2"],rownames(cor.mat),sep="-")
  colnames(cor.mat) <- paste(get.soc2.codes[colnames(cor.mat),"upSOC2"],colnames(cor.mat),sep="-")
  
  # Calculate positions of gaps in heatmap
  p.gaps <- cumsum(table(substr(rownames(cor.mat),1,2)))
  
  # Color palette
  p.color <- rev(colorRampPalette(brewer.pal(n=11,name="RdBu"))(650))
  
  # Heatmaps (with and without row/column names)
  p.lab <- Heatmap(cor.mat,name="Correlation",col=rev(brewer.pal(n=11,name="RdBu")),
                   cluster_rows=FALSE,cluster_columns=FALSE,
                   row_split=substr(rownames(cor.mat),1,2),row_gap=unit(1,"mm"),
                   column_split=substr(colnames(cor.mat),1,2),column_gap=unit(1,"mm"),
                   row_title_gp=gpar(fontsize=8),column_title_gp=gpar(fontsize=8),
                   width=unit(6.5,"in"),height=unit(6.5,"in"))
  p.pub <- Heatmap(cor.mat,name="Correlation",col=rev(brewer.pal(n=11,name="RdBu")),
                   cluster_rows=FALSE,cluster_columns=FALSE,
                   show_row_names=FALSE,show_column_names=FALSE,
                   row_split=substr(rownames(cor.mat),1,2),row_gap=unit(1,"mm"),
                   column_split=substr(colnames(cor.mat),1,2),column_gap=unit(1,"mm"),
                   row_title_gp=gpar(fontsize=8),column_title_gp=gpar(fontsize=8),
                   width=unit(6.5,"in"),height=unit(6.5,"in"))
  
  return(list(expectVal.mat=expectVal.matrix,cor.mat=cor.mat,
              plot.labeled=p.lab,plot.pub=p.pub))
}


## FUNCTION to compute overlap for two jobs A and B and generate plot
#' SUMMARY
#' 
#' DETAILED DESCRIPTION
#'
#' @param simulations
#' @param jobA
#' @param jobB
#' @return 
#' @export
computeOverlap <- function(simulations,jobA,jobB) {
  # Get mean predictions
  data <- getMeanPredictions(simulations)
  
  # Filter by jobs
  data.a <- droplevels(data[as.character(data$occupation_text)==jobA,
                            c("additive_group","requirement","occupation_text","frequency",
                              "intensity","upSOC2","upSOC3","upSOC4","known.val","mean.prediction")])
  data.b <- droplevels(data[as.character(data$occupation_text)==jobB,
                            c("additive_group","requirement","occupation_text","frequency",
                              "intensity","upSOC2","upSOC3","upSOC4","known.val","mean.prediction")])
  data.ab <- merge(data.a,data.b,
                   by=c("additive_group","requirement","frequency","intensity"),
                   sort=FALSE,suffixes=c(".a",".b"))
  rm(data.a)
  rm(data.b)
  
  # Compute overlap
  jobAB <- distinct(data.ab[,c("additive_group","requirement")])
  jobAB$overlap <- NA
  rownames(jobAB) <- jobAB$additive_group
  data.ab$o <- data.ab$mean.prediction.a * data.ab$mean.prediction.b
  for (i in levels(jobAB$additive_group)) {
    current.adg <- data.ab[which(as.character(data.ab$additive_group)==i),]
    jobAB[i,"overlap"] <- sum(current.adg$o)
  }
  
  # Plot
  m.ab <- mean(jobAB$overlap)
  sd.ab <- sd(jobAB$overlap)
  p.ab <- ggplot(jobAB) + geom_bar(aes(x=as.factor(additive_group),y=overlap),stat="identity") + 
    scale_x_discrete(name="Additive Group (Requirement)") + 
    scale_y_continuous(name="Overlap",breaks=seq(0,1,0.25),limits=c(0,1)) +
    geom_hline(yintercept=m.ab,linetype="dashed",color="red") +
    geom_hline(yintercept=m.ab+sd.ab,linetype="dashed",color="grey") +
    geom_hline(yintercept=m.ab-sd.ab,linetype="dashed",color="grey") +
    ggtitle(label=paste(jobA," vs. ",jobB,sep="")) + theme_bw()
  
  return(list(jobs.data=data.ab,
              overlap=jobAB,
              plot=p.ab))
}


### FUNCTION to compute 5-iteration rolling average MAE and STD (using missing observations only)
##' SUMMARY
##'
##' DETAILED DESCRIPTION
##'
##' @param simulation
##' @return
##' @export
# modelConvergence <- function(simulation) {
#   avg.mae <- c(NA,NA,NA,NA)
#   stddev <- c(NA,NA,NA,NA)
#   for (i in c(6:length(simulation))) {
#     select.missing <- which(simulation$iteration0$known.val!=1)
#     # print(length(select.missing))
#     maes <- c(mae(simulation[[i-4]]$results[select.missing,"previous.prediction"],
#                   simulation[[i-4]]$results[select.missing,"current.prediction"]),
#               mae(simulation[[i-3]]$results[select.missing,"previous.prediction"],
#                   simulation[[i-3]]$results[select.missing,"current.prediction"]),
#               mae(simulation[[i-2]]$results[select.missing,"previous.prediction"],
#                   simulation[[i-2]]$results[select.missing,"current.prediction"]),
#               mae(simulation[[i-1]]$results[select.missing,"previous.prediction"],
#                   simulation[[i-1]]$results[select.missing,"current.prediction"]),
#               mae(simulation[[i-0]]$results[select.missing,"previous.prediction"],
#                   simulation[[i-0]]$results[select.missing,"current.prediction"]))
#     avg.mae <- c (avg.mae,mean(maes))
#     stddev <- c(stddev,sd(maes))
# 
#   }
#   avg.mae.diff <- c(NA,abs(avg.mae[2:length(avg.mae)]-avg.mae[1:(length(avg.mae)-1)]))
#   return(list(rolling.mae.average=avg.mae,
#               rolling.mae.diff=avg.mae.diff,
#               standard.dev.mae=stddev,
#               convergence.iteration=which(avg.mae.diff<0.005 & stddev<0.01)[1]))
# }