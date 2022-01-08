# Post-imputation functions

#' Compute sums of occupational groups
#' 
#' This is an accessory function to ensure that the sum of the estimates for 
#' each of the occupational groups is 1. If the input data is the original (i.e.
#' incomplete) data, then all the sums should be <=1. If the input data is the 
#' results of smart guessing or imputing the missing values, the sums should all
#' be (approximately) 1.
#'
#' @param ors.data Data frame containing ORS data, including estimates
#' @param column.name Name of column containing estimates; default is 
#' "prediction"
#' @return A vector whose length is the number of occupational groups present in
#' the data, with each entry being the sum of the estimates associated with a 
#' given occupational group
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


#' Calculate confidence intervals for final predictions (missing values only)
#' 
#' Our procedure involves creating a set of simulated data based on the original 
#' known data. Thus, for each observation we end up with a distribution of 
#' predictions, rather than a single point estimate. This allows us to calculate 
#' confidence intervals around each prediction. This CI calculation is done only
#' for observations with missing values (i.e., those that were imputed).
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @param confidence.level Level at which to calculate confidence intervals; 
#' default is 0.95
#' @return A list of length two, containing a data frame describing distribution
#' of predictions (one prediction per simulation) for each observation that was 
#' imputed, and another data frame describing the confidence interval calculated
#' for each of these observations; note that the latter includes a column called 
#' "error.flag" that indicates whether there was an error in computing the CI 
#' (the detailed errors are printed during execution)
#' @seealso [blendImputations()]
#' @seealso [stats::t.test()]
#' @export
computeCIs <- function(blended.results,confidence.level=0.95) {
  # Get predictions from relevant iteration of each simulation to create distributions
  # Use MISSING OBSERVATIONS only
  p <- blended.results[blended.results$known.val==0,
                             grep(paste("sim[0-9]*\\.prediction",sep=""),colnames(blended.results))]
  
  # Get CIs
  returnConfInt <- function(distrib,con.lvl=confidence.level) {
    tryCatch(
      {
        return(c(stats::t.test(distrib,conf.level=con.lvl)$conf.int[1:2],mean(distrib),0))
        # return(c(stats::confint(stats::lm(unlist(distrib)~1),level=confidence.level)[1:2],mean(distrib),0))   # equivalent method
      } , error=function(e) {
        print(e)
        return(c(NA,NA,mean(distrib),1))
      }
    )
  }
  # p <- p[1:1000,]  # for debugging
  ci <- apply(p,1,returnConfInt)
  ci <- as.data.frame(t(ci))
  colnames(ci) <- c("lower","upper","mean.prediction","error.flag")
  
  return(list(distributions=p,CIs=ci))
}


#' Compute expected values for each observation
#' 
#' To compare occupations across requirements, we developed an expected value
#' measure we called the "Expected Level of Effort" (ELE). This measure is a 
#' weighted average of the frequency and intensity times the population estimate
#' for the various requirements. A low frequency/low intensity/low population 
#' estimate results in a low level of effort, and the converse for high.
#' 
#' For each occupational group, we calculate ELE as an expected value of 
#' frequency times intensity as follows, where \eqn{\mu_j} is the mean population 
#' prediction across all the simulations for the \eqn{j}th observation, and 
#' \eqn{F_j} and \eqn{I_j} are the frequency and intensity of the \eqn{j}th 
#' observation, respectively:
#' 
#' \deqn{E=\sum(\mu_j * F_j * I_j)}
#'  
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @return A data frame containing ELEs of each occupational group, arranged 
#' with requirement (additive groups) as rows and occupation as columns
#' @seealso [blendImputations()]
#' @export
computeEVs <- function(blended.results) {
  
  # Compute score
  blended.results$score <- blended.results$mean.prediction * blended.results$frequency * blended.results$intensity
  
  # Reformat and get aggregate scores to get expected values
  expectVal.matrix <- reshape2::dcast(blended.results,additive_group~occupation_text,
                                      value.var="score",fun.aggregate=sum)
  rownames(expectVal.matrix) <- expectVal.matrix$additive_group
  expectVal.matrix$additive_group <- NULL
  
  return(expectVal.matrix)
}


#' Compute correlation of ELEs, and plot relevant heatmap (grouped by SOC2 code)
#' 
#' Produces a heatmap displaying Pearson correlation between expected values 
#' (ELEs) by occupation. The heatmap is arranged according to SOC2 code, so that 
#' similar occupations are grouped together.
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @param plot.dim Length of each side of the plot (in inches); default is 6.5
#' @return A list of length four, containing a data frame of expected values, a 
#' correlation matrix, and two plot objects (one where the individual 
#' occupations are labeled in addition to the SOC2 group labels, and one where 
#' they are not); optionally produces a PDF of the heatmap
#' @seealso [blendImputations()]
#' @seealso [computeEVs()]
#' @seealso [stats::cor()]
#' @export
correlationPlot <- function(blended.results,plot.dim=6.5,print.plot=TRUE) {
  # Compute expected values
  expectVal.matrix <- imputeORS::computeEVs(blended.results)
  
  # Correlation
  cor.mat <- stats::cor(expectVal.matrix)
  
  # Get SOC2 codes
  get.soc2.codes <- blended.results[,c("upSOC2","occupation_text")]
  get.soc2.codes <- dplyr::distinct(get.soc2.codes)
  rownames(get.soc2.codes) <- get.soc2.codes$occupation_text
  rownames(cor.mat) <- paste(get.soc2.codes[rownames(cor.mat),"upSOC2"],rownames(cor.mat),sep="-")
  colnames(cor.mat) <- paste(get.soc2.codes[colnames(cor.mat),"upSOC2"],colnames(cor.mat),sep="-")
  
  # Calculate positions of gaps in heatmap
  p.gaps <- cumsum(table(substr(rownames(cor.mat),1,2)))
  
  # Color palette
  p.color <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(n=11,name="RdBu"))(650))
  
  # Heatmaps (with and without row/column names)
  p.lab <- ComplexHeatmap::Heatmap(cor.mat,name="Correlation",col=rev(brewer.pal(n=11,name="RdBu")),
                                   cluster_rows=FALSE,cluster_columns=FALSE,
                                   row_split=substr(rownames(cor.mat),1,2),row_gap=unit(1,"mm"),
                                   column_split=substr(colnames(cor.mat),1,2),column_gap=unit(1,"mm"),
                                   row_title_gp=gpar(fontsize=8),column_title_gp=gpar(fontsize=8),
                                   width=unit(plot.dim,"in"),height=unit(plot.dim,"in"))
  p.pub <- ComplexHeatmap::Heatmap(cor.mat,name="Correlation",col=rev(brewer.pal(n=11,name="RdBu")),
                                   cluster_rows=FALSE,cluster_columns=FALSE,
                                   show_row_names=FALSE,show_column_names=FALSE,
                                   row_split=substr(rownames(cor.mat),1,2),row_gap=unit(1,"mm"),
                                   column_split=substr(colnames(cor.mat),1,2),column_gap=unit(1,"mm"),
                                   row_title_gp=gpar(fontsize=8),column_title_gp=gpar(fontsize=8),
                                   width=unit(plot.dim,"in"),height=unit(plot.dim,"in"))
  
  if (print.plot) {
    grDevices::pdf("correlation-heatmap.pdf",width=(plot.dim+2),height=(plot.dim+0.5))
    print(p.pub)
    grDevices::dev.off()
  }
  
  return(list(expectVal.mat=expectVal.matrix,cor.mat=cor.mat,
              plot.labeled=p.lab,plot.pub=p.pub))
}


#' Compute overlap score between two jobs, and generate plot
#' 
#' Thanks to imputation, we arrive at a full distribution of population 
#' percentages for all jobs, across all requirements. Using this information, we 
#' can devise a way to measure the requirement "overlap" of any two occupations. 
#' We do this by taking the product of their estimates to yield a value on the 
#' interval \[0,1\]. This produces a way to measure similarity between different
#' occupations.
#' 
#' If \eqn{\omega_1ir} is the population mean of the \eqn{i}th level of the 
#' \eqn{r}th requirement for Job 1 (average of simulation predictions for 
#' missing values, and actual value for known estimates), and \eqn{\omega_2ir} 
#' is the same for Job 2, then we say that the overlap of \eqn{r}th requirement 
#' (\eqn{O_r}) for these two jobs is:
#' 
#' \deqn{O_r = \sum(\omega_1ir * \omega_2ir)}
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @param jobA First job to compare
#' @param jobB Second job to compare
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @return A list of length three, containing a data frame with all the data 
#' pertaining to the specified occupations, a data frame of the overlap values
#' for these occupations by requirement, and a plot object displaying these 
#' overlap values; optionally produces a plot (.png file) of the overlap
#' @seealso [blendImputations()]
#' @export
computeOverlap <- function(blended.results,jobA,jobB,print.plot=TRUE) {
  # Filter data by specified jobs
  data.a <- droplevels(blended.results[as.character(blended.results$occupation_text)==jobA,
                            c("additive_group","requirement","occupation_text","frequency",
                              "intensity","upSOC2","upSOC3","upSOC4","known.val","mean.prediction")])
  data.b <- droplevels(blended.results[as.character(blended.results$occupation_text)==jobB,
                            c("additive_group","requirement","occupation_text","frequency",
                              "intensity","upSOC2","upSOC3","upSOC4","known.val","mean.prediction")])
  data.ab <- merge(data.a,data.b,
                   by=c("additive_group","requirement","frequency","intensity"),
                   sort=FALSE,suffixes=c(".a",".b"))
  rm(data.a)
  rm(data.b)
  
  # Compute overlap
  jobAB <- dplyr::distinct(data.ab[,c("additive_group","requirement")])
  jobAB$overlap <- NA
  rownames(jobAB) <- jobAB$additive_group
  data.ab$o <- data.ab$mean.prediction.a * data.ab$mean.prediction.b
  for (i in levels(jobAB$additive_group)) {
    current.adg <- data.ab[which(as.character(data.ab$additive_group)==i),]
    jobAB[i,"overlap"] <- sum(current.adg$o)
  }
  
  # Plot
  m.ab <- mean(jobAB$overlap)
  sd.ab <- stats::sd(jobAB$overlap)
  p.ab <- ggplot2::ggplot(jobAB) + ggplot2::geom_bar(aes(x=as.factor(additive_group),y=overlap),stat="identity") + 
    ggplot2::scale_x_discrete(name="Additive Group (Requirement)") + 
    ggplot2::scale_y_continuous(name="Overlap",breaks=seq(0,1,0.25),limits=c(0,1)) +
    ggplot2::geom_hline(yintercept=m.ab,linetype="dashed",color="red") +
    ggplot2::geom_hline(yintercept=m.ab+sd.ab,linetype="dashed",color="grey") +
    ggplot2::geom_hline(yintercept=m.ab-sd.ab,linetype="dashed",color="grey") +
    ggplot2::ggtitle(label=paste(jobA," vs. ",jobB,sep="")) + ggplot2::theme_bw()
  
  if (print.plot) {
    ggplot2::ggsave(paste("overlap_",jobA,"-v-",jobB,".png",sep=""),plot=p.ab,width=10,height=4,units="in")
  }
  
  return(list(jobs.data=data.ab,
              overlap=jobAB,
              plot=p.ab))
}


#' Standardize expected values, and generate plot
#' 
#' The ELEs calculated by [computeEVs()] necessarily differ vastly in their 
#' range and distributions depending on the requirement. Standardizing them 
#' makes them more readily comparable.
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @param combine.adgs Option to combine the four Lifting/carrying requirements
#' and the two Reaching requirements into two aggregate requirements (LC and R, 
#' respectively); default is TRUE
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @return A list of length three, containing a data frame of expected values, a
#' matrix of standardized EVs, and a plot object of the standardized EVs; 
#' optionally produces a plot (.png file) of the standardized EVs
#' @seealso [blendImputations()]
#' @seealso [computeEVs()]
#' @seealso [robustHD::standardize()]
#' @export
standardizeEVs <- function(blended.results,combine.adgs=TRUE,print.plot=TRUE) {
  # Create combined Lifting/carrying and Reaching additive groups (for plotting)
  if (combine.adgs) {
    LC <- blended.results[blended.results$requirement=="Lifting/carrying",]
    LC$additive_group <- "LC"
    R <- blended.results[blended.results$requirement=="Reaching",]
    R$additive_group <- "R"
    blended.results <- rbind(blended.results,LC,R)
  }
  
  # Standardize expected values by requirement
  expected.vals <- imputeORS::computeEVs(blended.results)
  std.expectVal <- apply(expected.vals,1,robustHD::standardize,centerFun=mean,scaleFun=stats::sd)
  
  # Plot standardized expected values
  p <- ggplot2::ggplot(reshape2::melt(std.expectVal),aes(x=as.factor(Var2),y=value)) + 
    ggplot2::geom_point(position=position_jitter(width=0.4),color="green",alpha=0.05) + 
    ggplot2::geom_boxplot(outlier.colour="blue",alpha=0.6) + 
    ggplot2::geom_point(y=0,color="red") + 
    ggplot2::scale_x_discrete(name="Additive Group (Requirement)") + 
    ggplot2::scale_y_continuous(name="Standardized Expected Value") + 
    ggplot2::ggtitle(label="Standardized Expected Values, by Requirement") + 
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
  
  if (print.plot) {
    ggplot2::ggsave("std-EV-by-req.png",p,width=10,height=4,units="in")
  }
 
  return(list(expected.vals=expected.vals,
              std.expectVal=std.expectVal,
              stdEV.plot=p)) 
}


#' Calculate differences in standardized EVs between two jobs, and generate plot
#' 
#' Like the overlap metric, the standardized EVs can be used to compare two 
#' occupations, providing a measure of the similarity between jobs.
#'
#' @param blended.results Blended predictions from imputation models, calculated
#' at convergence iterations and blending proportions computed by 
#' [computeBlendingRatio()] (output of [blendImputations()])
#' @param jobA First job to compare
#' @param jobB Second job to compare
#' @param print.plot Should plot file (.png) be generated; default is TRUE
#' (create plot file)
#' @return A list of two objects, containing a data frame of the differences in 
#' standardized EV by requirement, and a plot object of these differences; 
#' optionally produces a plot (.png file) of the differences in standardized EVs
#' @seealso [blendImputations()]
#' @seealso [standardizeEVs()]
#' @export
computeStdEVdiff <- function(blended.results,jobA,jobB,print.plot=TRUE) {
  # Get standardized EVs
  stdEV <- imputeORS::standardizeEVs(blended.results,print.plot=FALSE)[[2]]
  
  # Compute differences in standardized EVs
  diff.std.ev <- stdEV[jobA,]-stdEV[jobB,]
  diff.std.ev <- data.frame(diff=diff.std.ev,adg=names(diff.std.ev))
  
  # Generate/plot
  p <- ggplot2::ggplot(diff.std.ev) + ggplot2::geom_bar(aes(y=diff,x=as.factor(adg)),stat="identity") + 
    ggplot2::scale_x_discrete(name="Additive Group (Requirement)") + 
    ggplot2::scale_y_continuous(name="Standardized Difference",breaks=seq(-7,7,1),limits=c(-7,7)) +
    ggplot2::ggtitle(label=paste(jobA," vs. ",jobB,sep="")) + ggplot2::theme_bw()
  if (print.plot) {
    ggplot2::ggsave(paste("stdEVdiff_",jobA,"-v-",jobB,".png",sep=""),plot=p,width=10,height=4,units="in")
  }
  
  return(list(diff.vals=diff.std.ev,diff.plot=p))
}

