# Post-imputation functions


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


#' Calculate confidence intervals for predictions (missing observations only)
#' 
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @param confidence.level
#' @return 
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


#' Compute expected values for each observation
#' 
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @return 
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


#' Compute correlation of EVs, and plot relevant heatmap (grouped by SOC2 code)
#' 
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @param print.plot 
#' @return 
#' @export
correlationPlot <- function(blended.results,print.plot=TRUE) {
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
                                   width=unit(6.5,"in"),height=unit(6.5,"in"))
  p.pub <- ComplexHeatmap::Heatmap(cor.mat,name="Correlation",col=rev(brewer.pal(n=11,name="RdBu")),
                                   cluster_rows=FALSE,cluster_columns=FALSE,
                                   show_row_names=FALSE,show_column_names=FALSE,
                                   row_split=substr(rownames(cor.mat),1,2),row_gap=unit(1,"mm"),
                                   column_split=substr(colnames(cor.mat),1,2),column_gap=unit(1,"mm"),
                                   row_title_gp=gpar(fontsize=8),column_title_gp=gpar(fontsize=8),
                                   width=unit(6.5,"in"),height=unit(6.5,"in"))
  
  if (print.plot) {
    grDevices::pdf("correlation-heatmap.pdf",width=8.5,height=7)
    print(p.pub)
    grDevices::dev.off()
  }
  
  return(list(expectVal.mat=expectVal.matrix,cor.mat=cor.mat,
              plot.labeled=p.lab,plot.pub=p.pub))
}


#' Compute overlap score between two jobs, and generate plot
#' 
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @param jobA
#' @param jobB
#' @param print.plot 
#' @return 
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
  sd.ab <- sd(jobAB$overlap)
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
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @param print.plot 
#' @return 
#' @export
standardizeEVs <- function(blended.results,print.plot=TRUE) {
  # Create combined Lifting/carrying and Reaching additive groups (for plotting)
  LC <- blended.results[blended.results$requirement=="Lifting/carrying",]
  LC$additive_group <- "LC"
  R <- blended.results[blended.results$requirement=="Reaching",]
  R$additive_group <- "R"
  blended.results.LCR <- rbind(blended.results,LC,R)
  
  # Standardize expected values by requirement
  expected.vals <- imputeORS::computeEVs(blended.results.LCR)
  std.expectVal <- apply(expected.vals,1,robustHD::standardize,centerFun=mean,scaleFun=sd)
  
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
#' DETAILED DESCRIPTION
#'
#' @param blended.results
#' @param jobA
#' @param jobB
#' @param print.plot 
#' @return 
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

