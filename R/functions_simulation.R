# Helper functions

#' Estimate beta distribution parameters
#'
#' Generates parameters (alpha and beta) to define a beta distribution based on 
#' the input data (mean and variance).
#'
#' @param mu Mean of distribution
#' @param var Variance (error) of distribution
#' @return Beta distribution parameters
#' @export
estBetaParams <- function(mu,var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#' Draw sample from specified beta distribution
#'
#' This function generates a sample of n values derived from a beta distribution 
#' specified by a mean and associated variance (error). The beta distribution is
#' generated using estBetaParams(). The default sample size is 10.
#'
#' @param mu Mean of distribution
#' @param stdev Variance (error) of distribution
#' @param n Number of values to draw; default is 10
#' @return A sample of size n
#' @export
stdErrMC<- function(mu,stdev,n=10){
  betaCoefs=estBetaParams(mu,stdev*stdev)
  sample=rbeta(n,betaCoefs$alpha,betaCoefs$beta)
  sample=sample*stdev/stats::sd(sample)
  sample=sample-mean(sample)+mu
  # if(anyNA(sample)){
  #          print(mu)
  #          print(stdev)
  #          print(betaCoefs)
  # }
  return(sample)
}


# Main function 

#' Generate simulated data
#'
#' Rather than simply running the models on the original data, we created 
#' simulations for each record to generate 10 unique datasets, where each 
#' (simulated) estimate was drawn from a beta distribution derived from the 
#' original data (mean and error estimates). This allows for generating a family
#' of models and a distribution of imputed values for each missing estimate, 
#' rather than a single point estimate. This distribution can be used to 
#' calculate a confidence level for each imputed value.
#'
#' @param df Data augmented with relevant predictors and default modeling 
#' weights (output of setDefaultModelingWeights())
#' @return Original input, with simulated data appended
#' @export
computeSimulations <-function(df){
  inCount=0
  outCount=0
  rownames(df)=NULL
  # sampledVals=df[,c("occupation_text","additive_group","requirement","frequency","intensity","value","std.error")]
  sampledVals <- df
  sampledVals[,c("valsim1","valsim2","valsim3","valsim4","valsim5",
                 "valsim6","valsim7","valsim8","valsim9","valsim10")] <- NA
  # print(head(sampledVals))
  
  if (max(df$value,na.rm=TRUE)==1) {
    df$value <- df$value*100
    df$std.error <- df$std.error*100
  }
  
  jobs=unique(df$occupation_text)
  for(job in jobs){
    jobData=df[df$occupation_text==job,]
    
    ag=unique(jobData$additive_group)
    for(g in ag){
      agData=jobData[jobData$additive_group==g,]
      #print(agData)
      iMax=which(agData$value==max(agData$value,na.rm = T))
      if(length(iMax)==0){
        # jump to next as values already set to NaN
        print("No non NaN")
        next
      }
      if(length(iMax)>1){print("iMax>1");print(iMax);iMax=c(iMax[1])}
      #print(iMax)
      iExist=which(!is.na(agData$value))
      #stddev=sqrt(10)*agData$std.error[iMax]/100
      stddev=agData$std.error[iMax]/100
      val=agData$value[iMax]/100
      
      if(stddev*stddev>(1-val)*val){stddev=sqrt((1-val)*val*.98);print("stddev over limit");outCount=outCount+1} 
      else { inCount=inCount+1}
      
      
      if(val>=.975){
        mainsample=rep(val,10)
      } else if(val<.015){
        mainsample=rep(val,10)
      } else if(stddev<=0.006){
        mainsample=rep(val,10)
      } else{
        mainsample=stdErrMC(val,stddev,n=10)
        # print(mainsample)
        # print(val)
        # print(stddev)
        mainsample=pmin(mainsample,1)
        mainsample=pmax(mainsample,0)
        if(any(mainsample<0)) print(mainsample)

      }
      if(anyNA(mainsample))stop()
      print(head(sampledVals))
      print(mainsample)
      sampledVals[agData$occupation_text[iMax]==sampledVals$occupation_text &
                    agData$requirement[iMax]==sampledVals$requirement &
                    agData$frequency[iMax]==sampledVals$frequency &
                    agData$intensity[iMax]==sampledVals$intensity &
                    agData$additive_group[iMax]==sampledVals$additive_group,
                  c((ncol(sampledVals)-9):ncol(sampledVals))]=t(mainsample[1:10])

      #print(mainsample)
      
      jExist=iExist[!iExist%in%iMax ]
      if(length(jExist)==0){next} 
      else{
        for( j in jExist){
          val2=agData$value[j]/100
          stddev2=max(stddev*val2/val,agData$std.error[j]/100)
          if(val2==0){
            smallpop=rep(0.0,10)
          } else{
            smallpop=(mean(mainsample)-mainsample)*val2/(1-val)+val2
            smallpop=pmin(smallpop,1)
            smallpop=pmax(smallpop,0)
          }
          if(any(smallpop<0)){print(mainsample);print(smallpop);print(val2);print(val);stop()}
          
          sampledVals[agData$occupation_text[j]==sampledVals$occupation_text &
                        agData$requirement[j]==sampledVals$requirement &
                        agData$frequency[j]==sampledVals$frequency &
                        agData$intensity[j]==sampledVals$intensity &
                        agData$additive_group[j]==sampledVals$additive_group
                      ,c((ncol(sampledVals)-9):ncol(sampledVals))]=t(smallpop[1:10])
          #print(smallpop)
        }
      }
    }
    
  }
  #sapply(df,getPerturbedAG)
  print(inCount)
  print(outCount)
  
  rownames(sampledVals) <-  gsub(" ","",paste(sampledVals$occupation_text,sampledVals$additive_group,sampledVals$data_type_text))
  sampledVals
}


