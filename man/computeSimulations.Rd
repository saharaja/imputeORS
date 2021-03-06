% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_simulation.R
\name{computeSimulations}
\alias{computeSimulations}
\title{Generate simulated data}
\usage{
computeSimulations(df)
}
\arguments{
\item{df}{Data augmented with relevant predictors and default modeling
weights (output of \code{\link[=setDefaultModelingWeights]{setDefaultModelingWeights()}})}
}
\value{
Original input, with simulated data appended
}
\description{
Rather than simply running the models on the original data, we created
simulations for each record to generate 10 unique datasets, where each
(simulated) estimate was drawn from a beta distribution derived from the
original data (mean and error estimates). This allows for generating a family
of models and a distribution of imputed values for each missing estimate,
rather than a single point estimate. This distribution can be used to
calculate a confidence level for each imputed value.
}
\seealso{
\code{\link[=stdErrMC]{stdErrMC()}}

\code{\link[=setDefaultModelingWeights]{setDefaultModelingWeights()}}
}
