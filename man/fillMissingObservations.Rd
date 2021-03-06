% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_preprocessing.R
\name{fillMissingObservations}
\alias{fillMissingObservations}
\title{Generate missing observations}
\usage{
fillMissingObservations(ors)
}
\arguments{
\item{ors}{Data frame of relevant observations (output of \code{\link[=syntheticAddGroups]{syntheticAddGroups()}})}
}
\value{
Data with newly generated observations added (completed data)
}
\description{
Observations that were missing in the data (i.e., no associated record or
estimate) are generated here.
This is done by "completing" the data using tidyr::complete.
}
\seealso{
\code{\link[=syntheticAddGroups]{syntheticAddGroups()}}

\code{\link[tidyr:complete]{tidyr::complete()}}
}
