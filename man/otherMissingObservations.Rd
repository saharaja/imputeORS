% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_preprocessing.R
\name{otherMissingObservations}
\alias{otherMissingObservations}
\title{Generate additional (optional) missing observations}
\usage{
otherMissingObservations(orsC)
}
\arguments{
\item{orsC}{Completed data (output of \code{\link[=fillMissingObservations]{fillMissingObservations()}})}
}
\value{
Completed data, with requested requirement/levels added
}
\description{
SSA analysts requested the addition of a few other missing observations,
based on requirement/level. These observations are appended to the data here.
}
\details{
The requested observations were as follows, listed by requirement and
constituent levels. Note that these were added to each occupation.

Humidity: CONSTANTLY

Heavy vibrations: CONSTANTLY, FREQUENTLY, SELDOM

Extreme cold: CONSTANTLY, FREQUENTLY

Climbing ladders, ropes, or scaffolds: CONSTANTLY, FREQUENTLY

Climbing ramps or stairs (work-related): CONSTANTLY, FREQUENTLY

Stooping: CONSTANTLY

Kneeling: CONSTANTLY, FREQUENTLY

Crawling: CONSTANTLY, FREQUENTLY
}
\seealso{
\code{\link[=fillMissingObservations]{fillMissingObservations()}}
}
