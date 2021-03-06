% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_preprocessing.R
\name{syntheticAddGroups}
\alias{syntheticAddGroups}
\title{Generate synthetic additive groups}
\usage{
syntheticAddGroups(orsRaw)
}
\arguments{
\item{orsRaw}{Original data (read in from CSV file)}
}
\value{
Data edited for new additive groups, and containing only relevant
records (mean and standard error estimates) and columns
}
\description{
Five synthetic additive groups were manually constructed based on the advice
of SSA analysts, using reference_group numbers:
(11) Literacy required,
(89) Pre-employment training: Certification,
(90) Pre-employment training: License,
(91) Pre-employment training: Educational Certificate, and
(78) Sitting or standing/walking.
}
\seealso{
\link{ors.from.csv} for original data used in our analysis
}
