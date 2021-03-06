% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ors.from.csv}
\alias{ors.from.csv}
\title{Data used for original analysis}
\format{
A data frame with 97,024 rows and 16 variables describing the
following for each observation:
\describe{
\item{estimate_series_new}{series identification code}
\item{reference_year}{year of data collection}
\item{ownership_text}{exclusively "Civilian" data}
\item{industry_text}{exclusively "All Workers" data}
\item{lower_soc_code}{SOC classification code}
\item{upper_soc_code}{SOC classification code}
\item{occupation_text}{describes occupation}
\item{subcell_text}{exclusively "All Workers" data}
\item{data_element_text}{describes requirement}
\item{data_type_text}{describes requirement level}
\item{estimate_type_text}{describes the metric \emph{value} represents (mean,
mode, percent, or percentile)}
\item{estimate_type}{is \emph{value} a mean or a standard error?}
\item{unit_of_measure}{describes units of \emph{value} (days, hours, percent
of day, percentage, or pounds)}
\item{reference_group}{code specifying reference group}
\item{additive_group}{code specifying additive group}
\item{value}{numerical estimate}
}
}
\source{
See \url{https://download.bls.gov/pub/time.series/or/or.txt} and
\url{https://www.bls.gov/help/hlpforma.htm#OR} for further insight into the
format of the original data.
}
\usage{
ors.from.csv
}
\description{
Raw data that was used to conduct the original analysis described in the full
paper.
}
\details{
This dataset was derived from the raw data available on the BLS-ORS  website,
and formatted to its current state using SAS code.
}
\keyword{datasets}
