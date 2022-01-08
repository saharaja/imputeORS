# Documentation of external data (imputeORS/data)

#' Data used for original analysis
#'
#' Raw data that was used to conduct the original analysis described in the full
#' paper. 
#' 
#' This dataset was derived from the raw data available on the BLS-ORS  website, 
#' and formatted to its current state using SAS code.
#'
#' @format A data frame with 97,024 rows and 16 variables describing the 
#' following for each observation:
#' \describe{
#'   \item{estimate_series_new}{series identification code}
#'   \item{reference_year}{year of data collection}
#'   \item{ownership_text}{exclusively "Civilian" data}
#'   \item{industry_text}{exclusively "All Workers" data}
#'   \item{lower_soc_code}{SOC classification code}
#'   \item{upper_soc_code}{SOC classification code}
#'   \item{occupation_text}{describes occupation}
#'   \item{subcell_text}{exclusively "All Workers" data}
#'   \item{data_element_text}{describes requirement}
#'   \item{data_type_text}{describes requirement level}
#'   \item{estimate_type_text}{describes the metric [[value]] represents (mean, 
#'   mode, percent, or percentile)}
#'   \item{estimate_type}{is [[value]] a mean or a standard error?}
#'   \item{unit_of_measure}{describes units of [[value]] (days, hours, percent 
#'   of day, percentage, or pounds)}
#'   \item{reference_group}{code specifying reference group}
#'   \item{additive_group}{code specifying additive group}
#'   \item{value}{numerical estimate}
#' }
#' @source See \url{https://download.bls.gov/pub/time.series/or/or.txt} and 
#' \url{https://www.bls.gov/help/hlpforma.htm#OR} for further insight into the 
#' format of the original data.
"ors.from.csv"
#"estim_ssa_output_2018"



#' Predictor mapping
#'
#' A data frame describing the refining of textual (categorical) information in 
#' the ORS dataset, along with the mapping of this refined information to three 
#' new predictors: Frequency (numerical), Intensity (numerical), and 
#' Requirement_category (categorical). This mapping is used to define the 
#' predictor columns for the data in preparation for modeling.
#'
#' @format A data frame with 299 rows and 7 variables:
#' \describe{
#'   \item{data_element_text}{original data field describing the requirement}
#'   \item{data_type_text}{original data field describing the requirement level}
#'   \item{additive_group}{original data field specifying the additive group}
#'   \item{Requirement}{new, refined data field describing the requirement}
#'   \item{Frequency}{new data field describing the requirement level (1 of 2)}
#'   \item{Intensity}{new data field describing the requirement level (2 of 2)}
#'   \item{Requirement_category}{new data field describing the requirement 
#'   category}
#' }
#' @source See \url{https://www.bls.gov/ors/information-for-survey-participants/pdf/occupational-requirements-survey-collection-manual-082020.pdf}
#' for the guidance used to arrive at Frequency/Intensity/Requirement_category 
#' mappings.
#' @seealso [getPredictors()]
"predictors.data"
#"data_mapping_for_predictors"

