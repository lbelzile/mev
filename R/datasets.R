#' Eskdalemuir Observatory daily rainfall
#'
#' This dataset contains exceedances of 30mm for daily
#' cumulated rainfall observations over the period 1970-1986.
#' These data were aggregated from hourly series.
#'
#' The station is one of the rainiest of the whole UK, with an average 1554m of cumulated rainfall per year.
#' The data consisted of 6209 daily observations, of which 4409 were non-zero.
#' Only the 93 largest observations are provided.
#' @name eskrain
#' @docType data
#' @format a vector with 93 daily cumulated rainfall measurements exceeding 30mm.
#' @source Met Office.
NULL


#' Best 200 times of women 1500m track
#'
#' 200 all-time best performance (in seconds) of women 1500-meter run.
#' @format a vector of size 200
#' @source <http://www.alltime-athletics.com/w_1500ok.htm>, accessed 14.08.2018
#' @name w1500ml
#' @docType data
NULL


#' River Nidd dataset
#'
#' The data consists of exceedances over the threshold 65 cubic meter per second of the River Nidd at Hunsingore Weir, for 35 years of data between 1934 and 1969.
#'
#' @source Natural Environment Research Council (1975). \emph{Flood Studies Report}, volume 4.  pp. 235--236.
#' @references Davison, A.C. and R.L. Smith (1990). Models for Exceedances over High Thresholds, \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{52}(3), 393--442. With discussion.
#' @format a vector of size 154
#' @name nidd
#' @seealso \code{\link[evir]{nidd.thresh}}
#' @docType data
NULL
