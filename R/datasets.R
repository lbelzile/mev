#' Eskdalemuir Observatory Daily Rainfall
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

#' Venice Sea Levels
#'
#' The \code{venice} data contains the 10 largest yearly sea levels (in cm)
#' from 1887 until 2017. Only the yearly maximum is available for 1922
#' and the six largest observations for 1936.
#'
#' @format a data frame with 131 rows and 11 columns containing the year of the measurement (first column)
#' and ordered 10-largest yearly observations, reported in decreasing order from largest (\code{r1}) to smallest (\code{r10}).
#'
#' @note Smith (1986) notes that the annual maxima seems to fluctuate around a constant sea level
#' up to 1930 or so, after which there is potential linear trend. Records of threshold exceedances above
#' 80 cm (reported on the website) indicate that observations are temporally clustered.
#'
#' The observations from 1931 until 1981 can be found in
#' Table 1 in Smith (1986), who reported data from Pirazzoli (1982).
#' The values from 1983 until 2017 were extracted by Anthony Davison from the City
#' of Venice website (accessed October 2018) and are licensed under the CC BY-NC-SA 3.0 license.
#' The Venice City website indicates
#' that later measurements were recorded by an instrument located in Punta Salute.
#'
#'
#' @references Smith, R. L. (1986) Extreme value theory based on the \emph{r}
#' largest annual events. \emph{Journal of Hydrology} \bold{86}, 27–43.
#' @references Pirazzoli, P., 1982. Maree estreme a Venezia (periodo 1872-1981). \emph{Acqua Aria} \bold{10}, 1023-1039.
#' @references  Coles, S. G. (2001) \emph{An Introduction to Statistical Modelling of Extreme Values}. London: Springer.
#' @source City of Venice, Historical archive <http://archive.comune.venezia.it/flex/cm/pages/ServeBLOB.php/L/EN/IDPagina/3045>. Last accessed October 2018.
#' @name venice
#' @seealso \code{\link[ismev]{venice}}
#' @docType data
NULL


#' Best 200 times of Women 1500m Track
#'
#' 200 all-time best performance (in seconds) of women 1500-meter run.
#' @format a vector of size 200
#' @source <http://www.alltime-athletics.com/w_1500ok.htm>, accessed 14.08.2018
#' @name w1500m
#' @docType data
NULL

#' Maiquetia Daily Rainfall
#'
#' Daily cumulated rainfall (in mm) at Maiquetia airport, Venezuela.
#' The observations cover the period from January 1961 to December 1999.
#' The original series had missing days in February 1996 (during which there were
#' 2 days with 1hr each of light rain) and January 1998 (no rain). These were replaced by zeros.
#'
#' @format a vector of size 14244 containing daily rainfall (in mm),
#' @source J.R. Cordova and M. González, accessed 25.11.2018 from <https://rss.onlinelibrary.wiley.com/hub/journal/14679876/series-c-datasets>
#' @name maiquetia
#' @references Coles, S. and L.R. Pericchi (2003). Anticipating Catastrophes through Extreme Value Modelling, \emph{Applied Statistics}, \bold{52}(4), 405-416.
#' @references  Coles, S., Pericchi L.R. and S. Sisson (2003). A fully probabilistic approach to extreme rainfall modeling, \emph{Journal of Hydrology}, \bold{273}, 35-50.
#' @docType data
#' @examples
#' \dontrun{
#' data(maiquetia, package = "mev")
#' day <- seq.Date(from = as.Date("1961-01-01"), to = as.Date("1999-12-31"), by = "day")
#' nzrain <- maiquetia[substr(day, 1, 4) < 1999 & maiquetia > 0]
#' fit.gpd(nzrain, threshold = 30, show = TRUE)
#'
#' }
NULL


#' River Nidd Flow
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
