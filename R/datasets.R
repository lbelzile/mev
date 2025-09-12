#' French wind data
#'
#' Daily mean wind speed (in km/h) at four stations in the south of France, namely Cap Cepet (\code{S1}), Lyon St-Exupery (\code{S2}), Marseille Marignane (\code{S3}) and Montelimar (\code{S4}).
#' The data includes observations from January 1976 until April 2023; days containing missing values are omitted.
#' @source European Climate Assessment and Dataset project \url{https://www.ecad.eu/}
#' @references
#' Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the
#' European Climate Assessment. Int. J. of Climatol., 22, 1441-1453.
#' @name frwind
#' @docType data
#' @format A data frame with 17209 observations and 8 variables:
#' \describe{
#' \item{\code{date}}{date of measurement}
#' \item{\code{S1}}{wind speed (in km/h) at Cap Cepet}
#' \item{\code{S2}}{wind speed (in km/h) at Lyon Saint-Exupery}
#' \item{\code{S3}}{wind speed (in km/h) at Marseille Marignane}
#' \item{\code{S4}}{wind speed (in km/h) at Montelimar}
#' \item{\code{H2}}{humidity (in percentage) at Lyon Saint-Exupery}
#' \item{\code{T2}}{mean temperature (in degree Celcius) at Lyon Saint-Exupery}
#' }
#' The \code{metadata} attribute includes latitude and longitude (in degrees, minutes, seconds), altitude (in m), station name and station id.
#' @examples
#' data(frwind, package = "mev")
#' head(frwind)
#' attr(frwind, which = "metadata")
NULL

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
#' from 1887 until 2019. Only the yearly maximum is available for 1922
#' and the six largest observations for 1936.
#'
#' @format a data frame with 133 rows and 11 columns containing the year of the measurement (first column)
#' and ordered 10-largest yearly observations, reported in decreasing order from largest (\code{r1}) to smallest (\code{r10}).
#'
#' @note Smith (1986) notes that the annual maxima seems to fluctuate around a constant sea level
#' up to 1930 or so, after which there is potential linear trend. Records of threshold exceedances above
#' 80 cm (reported on the website) indicate that observations are temporally clustered.
#'
#' The observations from 1931 until 1981 can be found in
#' Table 1 in Smith (1986), who reported data from Pirazzoli (1982).
#' The values from 1983 until 2019 were extracted by Anthony Davison from the City
#' of Venice website (accessed in May 2020) and are licensed under the CC BY-NC-SA 3.0 license.
#' The Venice City website indicates
#' that later measurements were recorded by an instrument located in Punta Salute.
#'
#'
#' @references Smith, R. L. (1986) Extreme value theory based on the \emph{r}
#' largest annual events. \emph{Journal of Hydrology} \bold{86}, 27–43.
#' @references Pirazzoli, P., 1982. Maree estreme a Venezia (periodo 1872-1981). \emph{Acqua Aria} \bold{10}, 1023-1039.
#' @references  Coles, S. G. (2001) \emph{An Introduction to Statistical Modelling of Extreme Values}. London: Springer.
#' @source City of Venice, Historical archive <https://www.comune.venezia.it/node/6214>. Last accessed November 5th, 2020.
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
#' @references Davison, A.C. and R.L. Smith (1990). Models for Exceedances over High Thresholds (with discussion), \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{52}(3), 393--442.
#' @format a vector of size 154
#' @name nidd
#' @seealso \code{nidd.thresh} from the \code{evir} package
#' @docType data
NULL

#' Magnetic storms
#'
#' Absolute magnitude of 373 geomagnetic storms lasting more than 48h with absolute magnitude (dst) larger than 100 in 1957-2014.
#'
#' @source Aki Vehtari
#' @references World Data Center for Geomagnetism, Kyoto, M. Nose, T. Iyemori, M. Sugiura, T. Kamei (2015), \emph{Geomagnetic Dst index}, <doi:10.17593/14515-74000>.
#' @docType data
#' @note For a detailed article presenting the derivation of the Dst index, see \code{http://wdc.kugi.kyoto-u.ac.jp/dstdir/dst2/onDstindex.html}
#' @format a vector of size 373
#' @name geomagnetic
NULL

#' Abisko rainfall
#'
#' Daily non-zero rainfall measurements in Abisko (Sweden) from January 1913 until December 2014.
#' @param date \code{Date} of the measurement
#' @param precip rainfall amount (in mm)
#' @format a data frame with 15132 rows and two variables
#' @name abisko
#' @source Abisko Scientific Research Station
#' @references A. Kiriliouk, H. Rootzen, J. Segers and J.L. Wadsworth (2019), \emph{Peaks over thresholds modeling with multivariate generalized Pareto distributions},  Technometrics, \bold{61}(1), 123--135, <doi:10.1080/00401706.2018.1462738>
NULL

#' Nutrient data
#'
#' Interview component of survey 'What we eat in
#' America'. These are extracted from the 2015–2016 National Health and Nutrition Examination Survey (NHANES, \url{https://wwwn.cdc.gov/nchs/nhanes/Default.aspx}) report and consist of the total nutrients for all food and beverage intake ingested over a 24 hours period.
#' @source National Center for Health Statistics, now available from the Wayback Machine via \url{https://web.archive.org/web/20201029113801/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DR1TOT_I.XPT}
#'
#' @details Note that the sample design oversampled specific population targets and that only respondants are provided. The website contains more information about sampling weights. There are multiple missing records.
#' @format A data frame with 9544 rows and 38 variables:
#' \describe{
#'   \item{\code{prot}}{proteins (in grams)}
#'   \item{\code{carb}}{carbonhydrate (in gram)}
#'   \item{\code{sugr}}{total sugars (in gram)}
#'   \item{\code{fibe}}{dietary fibers (in grams)}
#'   \item{\code{tfat}}{total fat (in grams)}
#'   \item{\code{sfat}}{saturated fat (in grams)}
#'   \item{\code{mfat}}{monounsaturated fat (in grams)}
#'   \item{\code{pfat}}{polyunsaturated fat (in grams)}
#'   \item{\code{chol}}{cholesterol (in milligrams)}
#'   \item{\code{atoc}}{vitamin E as alpha-tocopherol (in milligrams)}
#'   \item{\code{ret}}{retinol (in micrograms)}
#'   \item{\code{vara}}{Vitamin A as retinol activity equivalents (in micrograms).}
#'   \item{\code{acar}}{alpha-carotene (in micrograms)}
#'   \item{\code{bcar}}{beta-carotene (in micrograms)}
#'   \item{\code{cryp}}{beta-cryptoxanthin (in micrograms)}
#'   \item{\code{lyco}}{lycopene (in micrograms)}
#'   \item{\code{lz}}{lutein and zeaxanthin (in micrograms).}
#'   \item{\code{vb1}}{thiamin (vitamin B1, in milligrams)}
#'   \item{\code{vb2}}{riboflavin (vitamin B2, in milligrams)}
#'   \item{\code{niac}}{niacin (in milligrams)}
#'   \item{\code{vb6}}{vitamin B5 (in milligrams)}
#'   \item{\code{fola}}{total folate (in micrograms)}
#'   \item{\code{fa}}{folic acid (in micrograms)}
#'   \item{\code{ff}}{food folate (in micrograms)}
#'   \item{\code{chl}}{total choline (in milligrams)}
#'   \item{\code{vb12}}{vitamin B12 (in micrograms)}
#'   \item{\code{vc}}{vitamin C (in milligrams)}
#'   \item{\code{vd}}{vitamin D (comprising D2 and D3, in micrograms)}
#'   \item{\code{vk}}{vitamin K (in micrograms)}
#'   \item{\code{calc}}{calcium (in milligrams)}
#'   \item{\code{phos}}{phosphorus (in milligrams)}
#'   \item{\code{magn}}{magnesium (in milligrams)}
#'   \item{\code{iron}}{iron (in milligrams)}
#'   \item{\code{zinc}}{zinc (in milligrams)}
#'   \item{\code{copp}}{copper (in milligrams)}
#'   \item{\code{sodi}}{sodium (in milligrams)}
#'   \item{\code{pota}}{potassium (in milligrams)}
#'   \item{\code{sele}}{selenium (in micrograms)}
#'}
#' @name nutrients
#' @note These data are subject to a data user agreement, available at \url{https://www.cdc.gov/nchs/policy/data-user-agreement.html}
"nutrients"


#' Deaths from pandemics
#'
#' The data base contains estimated records of the number of deaths from pandemics.
#' @format A data frame with 72 rows and 8 variables:
#' \describe{
#'   \item{\code{event}}{name of the event}
#'   \item{\code{startyear}}{start year of the event}
#'   \item{\code{endyear}}{end year of the event}
#'   \item{\code{lower}}{lower bound on estimated deaths (in thousands)}
#'   \item{\code{average}}{average estimated deaths (in thousands)}
#'   \item{\code{upper}}{upper bound on estimated deaths (in thousands)}
#'   \item{\code{saverage}}{scaled average of estimated deaths (in thousands)}
#'   \item{\code{population}}{estimated population at risk (in thousands)}
#' }
#' @name pandemics
#' @source Cirillo, P. and N.N. Taleb (2020). \emph{Tail risk of contagious diseases}. Nat. Phys. \bold{16}, 606–613 (2020). <doi:10.1038/s41567-020-0921-x>
"pandemics"


#' Leeds air pollution
#'
#' Daily maximum data (hourly for PM10) on air pollution for the Leeds Centre station in Yorkshire and Humberside station. The data goes from January 1st, 1993, until December 31st, 2024. Data show seasonality and there are some outliers. From December 2nd, 2008 onwards, particulate matters (PM10 and PM2.5) are measured using  a tapered element oscillating microbalance (TEOM) and Filter Dynamics Measurement System (FDMS). The data for PM2.5 is missing before the change of instrumentation. A total of 231 daily measurements with only missing values were removed during preprocessing.
#'
#' @format A data frame with 11455 rows and 8 variables:
#' \describe{
#'   \item{\code{date}}{[character] a date with format yyy-mm-dd}
#'   \item{\code{O3}}{[integer] ozone (in nanograms per cubic meter)}
#'   \item{\code{NO}}{[integer] nitrogen oxyde (in nanograms per cubic meter)}
#'   \item{\code{CO}}{[double] carbon monoxyde (in micrograms per cubic meter)}
#'   \item{\code{NO2}}{nitrogen dioxyde  (in nanograms per cubic meter)}
#'   \item{\code{SO2}}{sulphur dioxide  (in nanograms per cubic meter)}
#'   \item{\code{PM10}}{[integer] particulate matter 10, (in nanograms per cubic meter)}
#'   \item{\code{PM2.5}}{[integer] particulate matter 2.5, (in nanograms per cubic meter)}
#'}
#'
#' @source Crown 2025 copyright Defra via \code{uk-air.defra.gov.uk}, licenced under the Open Government Licence (OGL).
"leedspollution"
