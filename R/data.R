#' The serosurvey data in Morrison 2010
#'
#' From "Epidemiology of Dengue Virus in Iquitos, Peru 1999 to 2005: Interepidemic and Epidemic Patterns of Transmission",
#' combining information from Fig. 2 and Fig. 3.  The data from Fig. 3 were extracted using \url{https://automeris.io/WebPlotDigitizer/}
#'
#' @format a data.frame (data.table, if installed) with 13 rows and 4 columns:
#' \describe{
#'   \item{Country}{character, common country name (all Peru for this data)}
#'   \item{Age}{integer, the age category}
#'   \item{Number}{integer, the number of samples}
#'   \item{Seropositive}{integer, the number of seropositive samples}
#' }
#'
#' @source \url{https://doi.org/10.1371/journal.pntd.0000670}
#'
#' @examples
#' require(denvax)
#' data(morrison2010)
#' with(morrison2010, plot(Age, Seropositive/Number*100, ylab="% Seropositive", ylim=c(0,100)))
"morrison2010"

#' The serosurvey data in L'Azou 2016
#'
#' From "Symptomatic Dengue in Children in 10 Asian and Latin American Countries", Table 4.
#'
#' @format a data.frame (data.table, if installed) with 20 rows and 4 columns:
#' \describe{
#'   \item{Country}{character, common country name (all Peru for this data)}
#'   \item{Age}{character, the bounding ages for the sample; format: lower age '-' upper age}
#'   \item{Number}{integer, the number of samples}
#'   \item{Seropositive}{integer, the number of seropositive samples}
#' }
#'
#' @source \url{https://doi.org/10.1056/NEJMoa1503877}
#' @examples
#' require(denvax); require(ggplot2)
#' data(lazou2016)
#' ggplot(lazou2016) + aes(Age, Seropositive/Number*100, color = Country) +
#'   geom_point() + labs(y="Seropositive %", x="Age Group") + lims(y=c(0,100)) +
#'   theme_minimal()
"lazou2016"
