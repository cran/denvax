% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{morrison2010}
\alias{morrison2010}
\title{The serosurvey data in Morrison 2010}
\format{a data.frame (data.table, if installed) with 13 rows and 4 columns:
\describe{
  \item{Country}{character, common country name (all Peru for this data)}
  \item{Age}{integer, the age category}
  \item{Number}{integer, the number of samples}
  \item{Seropositive}{integer, the number of seropositive samples}
}}
\source{
\url{https://doi.org/10.1371/journal.pntd.0000670}
}
\usage{
morrison2010
}
\description{
From "Epidemiology of Dengue Virus in Iquitos, Peru 1999 to 2005: Interepidemic and Epidemic Patterns of Transmission",
combining information from Fig. 2 and Fig. 3.  The data from Fig. 3 were extracted using \url{https://automeris.io/WebPlotDigitizer/}
}
\examples{
require(denvax)
data(morrison2010)
with(morrison2010, plot(Age, Seropositive/Number*100, ylab="\% Seropositive", ylim=c(0,100)))
}
\keyword{datasets}
