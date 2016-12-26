\name{iran.quakes}
\alias{iran.quakes}
\docType{data}
\title{An Iranian Earthquake Catalog}
\description{
  A data frame with 5970 rows and 5 columns giving occurrence date, time,
  longitude, latitude and magnitude of shallow earthquakes (depth < 100 km)
  occurred since 1973-01-01 till 2016-01-01 in Iran and its vicinity
  (40-65E and 22-42N). Only earthquakes with magnitude greater than
  or equal to 4.5 are included.
}
\format{
  An object of class \code{"data.frame"} containing the following columns:
  \itemize{
  \item \code{date}
  Occurrence date in the format "yyyy-mm-dd"
  \item \code{time}
  Occurrence time (UTC) in the format "hh:mm:ss"
  \item \code{long}
  Latitude of epicenter in decimal degrees
  \item \code{lat}
  Latitude of epicenter in decimal degrees
  \item \code{mag}
  Magnitude in body-wave magnitude scale (mb)
  }
}
\usage{data(iran.quakes)}
\source{
  The ANSS Comprehensive Catalog (ComCat):
  \url{http://earthquake.usgs.gov/earthquakes/search/}
}
%\references{
%  Jalilian, A. (2017). Modeling Earthquakes of Iran by ETAS model. \emph{submitted}.
%}
\examples{
  summary(iran.quakes)

  gregion <- list(lat = c(26, 25, 29, 38, 35), long = c(52, 59, 58, 45, 43))
  # creat an earthquake catalog
  iran.cat <- catalog(iran.quakes, study.start = "1991/01/01",
     study.end = "2011/01/01", region.poly = gregion, mag.threshold = 4.5)


  \dontrun{
  plot(iran.cat)
  iran.fit <- etas(iran.cat)}

  zagros <- list(lat = c(27, 26, 29, 29, 35, 33),
     long = c(52, 58, 58, 54, 48, 46))
  iran.cat <- catalog(iran.quakes, study.start = "1991/01/01",
     study.end = "2011/01/01", region.poly = zagros, mag.threshold = 4)
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
