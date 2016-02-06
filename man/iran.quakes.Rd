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
\references{
  Jalilian, A. (2016). Modeling Earthquakes of Iran.
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
