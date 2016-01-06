\name{jap.quakes}
\alias{jap.quakes}
\docType{data}
\title{Earthquakes of Japan}
\description{
  A data frame with 9189 rows and 6 columns giving occurrence date,
  time (in decimal days after the first earquake),
  longitude, latitude, magnitude and depth of shalow earthquakes
  (depth < 100 km) occurred since 1966-01-01 till
  2016-01-01 in Japan and its vicinity (128-145E and 27-45N).
  Only earthquakes with magnitude greater than or equal to 4.5
  are included.
}
\format{
  An object of class \code{"data.frame"} containing the following columns:
  \itemize{
  \item \code{date}
  Occurrence date in the format "yyyy/mm/dd"
  \item \code{time}
  Occurrence time in decimal days after the first earquake
  \item \code{long}
  Latitude of epicenter in decimal degrees
  \item \code{lat}
  Latitude of epicenter in decimal degrees
  \item \code{mag}
  Magnitude of each earthquake by JMA (Japan Meteorological Agency)
  \item \code{depth}
  Depth of each earthquake
  }
}
\usage{data(jap.quakes)}
\source{
  Data retrieved from the ISC Bulletin event catalogue
}

\references{
  International Seismological Centre, On-line Bulletin,
  \url{http://www.isc.ac.uk}, Internatl. Seis. Cent.,
  Thatcham, United Kingdom, 2013.
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
