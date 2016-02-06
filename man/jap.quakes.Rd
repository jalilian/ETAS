\name{jap.quakes}
\alias{jap.quakes}
\docType{data}
\title{Earthquakes of Japan}
\description{
  A data frame with 9189 rows and 6 columns giving occurrence date,
  time (in decimal days after the first earquake),
  longitude, latitude, magnitude and depth of shallow earthquakes
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

   Di Giacomo, D., Harris, J., Villasenor, A., Storchak, D. A.,
   Engdahl, E. R., Lee, W. H., & Team, D. E. (2015).
   ISC-GEM: Global Instrumental Earthquake Catalogue (1900-2009),
   I. Data collection from early instrumental seismological bulletins.
   \emph{Physics of the Earth and Planetary Interiors}, \bold{239}, 14-24.
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
