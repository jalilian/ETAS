\name{catalog}
\alias{catalog}
\title{Create an Earthquake Catalog}
\description{
  Creates an object of class \code{"catalog"} representing
  an earthquake catalog dataset. An  earthquake catalog is a
  chronologically ordered list of time, epicenter and magnitude
  of all recorded earthquakes  in geographical region during
  a specific time period.
}
\usage{
  catalog(data, time.begin=NULL, study.start=NULL,
          study.end=NULL, study.length=NULL,
          lat.range=NULL, long.range=NULL,
          region.poly=NULL, mag.threshold=NULL,
          flatmap=TRUE, tz="GMT")
}
\arguments{
  \item{data}{A \code{data.frame} contaning date, time,
  latitude, longitude and magnitude of earthquakes.}
  \item{time.begin}{The beginning of time span of the catalog.
  A character string or an object that can be converted to
  date-time (calendar dates plus time to the nearest second) by
  \code{as.POSIXlt}. The default \code{NULL} sets it
  to the date-time of the first event.}
  \item{study.start}{The start of the study period.
  A character string or an object that can be converted to
  date-time by \code{as.POSIXlt}. If not specified (\code{NULL}),
  then \code{time.begin} is used.}
  \item{study.end}{The end of the study period.
  A character string or an object that can be converted to
  date-time by \code{as.POSIXlt}. The default \code{NULL} sets it
  to the date-time of the last event.}
  \item{study.length}{A single numeric value specifying the length
  of the study period in decimal days. Incompatible with
  \code{study.end}: either \code{study.end} or \code{study.length}
  can be specified, but not both.}
  \item{lat.range}{The latitude range of a rectangular study region.
  A numeric vector of size 2 giving (latmin, latmax). By default
  (\code{NULL}) the range of the latitudes of events is used.}
  \item{long.range}{The longitude range of a rectangular study region.
  A numeric vector of size 2 giving (longmin, longmax). By default
  (\code{NULL}) the range of the longitudes of events is used.}
  \item{region.poly}{Polygonal boundary of a non-rectangular
  study region. A list with components \bold{lat} and \bold{long}
  of equal length specifying the coordinates of the vertices of
  a polygonal study region. The vertices must be listed in
  \bold{anticlockwise} order.}
  \item{mag.threshold}{The magnitude threshold of the catalog.
  A positive numeric value. The default (\code{NULL}) sets it to
  the minimum magnitude of all events.}
  \item{flatmap}{Logical flag indicating whether to adjust the
  longitudes of events and the study region in order to have a
  flat map coordinates.}
  \item{tz}{A character string specifying the time zone to be used
  for the date-time conversion in \code{as.POSIXlt}.
  The default \code{"GMT"} is the UTC (Universal Time, Coordinated).}
}
\value{
  An object of class \code{"catalog"} containing an earthquake
  catalog dataset.
}
\details{
  The \code{data} is required to have at least 5 columns with names
  \code{date}, \code{time}, \code{lat}, \code{long} and \code{mag}
  containing, respectively, the date, time, latitude, longitude
  and magnitude of each event in the catalog.

  The geographical study region can be rectangular or polygonal:
 \itemize{
    \item
    \bold{rectangular study region} can be specified by \code{lat.range}
    and \code{long.range} which must be numeric vectors of length 2.
    \item
    \bold{polygonal study region} can be specified by \code{region.poly}
    which containst coordinates of the vertices of the polygon. It must
    be eighter a \code{list} with components \bold{lat} and \bold{long}
    of equal length or a \code{data.frame} with columns \bold{lat}
    and \bold{long}. The vertices must be listed in
    \emph{anticlockwise} order and no vertex should be repeated
    (i.e. do not repeat the first vertex).
  }

  The function \code{\link{inside.owin}} in the \code{spatstat}
  is used to indicate whether events lie inside the study region.
  Only events inside the study region and the study period
  (\code{study.start}, \code{study.end}) are considered as
  \emph{target} events. Other events are assumed to be
  \emph{complementary} events.

  If the events in \code{data} are not chronologically sorted,
  then a warning will be produced and the events will be sorted
  in ascending order with respect to time of occurrence.

  If \code{flatmap=TRUE}, longitude-latitude coordinates convert to
  flat map coordinates. The algorithm does not change the latitude
  of points but adjust the longitude of points with respect to the
  centroid of the target geographical region in order to
  get coordinates on a flat surface.
}
\seealso{
  \code{etas}.
}
\references{
  Zhuang, J. (2012). Long-term earthquake forecasts based on the
  epidemic-type aftershock sequence (ETAS) model for short-term
  clustering. \emph{Research in Geophysics}, \bold{2}(1), 8.
}

\examples{
  data(iran.quakes)
  summary(iran.quakes)

  # creating a catalog with rectangular study region
  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1985/01/01", study.end="2016/01/01",
     lat.range=c(25, 42), long.range=c(42, 63), mag.threshold=4.5)

  print(iran.cat)
  \dontrun{
  plot(iran.cat)
  }

  # equivalently, specifying the length of the study period
  iran.cat2 <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1985/01/01", study.length=11322,
     lat.range=c(25, 42), long.range=c(42, 63), mag.threshold=4.5)

  print(iran.cat2)



  data(jap.quakes)
  summary(jap.quakes)

  # specifying a polygonal geographical region
  jpoly <- list(long=c(134.0, 137.9, 143.1, 144.9, 147.8,
      137.8, 137.4, 135.1, 130.6), lat=c(31.9, 33.0, 33.2,
      35.2, 41.3, 44.2, 40.2, 38.0, 35.4))
  # creating a catalog with polygonal study region
  jap.cat <- catalog(jap.quakes, time.begin="1966-01-01",
      study.start="1970-01-01", study.end="2010-01-01",
      region.poly=jpoly, mag.threshold=4.5)

  print(jap.cat)
  \dontrun{
  plot(jap.cat)
  }
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
