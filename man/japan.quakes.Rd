\name{japan.quakes}
\alias{japan.quakes}
\alias{jap.quakes}
\docType{data}
\title{A Japanese Earthquake Catalog}
\description{
  A data frame with 13724 rows and 6 columns giving occurrence date,
  time, longitude, latitude, magnitude and depth of shallow earthquakes
  (depth < 100 km) occurred since 1926-01-08 till
  2007-12-29 in Japan and its vicinity (128-145E and 27-45N).
  Only earthquakes with magnitude greater than or equal to 4.5
  are included.
}
\format{
  An object of class \code{"data.frame"} containing the following columns:
  \itemize{
  \item \code{date}
  Occurrence date in the format "yyyy/mm/dd"
  \item \code{time}
  Occurrence time in the format "hh:mm:ss.ss"
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
\usage{data(japan.quakes)}
\source{
  Data are retrieved from the Japan Meteorological Agency (JMA) by
  Jiancang Zhuang and acoppany the \bold{Fortran} code at
  \url{http://bemlar.ism.ac.jp/zhuang/software.html}.
}

\references{
  Zhuang J, Ogata Y, Vere-Jones D (2002). Stochastic Declustering of
  Space-Time Earthquake Occurrences.
  \emph{Journal of the American Statistical Association},
  \bold{97}(458), 369--380.
  \href{http://dx.doi.org/10.1198/016214502760046925}{doi:10.1198/016214502760046925}.

  Zhuang J, Ogata Y, Vere-Jones D (2006). Diagnostic Analysis of Space-Time
  Branching Processes for Earthquakes. In \emph{Case Studies in Spatial Point
  Process Modeling}, pp. 275--292. Springer Nature.
  \href{http://dx.doi.org/10.1007/0-387-31144-0_15}{doi:10.1007/0-387-31144-0_15}.
}
\examples{
  # set the target polygonal study region
  jpoly <- list(long=c(134.0, 137.9, 143.1, 144.9, 147.8,
      137.8, 137.4, 135.1, 130.6), lat=c(31.9, 33.0, 33.2,
      35.2, 41.3, 44.2, 40.2, 38.0, 35.4))

  # creat an earthquake catalog
  japan.cat <- catalog(japan.quakes, study.start="1953-05-26",
      study.end="1990-01-08", region.poly=jpoly, mag.threshold=4.5)

  \dontrun{
  plot(japan.cat)}

  param0 <- c(0.592844590, 0.204288231, 0.022692883, 1.495169224,
  1.109752319, 0.001175925, 1.860044210, 1.041549634)

  \dontrun{
  nthreads <- parallel::detectCores()
  japan.fit <- etas(japan.cat, param0, nthreads=nthreads)}

}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
