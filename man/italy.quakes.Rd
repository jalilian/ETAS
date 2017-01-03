\name{italy.quakes}
\alias{italy.quakes}
\docType{data}
\title{An Italian Earthquake Catalog}
\description{
  A data frame with 2158 rows and 6 columns giving occurrence date,
  time, longitude, latitude, magnitude and depth of earthquakes
  occurred since 2005-04-16 till 2013-11-01 in Italy and its vicinity
  (6.15-19E and 35-48N) with magnitude greater than or equal to 3.
}
\format{
  An object of class \code{"data.frame"} containing the following columns:
  \itemize{
  \item \code{date}
  Occurrence date in the format "yyyy/mm/dd"
  \item \code{time}
  Occurrence time in decimal days after the first earthquake
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
\usage{data(italy.quakes)}
\source{
  Data are retrieved from the Italian Seismological Instrumental and
  Parametric Data Base (ISIDE) available at \url{http://iside.rm.ingv.it}
  by Marcello Chiodi and Giada Adelfio in the \pkg{etasFLP} package.
}

\references{
  Marcello Chiodi and Giada Adelfio (2015). etasFLP: Mixed FLP and ML
  Estimation of ETAS Space-Time Point Processes. R package version
  1.3.0. \url{https://CRAN.R-project.org/package=etasFLP}.

}
\examples{
  # creat an earthquake catalog
  italy.cat <- catalog(italy.quakes, dist.unit="km")

  \dontrun{
  plot(italy.cat)}

  # set initial parameter values
  mu <- 1
  k0 <- 0.005
  c <- 0.005
  alpha <- 1.05
  p <- 1.01
  D <- 1.1
  q <- 1.52
  gamma <- 0.6
  # reparametrization: transform k0 to A
  A <- pi * k0 / ((p - 1) * c^(p - 1) * (q - 1) * D^(q - 1))
  param0 <- c(mu, A, c, alpha, p, D, q, gamma)

  \dontrun{
  italy.fit <- etas(italy.cat, param0)}
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
