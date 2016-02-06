\name{rates}
\alias{rates}
\title{Clustering Part of Conditional Intensity Function of the ETAS Model}
\description{
  A function to estimate the background seismicity rate and clustering
  (triggering) coefficient for a fitted ETAS model.
}
\usage{
  rates(fit, dimyx=NULL, method="zhuang", plot.it=TRUE)
}
\arguments{
  \item{fit}{A fitted ETAS model. An object of class \code{"etas"}.}
  \item{dimyx}{Dimensions of the rectangular discretization
           grid for the goegraphical study region.
           A numeric vector of length 2.}
  \item{method}{A character string specifying the method of smoothing.
           Currently methods "zhuang" and "spatstat" are implemented.}
  \item{plot.it}{Logical flag indicating whether to plot the rates or
           return them as pixel images.}
}
\value{
  If \code{plot.it=TRUE}, the function produces plots of the
  background seismicity rate and clustering coefficient.

  If \code{plot.it=FALSE}, it returns a list of length 3,
  with the total spatial intensity, background seismicity rate
  and clustering as components (objects of
  \code{im} class in the package \code{spatstat}).
}
\details{

  The argument \code{dimyx} determines the rectangular discretization
  grid dimensions. If it is given, then it must be a numeric vector
  of length 2 where the forst component \code{dimyx[1]} is the
  number of subdivisions in the y-direction (latitude) and the
  second component \code{dimyx[2]} is the number of subdivisions
  in the x-direction (longitude). The default (\code{NULL}) sets
  it to be \code{dimyx=c(128, 128)}.
}
\references{
  Zhuang, J., Ogata, Y. and Vere-Jones, D. (2005).
  Diagnostic analysis of space-time branching processes for earthquakes.
  Lecture Note in Statistics: \emph{Case Studies in Spatial Point Process Models}
  (Baddeley, A., Gregori, P., Mateu, J., Stoica, R. and Stoyan, D.),
  Springer-Verlag, New York, \bold{185}, 276--292.

  Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002).
  Stochastic declustering of space-time earthquake occurrences.
  \emph{Journal of the American Statistical Association},
  \bold{97}, 369--380.
}

\seealso{
  \code{\link{etas}}
}
\examples{

  # preparing the catalog
  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1996/01/01", study.end="2016/01/01",
     lat.range=c(25, 42), long.range=c(42, 63), mag.threshold=4.5)

  print(iran.cat)
  \dontrun{
  plot(iran.cat)}

  # initial parameters values
  param01 <- c(0.46, 0.23, 0.022, 2.8, 1.12, 0.012, 2.4, 0.35)

  # fitting the model and estimating the rates
  \dontrun{
  iran.fit <- etas(iran.cat, param0=param01, no.itr=5)
  rates(iran.fit, dimyx=c(200, 250))
  iran.rates <- rates(iran.fit, dimyx=c(200, 250), plot.it=FALSE)
  summary(iran.rates$background)}
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
