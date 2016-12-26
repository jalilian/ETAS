\name{lambda}
\alias{lambda}
\title{Clustering Part of Conditional Intensity Function of the ETAS Model}
\description{
  A function to compute the clustering part of the conditional
  intensity function of the ETAS model at specified time and location.
}
\usage{
  lambda(t, x, y, param, object)
}
\arguments{
  \item{t}{A numeric value. The time that the conditional intensity is
            to be computed at.}
  \item{x}{A numeric value. The x-coordinate of the location that the
           conditional intensity is to be computed at.}
  \item{y}{A numeric value. The y-coordinate of the location that the
           conditional intensity is to be computed at.}
  \item{param}{Vector of model parameters.}
  \item{object}{An object of class \code{"catalog"} containing an
  earthquake catalog dataset.}
}
\value{
  A numeric value.
}
\details{
  For a given \eqn{t}{t}, \eqn{x}{x} and \eqn{y}{y}, this function
  computes
    \deqn{
      \sum_{t_i < t} k(m_i)g(t - t_i)f(x - x_i, y - y_i|m_i).}{sum[t[i] < t] k(m[i]) g(t - t[i]) f(x - x[i], y - y[i]|m[i]).}
}
\references{
  Zhuang J, Ogata Y, Vere-Jones D (2002).
  Stochastic Declustering of Space-Time Earthquake Occurrences.
  \emph{Journal of the American Statistical Association},
  \bold{97}(458), 369--380.
  \href{http://dx.doi.org/10.1198/016214502760046925}{doi:10.1198/016214502760046925}.

  Zhuang J, Ogata Y, Vere-Jones D (2006).
  Diagnostic Analysis of Space-Time Branching Processes for Earthquakes.
  In \emph{Case Studies in Spatial Point Process Modeling},
  pp. 275--292. Springer Nature.
  \href{http://dx.doi.org/10.1007/0-387-31144-0_15}{doi:10.1007/0-387-31144-0_15}.
}

\seealso{
  \code{\link{etas}}
  \code{\link{catalog}}
}
\examples{

  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1996/01/01", study.end="2016/01/01",
     lat.range=c(25, 42), long.range=c(42, 63), mag.threshold=4.5)

  param <- c(0.46, 0.23, 0.022, 2.8, 1.12, 0.012, 2.4, 0.35)

  \dontrun{
  lambda(15706, 40.12, 34.5, param, iran.cat)}

}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
