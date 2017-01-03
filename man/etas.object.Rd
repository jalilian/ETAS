\name{etas.object}
\alias{etas.object}
\alias{methods.etas}
\title{Class of Fitted ETAS Models}
\description{
  A class \code{etas} to represent a fitted ETAS model.
  The output of \code{\link{etas}}.
}
\details{
  An object of class \code{etas} represents an ETAS model
   that has been fitted to a spatio-temporal point pattern (catalog)
  of earthquakes. It is the output of the model fitter, \code{\link{etas}}.

  The class \code{etas} has methods for the following
  standard generic functions:

  \tabular{lll}{
    generic \tab method \tab description \cr
    \code{print} \tab \code{\link{print.etas}}
    \tab print details \cr
  }
}

\seealso{
  \code{\link{etas}},
}
\examples{
  # fitting the ETAS model to an Iranian catalog

  data(iran.quakes)
  summary(iran.quakes)

    # fitting the ETAS model to an Iranian catalog
  # preparing the catalog
  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1986/01/01", study.end="2016/01/01",
     lat.range=c(26, 40), long.range=c(44, 63), mag.threshold=5)
  print(iran.cat)
  \dontrun{
  plot(iran.cat)}

  # setting initial parameter values
  param0 <- c(0.46, 0.23, 0.022, 2.8, 1.12, 0.012, 2.4, 0.35)

  # fitting the model
  \dontrun{
  iran.fit <- etas(iran.cat, param0=param0)}

}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{attribute}

