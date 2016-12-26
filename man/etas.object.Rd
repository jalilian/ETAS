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

  # preparing the catalog
  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1985/01/01", study.end="2016/01/01",
     lat.range=c(25, 42), long.range=c(42, 63), mag.threshold=4.5)

  \dontrun{
  plot(ir.cat)
  }

  # initial parameters values
  param01 <- c(0.4339678,
   0.1988628,
   0.0345206,
   1.6290137,
   1.1286776,
   0.0072539,
   2.1705884,
   0.5706402)

  # fitting the model
  \dontrun{
  res <- etas(ir.cat, param0=param01)
  }
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{attribute}

