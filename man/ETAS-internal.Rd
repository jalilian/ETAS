\name{ETAS-internal}
\title{Internal ETAS functions}
\alias{rates.inter}
\alias{lambdax}
\alias{etasfit}
\alias{decluster}
\alias{print.etas}
\alias{plot.etas}
\alias{print.catalog}
\alias{plot.catalog}
\description{
  Internal ETAS functions.
}
\usage{
rates.inter(theta, object, bwd, dimyx=NULL,
            method="zhuang", plot.it=TRUE)
lambdax(rt, rx, ry, theta, revents)
etasfit(theta, revents, rpoly, tperiod, integ0, verbose)
decluster(theta, rbwd, revents, rpoly, tperiod)
\method{print}{etas}(x,\dots)
\method{plot}{etas}(x, which="est", dimyx=NULL, \dots)
\method{print}{catalog}(x,\dots)
\method{plot}{catalog}(x,\dots)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
