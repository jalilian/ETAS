\name{ETAS-internal}
\title{Internal ETAS functions}
\alias{rates.inter}
\alias{lambdax}
\alias{etasfit}
\alias{decluster}
\alias{cxxfit}
\alias{cxxdeclust}
\alias{cxxrates}
\alias{cxxtimetrans}
\alias{cxxlambdtemp}
\alias{cxxlambspat}
\alias{timetransform}
\alias{lambdaspatial}
\alias{lambdatemporal}
\alias{longlat2xy}
\alias{xy2longlat}
\alias{print.etas}
\alias{plot.etas}
\alias{print.catalog}
\alias{plot.catalog}
\description{
  Internal ETAS functions.
}
\usage{
rates.inter(theta, object, bwd, lat.range=NULL, long.range=NULL,
            dimyx=NULL, plot.it=TRUE)
lambdax(rt, rx, ry, theta, revents)
etasfit(theta, revents, rpoly, tperiod, integ0, ihess, verbose,
    ndiv, eps, cxxcode, nthreads)
decluster(theta, rbwd, revents, rpoly, tperiod, ndiv, cxxcode)
cxxfit(tht, revents, rpoly, tperiod, rinteg0, ihess, ndiv, eps, verbose, nthreads)
cxxdeclust(param, revents, rpoly, bwd, tperiod, ndiv)
cxxrates(param, revents, bwd, tperiod, gx, gy)
cxxtimetrans(theta, revents, rpoly, tperiod, integ0, ndiv)
cxxlambdtemp(tg, theta, revents, rpoly, tperiod, integ0, ndiv)
cxxlambspat(xg, yg, theta, revents, rpoly, tperiod, bwd)
timetransform(fit)
lambdaspatial(x, y, fit)
lambdatemporal(t, fit)
longlat2xy(long, lat, region.poly, dist.unit="degree")
xy2longlat(x, y, region.poly, dist.unit="degree")
\method{print}{etas}(x,\dots)
\method{plot}{etas}(x, which="est", dimyx=NULL, \dots)
\method{print}{catalog}(x,\dots)
\method{plot}{catalog}(x,\dots)

}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
