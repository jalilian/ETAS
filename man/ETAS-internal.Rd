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
\alias{cxxSmooth}
\alias{cxxstpoisstest}
\alias{cxxstpoisstestMP}
\alias{cdeclust}
\alias{cfit}
\alias{clambdax}
\alias{timetransform}
\alias{lambdaspatial}
\alias{lambdatemporal}
\alias{longlat2xy}
\alias{xy2longlat}
\alias{print.etas}
\alias{plot.etas}
\alias{print.catalog}
\alias{plot.catalog}
\alias{poiss.test}
\alias{Smoothcatalog}
\alias{allanfactor}
\alias{morisitaindex}
\alias{decimalplaces}
\alias{roundoffErr}
\alias{ttrend}
\alias{simetas}
\description{
  Internal ETAS functions.
}
\usage{
rates.inter(theta, object, bwd, lat.range=NULL, long.range=NULL,
            mver=1, dimyx=NULL, plot.it=TRUE)
lambdax(rt, rx, ry, theta, revents)
etasfit(theta, revents, rpoly, tperiod, integ0, ihess, verbose,
    ndiv, eps, cxxcode, nthreads, mver)
decluster(theta, rbwd, revents, rpoly, tperiod, ndiv, cxxcode, mver)
cxxfit(tht, revents, rpoly, tperiod, rinteg0, ihess, ndiv, eps, verbose, nthreads, mver)
cxxdeclust(param, revents, rpoly, bwd, tperiod, ndiv, mver)
cxxrates(param, revents, bwd, tperiod, gx, gy, mver)
cxxtimetrans(theta, revents, rpoly, tperiod, integ0, ndiv, mver)
cxxlambdtemp(tg, theta, revents, rpoly, tperiod, integ0, ndiv, mver)
cxxlambspat(xg, yg, theta, revents, rpoly, tperiod, bwd, mver)
cxxSmooth(x, y, bwd, gx, gy, expand)
cxxstpoisstest(xrank, yrank, M)
cxxstpoisstestMP(xrank, yrank, M, nthreads)
cdeclust(theta, rbwd,revents, rpoly, tperiod)
cfit(theta, rdata, ihess, rverbose)
clambdax(rt, rx, ry, theta, revents)
timetransform(fit)
lambdaspatial(x, y, fit)
lambdatemporal(t, fit)
longlat2xy(long, lat, region.poly, dist.unit="degree")
xy2longlat(x, y, region.poly, dist.unit="degree")
poiss.test(object, which="joint", r=NULL, lambda=NULL, bwd=NULL, dimyx=NULL, 
  nsim=299, n.perm=1000, verbose=TRUE, nthreads=1)
Smoothcatalog(object, type="spatial", bwd=NULL, bwm=NULL, nnp=NULL, 
  dimyx=NULL, convert=FALSE)
allanfactor(object, K=100, nsim=500, cat.name=NULL)
morisitaindex(object, K=11, bwd=NULL, dimyx=NULL, cat.name=NULL)
decimalplaces(x) 
roundoffErr(x)
ttrend(obj, tmax=NULL, tmin=NULL)
simetas(param, bkgd, sim.start, sim.end=NULL, sim.length=NULL, 
        lat.range=NULL, long.range=NULL, region.poly=NULL, 
        mag.threshold=NULL, flatmap=TRUE, dist.unit="degree", 
        cat0=NULL, roundoff=FALSE, tz="GMT")
\method{print}{etas}(x,\dots)
\method{plot}{etas}(x, which="est", dimyx=NULL, \dots)
\method{print}{catalog}(x,\dots)
\method{plot}{catalog}(x,\dots)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
