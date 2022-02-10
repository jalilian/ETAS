##############################
# temporal trends: cumulative sums of events

ttrend <- function(obj, tmax=NULL, tmin=NULL)
{
  tt <- obj$revents[obj$revents[, 5] == 1, 1]
  if (is.null(tmin))
    tmin <- obj$rtperiod[1]
  if (is.null(tmax))
    tmax <- obj$rtperiod[2]
  tt <- tt[tt >= tmin]
  tvs <- seq(tmin, tmax, length.out=500)
  lamb <- length(tt)/(tmax - tmin)
  ffun <- function(tmp)
  {
    colSums(outer(tmp, tvs, "<=")) / sqrt(tvs - tmin) - lamb * sqrt(tvs - tmin)
  }
  vv <- ffun(tt)
  if (tt[1] == 0)
    vv[1] <- NA
  
  simfun <- function(i) 
  {
    N <- stats::rpois(1, lamb * (tmax - tmin))
    ttsim <- runif(N, tmin, tmax)
    ffun(ttsim)
  }
  vvsim <- matrix(unlist(lapply(1:99, simfun)), ncol=99)
  lo <- apply(vvsim, 1, min)
  up <- apply(vvsim, 1, max)
  
  plot(tvs - tmin, vv, type="n", axes=FALSE, xlab=expression(t), 
       ylab=expression(Y(t)), ylim=range(c(lo, up, vv), na.rm = TRUE))
  axis(1); axis(2)
  polygon(c(tvs - tmin, rev(tvs - tmin)), c(lo, rev(up)), 
          col="grey80", border="grey80")
  lines(tvs - tmin, rep(0, length(tvs)), col=3, lty=3)
  lines(tvs - tmin, vv)
}
