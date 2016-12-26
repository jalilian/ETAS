
rates <- function(fit, dimyx=NULL, plot.it=TRUE)
{
  spatstat::verifyclass(fit, "etas")
  rates.inter(fit$param, fit$object, fit$bwd, dimyx=dimyx, plot.it=plot.it)
}

rates.inter <- function(theta, object, bwd, dimyx=NULL, plot.it=TRUE)
{
  tperiod <- object$rtperiod
  xx <- object$revents[, "xx"]
  yy <- object$revents[, "yy"]
  flag <- object$revents[, 5]
  bk <- object$revents[, 6]
  pb <- object$revents[, 7]

  if (is.null(dimyx))
  {
    rv <- diff(range(xx[flag == 1])) / diff(range(yy[flag == 1]))
    if (rv > 1)
    {
      dimyx <- 256 * c(1, rv)
    } else
    {
      dimyx <- 256 * c(1 / rv, 1)
    }
  }
  if (!is.numeric(dimyx) || length(dimyx) != 2)
    stop(paste(sQuote(dimyx), "must be a numeric vector of length 2."))

  gr <- spatstat::gridcenters(object$region.win, dimyx[2], dimyx[1])
  gx <- gr$x
  gy <- gr$y

  r2 <- outer(xx, gx, FUN="-")^2 + outer(yy, gy, FUN="-")^2
  s1 <- exp(-r2/(2 * bwd^2)) / (2 * pi * bwd^2)
  s2 <- pb *  s1
  total <- colSums(s1)/diff(tperiod)
  bkgd <- colSums(s2)/diff(tperiod)
  clust <- 1 - bkgd/total

  proj <- xy2longlat(gx, gy, object$region.poly, object$dist.unit)
  gx <- proj$long
  gy <- proj$lat

  out <- data.frame(x=proj$long, y=proj$lat, total=total,
                    bkgd=bkgd, clust=clust)

  if (plot.it)
  {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(1, 2), mar=c(2, 2, 2, 4))
    fields::quilt.plot(out$x, out$y, out$bkgd,
                       main="background seismicity rate", asp=TRUE)
    maps::map('world', add=TRUE, col="grey50")
    fields::quilt.plot(out$x, out$y, out$clust,
                       main="clustering coefficient", asp=TRUE)
    maps::map('world', add=TRUE, col="grey50")
    par(oldpar)
    invisible(out)
  }
  else
    return(out)
}

probs <- function(fit)
{
  spatstat::verifyclass(fit, "etas")
  object <- fit$object

  xx <- object$longlat$long
  yy <- object$longlat$lat

  pb <- object$revents[, 7]
  return(data.frame(long=xx, lat=yy, prob=pb))
}
