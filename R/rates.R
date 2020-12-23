
probs <- function(fit)
{
  spatstat.geom::verifyclass(fit, "etas")
  object <- fit$object

  xx <- object$longlat$long
  yy <- object$longlat$lat

  # probability of each event being a triggered event
  pb <- 1 - object$revents[, 7]
  return(data.frame(long=xx, lat=yy, prob=pb, target=object$revents[, 5] == 1))
}


rates <- function(fit, lat.range=NULL, long.range=NULL,
                  dimyx=NULL, plot.it=TRUE)
{
  spatstat.geom::verifyclass(fit, "etas")
  rates.inter(fit$param, fit$object, fit$bwd, lat.range=lat.range,
               long.range=long.range, dimyx=dimyx, plot.it=plot.it)
}

rates.inter <- function(theta, object, bwd, lat.range=NULL, long.range=NULL,
                        dimyx=NULL, plot.it=TRUE)
{
  if (is.null(lat.range))
  {
    lat.range <- range(object$region.poly$lat)
  }
  if (is.null(long.range))
  {
    long.range <- range(object$region.poly$long)
  }
  xy.bnd <- longlat2xy(c(long.range, rev(long.range)),
                      rep(lat.range, each=2),
                      object$region.poly, object$dist.unit)
  if (is.null(dimyx))
  {
    rv <- diff(range(xy.bnd$x)) / diff(range(xy.bnd$y))
    if (rv > 1)
    {
      dimyx <- round(128 * c(1, rv))
    } else
    {
      dimyx <- round(128 * c(1 / rv, 1))
    }
  }

  gx <- seq(min(xy.bnd$x), max(xy.bnd$x), length.out=dimyx[2])
  gy <- seq(min(xy.bnd$y), max(xy.bnd$y), length.out=dimyx[1])
  out <- cxxrates(theta, object$revents, bwd, object$rtperiod, gx, gy)

  out <- list(x=seq(long.range[1], long.range[2], length.out=dimyx[2]),
              y=seq(lat.range[1], lat.range[2], length.out=dimyx[1]),
              bkgd=out$bkgd, total=out$total,
              clust=out$clust, lamb=out$lamb)

  if (plot.it)
  {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(2, 2), mar=c(2, 2, 2.25, 5))
    fields::image.plot(out$x, out$y, out$bkgd, asp=TRUE, xlab="",
                       ylab="", main="background seismicity rate",
                       legend.width=1, legend.shrink=0.9)
    maps::map('world', add=TRUE, col="grey50")
    polygon(object$region.poly$long, object$region.poly$lat, border=2)

    fields::image.plot(out$x, out$y, out$total, asp=TRUE, xlab="",
                       ylab="", main="total spatial seismicity rate",
                       legend.width=1, legend.shrink=0.9)
    maps::map('world', add=TRUE, col="grey50")
    polygon(object$region.poly$long, object$region.poly$lat, border=2)

    fields::image.plot(out$x, out$y, out$clust, asp=TRUE, xlab="",
                       ylab="", main="clustering coefficient",
                       legend.width=1, legend.shrink=0.9)
    maps::map('world', add=TRUE, col="grey50")
    polygon(object$region.poly$long, object$region.poly$lat, border=2)

    fields::image.plot(out$x, out$y, out$lamb, asp=TRUE, xlab="",
                       ylab="", main="",
                       legend.width=1, legend.shrink=0.9)
    title(main="conditional intensity function\nat the end of study period")
    maps::map('world', add=TRUE, col="grey50")
    polygon(object$region.poly$long, object$region.poly$lat, border=2)
    par(oldpar)
    invisible(out)
  }
  else
    return(out)
}
