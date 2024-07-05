
simetas <- function(param, bkgd, sim.start, sim.end=NULL, sim.length=NULL, 
                    lat.range=NULL, long.range=NULL, region.poly=NULL, 
                    mag.threshold=NULL,
                    flatmap=TRUE, dist.unit="degree", cat0=NULL, 
                    roundoff=FALSE, tz="GMT")
{
  sim.start <- as.POSIXlt(sim.start, tz=tz)
  if (!is.null(sim.length))
  {
    if (!is.null(sim.end))
      stop("eihter sim.end or sim.length needs to be specified, not both")
    else
    {
      if(!is.numeric(sim.length) || length(sim.length) > 1)
        stop(paste("study.length must be single numeric value: in deciaml days"))
      sim.end <- sim.start + sim.length * 24 * 60 * 60
    }
  } else {
    sim.end <- as.POSIXlt(sim.end, tz=tz)
    if (sim.end < sim.start)
      stop(paste("sim.end", sQuote(sim.end),
                 "can not be set before sim.start", sQuote(sim.start)))
  }
  simtt <- date2day(c(sim.start, sim.end), time.begin=NULL, tz=tz)
  
  if (is.null(region.poly))
  {
    region.poly <- list(long=c(long.range, rev(long.range)),
                        lat=rep(lat.range, each=2))
    region.win <- spatstat.geom::owin(xrange=long.range, yrange=lat.range)
  } else{
    if (is.data.frame(region.poly))
      region.poly <- as.list(region.poly)
    if (!is.list(region.poly) || !all(c("lat", "long") %in% names(region.poly)))
      stop("region.poly must be a list with components lat and long")
    if (any(is.na(region.poly$lat)) || any(is.na(region.poly$long)))
      stop("lat and long coordinates must not contain NA values")
    if (!is.numeric(region.poly$lat) || !is.numeric(region.poly$long) ||
        length(region.poly$lat) != length(region.poly$long))
      stop("lat and long coordinates must be numeric vectors of equal length")
    if (length(region.poly$lat) < 3)
      stop("region.poly needs at least 3 vertices")
    region.win <- spatstat.geom::owin(poly=list(x=region.poly$long, y=region.poly$lat))
    region.area <- spatstat.geom::area.owin(region.win) #Area.xypolygon(list(x=region.poly$long, y=region.poly$lat))
    if (region.area < 0)
      stop(paste("Area of polygon is negative -",
                 "maybe traversed in wrong direction?"))
  }
  
  U <- spatstat.geom::as.im(list(x=bkgd$x, y=bkgd$y, z=bkgd$bkgd))
  U <- U * diff(simtt)
  X0 <- spatstat.random::rpoispp(U[region.win, drop=FALSE])
  if (X0$n > 0)
  {
    tt <- runif(X0$n, min=simtt[1], max=simtt[2])
    mm <- rexp(n=X0$n, rate=param["beta"])
  } else {
    tt <- NULL
    mm <- NULL
  }
  # project long-lat coordinates to flat map coordinates
  if (flatmap)
  {
    proj <- longlat2xy(long=X0$x, lat=X0$y, region.poly=region.poly,
                       dist.unit=dist.unit)
    xx <- proj$x
    yy <- proj$y
    region.win <- proj$region.win
  }
  
  
  revents <- data.frame(tt=tt, xx=xx, yy=yy, mm=mm)
  
  l <- 1
  out <- vector("list", length=1000)
  out[[1]] <- revents 
  while (TRUE)
  {
    l <- l + 1
    out[[l]] <- NULL
    for (i in 1:nrow(out[[l - 1]]))
    {
      nl <- rpois(n=1, lambda=param["A"] * 
                    exp(param["alpha"] * out[[l - 1]]$mm[i]))
      if (nl > 0)
      {
        tl <- out[[l - 1]]$tt[i] + 
          param["c"] + param["c"] / ((1 - runif(nl))^(param["p"] - 1))
        sig <- param["D"] * exp(param["gamma"] * out[[l - 1]]$mm[i])
        rl <- sqrt(sig + sig / ((1 - runif(nl))^(param["q"] - 1)))
        th <- runif(nl, 0, 2 * pi)
        xl <- out[[l - 1]]$xx[i] + rl * cos(th)
        yl <- out[[l - 1]]$yy[i] + rl * sin(th)
        out[[l]] <- rbind(out[[l]], data.frame(tt=tl, xx=xl, yy=yl, 
                                               mm=rexp(nl, rate=param["beta"])))
        
      }
    }
    
    if (is.null(out[[l]]))
    {
      break
    }
  }
  
  L <- min(which(unlist(lapply(out, is.null))))
  if (L > 1)
  {
    for (i in 2:(L - 1))
      out[[1]] <- rbind(out[[1]], out[[i]])
  }
  out <- out[[1]]
  
  out <- out[order(out$tt), ]
  out$tt <- sim.start + out$tt * 24 * 60 * 60
  if (flatmap)
  {
    proj <- xy2longlat(x=out$xx, y=out$yy, region.poly=region.poly,
                       dist.unit=dist.unit)
    out$long <- proj$long
    out$lat <- proj$lat
  }
  
  if (is.null(mag.threshold))
    mag.threshold <- 0
  out$mag <- out$mm + mag.threshold
  out$date <- substr(as.character(out$tt), 1, 10)
  out$time <- substr(as.character(out$tt), 12, 19)
  data <- data.frame(date=out$date, time=out$time, 
                     long=out$long, lat=out$lat, mag=out$mag)
  simcat <- catalog(data, study.start=NULL, study.end=sim.end, 
                    lat.range=lat.range, long.range=long.range,
                    region.poly=region.poly, mag.threshold=mag.threshold, 
                    flatmap=flatmap, dist.unit=dist.unit, 
                    roundoff=roundoff, tz=tz)
  attr(simcat, "data") <- data
  return(simcat)
  
}
