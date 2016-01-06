
catalog <- function(data, time.begin=NULL, study.start=NULL,
                    study.end=NULL, study.length=NULL,
                    lat.range=NULL, long.range=NULL,
                    region.poly=NULL, mag.threshold=NULL,
                    flatmap=TRUE, tz="GMT")
{
  data <- as.data.frame(data)
  dnames <- tolower(names(data))
  names(data) <- dnames
  vnames <- c("date", "time", "long", "lat", "mag")
  if (!all(vnames %in% dnames))
    stop(paste("argument", sQuote(data),
               "must be a data frame with column names ",
               toString(sQuote(vnames))))
  if (any(is.na(data[, vnames])))
    stop(paste(sQuote(vnames), "must not contain NA values"))
  if (!is.numeric(data$lat) || !is.numeric(data$long) || !is.numeric(data$mag))
    stop("lat, long and mag columns must be numeric vectors")

  # extract spatial coordinates and magnitude
  xx <- data$long  # longitude of epicenter: coordinates
  yy <- data$lat   # latitude of epicenter : coordinates
  mm <- data$mag   # magnitude

  # extract date and time of events
  dt <- as.POSIXlt(paste(data$date, data$time), tz=tz)
  if (any(duplicated(dt)))
    stop(paste("no more than one event can occur simultaneously!",
               "check events", toString(which(duplicated(dt)))))
  if (is.unsorted(dt))
  {
    warning(paste("events were not chronologically sorted:",
                  "they have been sorted in ascending order"))
    data <- data[order(dt), ]
    dt <- sort(dt)
  }

  if (is.null(time.begin))
    time.begin <- min(dt)
  else
  {
    time.begin <- as.POSIXlt(time.begin, tz=tz)
    if (all(dt < time.begin))
      stop(paste("change time.begin: no event has occurred after:",
                 sQuote(time.begin)))
  }
  if (is.null(study.start))
    study.start <- time.begin
  else
  {
    study.start <- as.POSIXlt(study.start, tz=tz)
    if (study.start < time.begin)
      stop(paste("study.start", sQuote(study.start),
                 "can not be set before time.begin =", sQuote(time.begin)))
  }
  if (!is.null(study.length))
  {
    if (!is.null(study.end))
      stop("eihter study.end or study.length needs to be specified, not both")
    else
    {
      if(!is.numeric(study.length) || length(study.length) > 1)
        stop(paste("study.length must be single numeric value: in deciaml days"))
      study.end <- study.start + study.length * 24 * 60 * 60
    }
  }
  if (is.null(study.end))
    study.end <- max(dt)
  else
  {
    study.end <- as.POSIXlt(study.end, tz=tz)
    if (study.end < study.start)
      stop(paste("study.end", sQuote(study.end),
                 "can not be set before study.start", sQuote(study.start)))
  }
  tt <- date2day(dt, time.begin, tz=tz)

  # spatial region
  if (is.null(lat.range))
    lat.range <- range(yy)
  else if (!is.vector(lat.range) || length(lat.range) != 2 || lat.range[2] <= lat.range[1])
    stop("lat.range must be a vector of length 2 giving (lat.min, lat.max)")
  if (is.null(long.range))
    long.range <- range(xx)
  else if (!is.vector(long.range) || length(long.range) != 2 || long.range[2] <= long.range[1])
    stop("long.range must be a vector of length 2 giving (long.min, long.max)")

  if (is.null(region.poly))
  {
    region.poly <- list(long=c(long.range, rev(long.range)),
                        lat=rep(lat.range, each=2))
    region.win <- spatstat::owin(xrange=long.range, yrange=lat.range)
  }
  else
  {
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
    region.win <- spatstat::owin(poly=list(x=region.poly$long, y=region.poly$lat))
    region.area <- spatstat::area.owin(region.win) #Area.xypolygon(list(x=region.poly$long, y=region.poly$lat))
    if (region.area < 0)
      stop(paste("Area of polygon is negative -",
                 "maybe traversed in wrong direction?"))
  }

  # magnitude threshold
  if (is.null(mag.threshold))
    mag.threshold <- min(data$mag)
  else if (!is.numeric(mag.threshold) || length(mag.threshold) > 1)
    stop("mag.threshold must be a single numeric value")

  # project long-lat coordinates to flat map coordinates
  longlat.coord <- data.frame(long=xx, lat=yy)
  if (flatmap)
  {
    ymean <- spatstat::centroid.owin(region.win)$y
    theta <- cos(ymean * pi / 180)
    A <- matrix(c(theta, 0, 0, 1), ncol=2, nrow=2)
    region.win <- spatstat::affine(region.win, mat=A, vec=c(0, 0), rescue=TRUE)
    xx <- theta * xx
  }

  ok <- (dt <= study.end) & (dt > time.begin) & (mm >= mag.threshold)
  xx <- xx[ok]
  yy <- yy[ok]
  tt <- tt[ok]
  mm <- mm[ok] - mag.threshold
  flag <- as.integer(spatstat::inside.owin(xx, yy, region.win))
  flag[dt[ok] < study.start] <- -2
  revents <- cbind(tt, xx, yy, mm, flag, 0, 1, 0)
  longlat.coord <- longlat.coord[ok, ]
  longlat.coord$flag <- flag
  longlat.coord$dt <- dt[ok]
  X <- spatstat::ppx(data.frame(t=tt, x=xx, y=yy, m=mm),
           coord.type=c("t", "s", "s", "m"))

  switch(region.win$type, polygonal= {
    px <- region.win$bdry[[1]]$x
    py <- region.win$bdry[[1]]$y
  }, rectangle = {
    px <- c(region.win$xrange, rev(region.win$xrange))
    py <- rep(region.win$yrange, each=2)
  })
  # repeat the first vertex
  np <- length(px) + 1
  px[np] <- px[1]
  py[np] <- py[1]
  rpoly <- cbind(px, py)

  rtperiod <- c(date2day(study.start, time.begin, tz=tz),
                date2day(study.end, time.begin, tz=tz))
  out <- list(revents=revents, rpoly=rpoly, rtperiod=rtperiod, X=X,
              region.poly=region.poly, region.win=region.win,
              time.begin=time.begin, study.start=study.start,
              study.end=study.end, study.length=study.length,
              mag.threshold=mag.threshold, longlat.coord=longlat.coord)
  class(out) <- "catalog"
  return(out)
}


print.catalog <- function (x, ...)
{
  cat("earthquake catalog:\n  time begin", as.character(x$time.begin),
      "\n  study period:", as.character(x$study.start),
      " to ", as.character(x$study.end), "(", x$rtperiod[2], "days)")
  cat("\ngeographical region:\n  ")
  switch(x$region.win$type, rectangle={
    cat("  rextangular = [", x$region.poly$long[1], ",", x$region.poly$long[2],
        "] x [", x$region.poly$lat[1], x$region.poly$lat[2], "]\n")
  }, polygonal={
    cat("  polygonal with vertices:\n")
    print(cbind(lat=x$region.poly$lat, long=x$region.poly$long))
  })
  cat("threshold magnitude:", x$mag.threshold)
  cat("\nnumber of events:\n  total events", nrow(x$revents),
      ":", sum(x$revents[, 5] == 1), "target events, ",
      sum(x$revents[, 5] != 1), "complementary events\n  (",
      sum(x$revents[, 5] == 0), "events outside geographical region,",
      sum(x$revents[, 5] == -2), "events outside study period)")
}

plot.catalog <- function(x, ...)
{
  oldpar <- par(no.readonly = TRUE)
  par(mfrow=c(2, 2), mar=c(4, 4.1, 1, 1))
  plot(x$longlat.coord$long, x$longlat.coord$lat, xlab="long", ylab="lat", col=8,
       cex=2 * (x$revents[, 4] + 0.1)/max(x$revents[, 4]), asp=TRUE)
  ok <- x$revents[, 5] == 1
  points(x$longlat.coord$long[ok], x$longlat.coord$lat[ok], col=4,
         cex=2 * (x$revents[ok, 4] + 0.1)/max(x$revents[ok, 4]))
  polygon(x$region.poly$long, x$region.poly$lat, border=2)
  mbk <- seq(0, max(x$revents[, 4]), 0.1) + x$mag.threshold
  mct <- cut(x$revents[, 4] + x$mag.threshold, mbk)
  plot(mbk[-length(mbk)], as.numeric(table(mct)), type="b",
       xlab="mag", ylab="number of events", axes=FALSE)
  axis(1); axis(2)
  tbk <- seq(0, max(x$revents[, 1]), l=100)
  tct <- cut(x$revents[, 1], tbk)
  plot(tbk[-length(tbk)], cumsum(table(tct)), type="l",
       xlab="time", ylab="number of events", axes=FALSE)
  axis(1); axis(2)
  abline(v=x$rtperiod[1], col=2, lty=2)
  abline(v=x$rtperiod[2], col=2, lty=2)
  plot(x$revents[, 1], x$revents[, 3], xlab="time", ylab="lat",
       cex=2 * (x$revents[, 4] + 0.1)/max(x$revents[, 4]), col=8)
  points(x$revents[ok, 1], x$revents[ok, 3], col=4,
         cex=2 * (x$revents[ok, 4] + 0.1)/max(x$revents[ok, 4]))
  par(oldpar)
}
