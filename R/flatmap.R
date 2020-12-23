# Equirectangular projection

longlat2xy <- function(long, lat, region.poly, dist.unit="degree")
{
  if (length(long) != length(lat))
    stop(paste(sQuote("long"), "and", sQuote("lat"), "must be of equal length."))
  switch(dist.unit, degree={
    region.win <- spatstat.geom::owin(poly=list(x=region.poly$long,
                                           y=region.poly$lat))
    region.cnt <- spatstat.geom::centroid.owin(region.win)
    x <- cos(region.cnt$y * pi / 180) * (long - region.cnt$x)
    y <- lat - region.cnt$y
    px <- cos(region.cnt$y/180 * pi) * (region.poly$long - region.cnt$x)
    py <- region.poly$lat - region.cnt$y
  }, km={
    x <- 111.320 * cos(lat/180 * pi) * long
    y <- 110.574 * lat
    px <- 111.320 * cos(region.poly$lat/180 * pi) * region.poly$long
    py <- 110.574 * region.poly$lat
  }, stop(paste(sQuote("dist.unit"), "argument must be eighter",
                sQuote("degree"), "or", sQuote("km"), ".")))
  region.win <- spatstat.geom::owin(poly=list(x=px, y=py))
  list(x=x, y=y, region.win=region.win)
}

xy2longlat <- function(x, y, region.poly, dist.unit="degree")
{
  if (length(x) != length(y))
    stop(paste(sQuote("x"), "and", sQuote("y"), "must be of equal length."))
  switch(dist.unit, degree={
    region.win <- spatstat.geom::owin(poly=list(x=region.poly$long,
                                           y=region.poly$lat))
    region.cnt <- spatstat.geom::centroid.owin(region.win)
    long <- x / cos(region.cnt$y * pi / 180) + region.cnt$x
    lat <- y + region.cnt$y
  }, km={
    lat <- y / 110.574
    long <- x / (111.320 * cos(lat/180 * pi))
  }, stop(paste(sQuote("dist.unit"), "argument must be eighter",
                sQuote("degree"), "or", sQuote("km"), ".")))
  list(long=long, lat=lat)
}
