
search.isc <- function(start.year=1900, start.month=1, start.day=01,
                       end.year=2018, end.month=12, end.day=31,
                       searchshape="RECT",
                       lat.bot=NULL, lat.top=NULL,
                       long.left=NULL, long.right=NULL,
                       lat.ctr=NULL, long.ctr=NULL, 
                       radius=NULL, dist.units="deg",
                       dep.min=0, dep.max=100, nulldep=TRUE,
                       mag.min=4, mag.max=NULL, 
                       mag.type='MB', mag.agcy=NULL,
                       mirror=TRUE)
{
  servers <- c("http://www.isc.ac.uk",
               "http://isc-mirror.iris.washington.edu")
  dlink <- ifelse(mirror, servers[2], servers[1])
  dlink <- paste(dlink, "/cgi-bin/web-db-v4?request=COMPREHENSIVE",
                 "&out_format=CATCSV", sep="")
  
  if (!any(searchshape %in% c("RECT", "CIRC")))
    stop("searchshape must be either \'RECT\' or \'CIRC\'")
  
  if (!any(dist.units %in% c("deg", "km")))
    stop("dist.units must be either \'deg\' or \'km\'")
  
  dlink <- paste(dlink, "&searchshape=", searchshape,
                 "&bot_lat=", lat.bot, "&top_lat=", lat.top, 
                 "&left_lon=", long.left, "&right_lon=", long.right,
                 "&ctr_lat=", lat.ctr, "&ctr_lon=", long.ctr,
                 "&radius=", radius, 
                 "&max_dist_units=", dist.units,
                 "&srn=&grn=",
                 "&start_year=", start.year, 
                 "&start_month=", start.month, 
                 "&start_day=", start.day, 
                 "&start_time=00%3A00%3A00",
                 "&end_year=", end.year, 
                 "&end_month=", end.month, 
                 "&end_day=", end.day, 
                 "&end_time=00%3A00%3A00",
                 "&min_dep=", dep.min, "&max_dep=", dep.max, sep='')

  if (nulldep)
    dlink <- paste(dlink, "&null_dep=on", sep='')
  
  dlink <- paste(dlink, "&min_mag=", mag.min, 
                 "&max_mag=", mag.max, 
                 "&req_mag_type=", mag.type, 
                 "&req_mag_agcy=", mag.agcy,  sep='')
  
  download.file(dlink, destfile="~/iscdata.txt")
  allcontent <- readLines("~/iscdata.txt")
  if (length(allcontent) <= 27)
    stop("could not download the data: please try agian")
  lstart <- grep("--EVENT--", allcontent)
  lend <- grep("STOP", allcontent)
  skip <- allcontent[(lstart +1):(lend - 2)]
  vmax <- max(unlist(lapply(strsplit(skip, ","), length)))
  cnames <- c("EVENTID", "AUTHOR", "DATE", "TIME", "LAT", "LON", 
              "DEPTH", "DEPFIX", "AUTHOR", "TYPE", "MAG")
  cnames <- c(cnames, rep(c("AUTHOR", "TYPE", "MAG"), (vmax - 11)/3))
  mydata <- read.table(textConnection(skip[-1]), sep=",",
                       col.names=cnames, fill=TRUE)
  idx <- (1:length(cnames))[-c(3:6, 11)]
  
  qdata <- data.frame(date=mydata$DATE, time=mydata$TIME, 
                      lat=mydata$LAT, long=mydata$LON, 
                      mag=mydata$MAG, mydata[, idx])
  na <- apply(qdata[, 1:5], 1, function(a){ any(is.na(a)) })
  qdata <- qdata[!na, ]
  if (searchshape == "RECT")
  {
    ok <- (qdata$long >= long.left) & (qdata$long <= long.right)
    ok <- ok & (qdata$lat >= lat.bot) & (qdata$lat <= lat.top)
    qdata <- qdata[ok, ]
  }
  
  return(qdata)
}
