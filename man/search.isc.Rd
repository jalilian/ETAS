\name{search.isc}
\alias{search.isc}
\title{Retrieving data from the ISC Bulletin}
\description{
  Searches the International Seismological Centre (ISC) bulletin 
  for recorded eqrthqukes ocurred during a specified time 
  period in a specified geographical region (rectangular or circular) 
  with a specified magnitude type and magnitude threshould.
}
\usage{
  search.isc(start.year=1900, start.month=1, start.day=01,
             end.year=2018, end.month=12, end.day=31,
             searchshape="RECT", lat.bot=NULL, lat.top=NULL,
             long.left=NULL, long.right=NULL, lat.ctr=NULL, 
             long.ctr=NULL, radius=NULL, dist.units="deg", 
             dep.min=0, dep.max=100, nulldep=TRUE,
             mag.min=4, mag.max=NULL, mag.type='MB', 
             mag.agcy=NULL, mirror=TRUE)
}
\arguments{
  \item{start.year}{A single numeric value specifying the start 
  year of the time period.}
  \item{start.month}{A single numeric value (1 to 12) specifying 
  the start month of the time period.}
  \item{start.day}{A single numeric value (1 to 31) specifying 
  the start day of the time period.} 
  \item{end.year}{A single numeric value specifying the end 
  year of the time period.}
  \item{end.month}{A single numeric value (1 to 12) specifying 
  the end month of the time period.}
  \item{end.day}{A single numeric value (1 to 31) specifying 
  the end day of the time period.} 
  \item{searchshape}{A character string specifying the shape of
  the search region. Rectangular search with 
  \code{searchshape="RECT"} (default) and circular search 
  \code{searchshape="CIRC"} are possible.}
  \item{lat.bot}{A single numeric value (-90 to 90) specifying 
  the bottom latitude of the rectangular search region 
  (if \code{searchshape="RECT"}).}
  \item{lat.top}{A single numeric value (-90 to 90) specifying 
  the top latitude of the rectangular search region
  (if \code{searchshape="RECT"}).}
  \item{long.left}{A single numeric value (-180 to 180) specifying 
  the left longitude of the rectangular search region
  (if \code{searchshape="RECT"}).}
  \item{long.right}{A single numeric value (-180 to 180) specifying 
  the right longitude of the rectangular search region
  (if \code{searchshape="RECT"}).}
  \item{lat.ctr}{A single numeric value (-90 to 90) specifying 
  the central latitude of the circular search region
  (if \code{searchshape="CIRC"}).}
  \item{long.ctr}{A single numeric value (-180 to 180) specifying 
  the central longitude of the circular search region
  (if \code{searchshape="CIRC"}).}
  \item{radius}{A single numeric value specifying the search radius
  of the circular search region
  (if \code{searchshape="CIRC"}).}
  \item{dist.units}{A character string specifying the unit of
  the search radius: degrees with \code{"deg"} and kilometers 
  with \code{"km"} are the available options.}
  \item{dep.min}{A single numeric value specifying the minumum 
  depth (km).}  
  \item{dep.max}{A single numeric value specifying the maximum 
  depth (km).}
  \item{nulldep}{Logical flag indicating whether to include 
  events with unknown depths.}
  \item{mag.min}{A single numeric value specifying the minumum 
  magnitude.}
  \item{mag.max}{A single numeric value specifying the maximum 
  magnitude.}
  \item{mag.type}{A character string specifying the magnitude 
  type. Due to the large number of magnitude types and the 
  various notation used, the selected type will search all 
  magnitudes of that type. E.g. the default \code{"MB"} will 
  match all body-wave magnitudes: mb, mB, MB, mb1mx, etc.}
  \item{mag.agcy}{A character string specifying the magnitude 
  author.  \code{"ISC"}, \code{"NEIC"}, \code{"GCMT/HRVD"}, 
  \code{"JMA"}, amongst others are considered reliable magnitude
  authors. The default uses \code{"Any"}, which will match any 
  author that has a magnitude for both prime and secondary 
  hypocentres included in an event. Selecting 
  \code{"prime author"} will limit the search by magnitudes 
  directly associated with the prime hypocentre (usually the 
  same author as the prime hypocentre) - i.e. for events with 
  an ISC authored prime hypocentre, selecting 'Prime author', 
  will limit events to those with ISC determined magnitudes.}
  \item{mirror}{Logical flag indicating whether to use the ISC 
  mirror server hosted at IRIS DMC (\url{http://isc-mirror.iris.washington.edu}) as it may be faster.}
}
\value{
  An object of class \code{"data.frame"} which can be passed to 
  the \code{"catalog"} function in order to create an earthquake
  catalog dataset.
}
\details{

  Any use of data from the ISC by academic or commercial 
  organisations, as well as individuals should always be cited. 
  The correct format for citations may be found on 
  \url{http://www.isc.ac.uk/citations/}.

  The ISC is named as a valid data centre for citations within 
  American Geophysical Union (AGU) publications. As such, please 
  follow the AGU guidelines when referencing ISC data in one of 
  their journals. The ISC may be cited as both the institutional 
  author of the Bulletin and the source from which the data 
  were retrieved. 
}
\seealso{
  \code{catalog}.
}
\references{
  International Seismological Centre (2015). 
  \emph{On-line Bulletin}, 
  \url{http://www.isc.ac.uk}, 
  Internatl. Seismol. Cent., 
  Thatcham, United Kingdom.
  
  International Seismological Centre (2015). 
  \emph{Reference Event Bulletin}, 
  \url{http://www.isc.ac.uk}, 
  Internatl. Seismol. Cent., 
  Thatcham, United Kingdom.
  
  International Seismological Centre (2015). 
  \emph{EHB Bulletin}, 
  \url{http://www.isc.ac.uk}, 
  Internatl. Seismol. Cent., 
  Thatcham, United Kingdom.
}

\examples{

  
  # rectangular search
  \dontrun{
  mydata <- search.isc(start.year=1980, start.month=1, start.day=01,
                   end.year=2017, end.month=11, end.day=11, 
                   lat.bot=33, lat.top=37,
                   long.left=44, long.right=48,
                   mag.min=3.5, mag.type='MB')  
  mycatalog <- catalog(mydata, study.start = "1990-01-01")
  plot(mycatalog)
  }
  # circular search
  \dontrun{
  mydata2 <- search.isc(start.year=2017, start.month=11, start.day=12,
                   end.year=2018, end.month=3, end.day=01,
                   searchshape="CIRC", 
                   lat.ctr=34.905, long.ctr=45.956,
                   radius=150, dist.units="km",
                   mag.min=3.5, mag.type='ML') 
  mycatalog2 <- catalog(mydata2)
  plot(mycatalog2)
  }

}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
