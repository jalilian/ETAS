\name{date2day}
\alias{date2day}
\title{Convert date-time data to numeric data in decimal days}
\description{
  A function to convert date-time data to decimal days with respect to a date-time origin.
}
\usage{
  date2day(dates, start=NULL, tz="", \dots)
}
\arguments{
  \item{dates}{A date-time or date object. Typically, it is a character
    vector containing date-time information.}
  \item{start}{A date-time or date object. Determines the origin of the conversion.}
  \item{tz}{Optional. Timezone specification to be used for the conversion.}
  \item{\dots}{Arguments to be passed to \code{as.POSIXlt}.}
}
\value{
  A numeric vector of the same length as \code{dates}.
}
\details{
  The arguments \code{dates} and \code{start} must be of
  appropriate format to be passed to \code{as.POSIXlt} function.
}
\seealso{

  \code{\link{as.POSIXlt}} and \code{\link{difftime}} for appropriate format of the data
  to be converted.
}
\examples{
  # date-time data of Iran's earthquakes between 1973/01/01 and 2016/01/01
  dt <- paste(iran.quakes$date, iran.quakes$time)
  # origin of the conversion
  start <- "1973/01/01 00:00:00"
  # time in days since 1973-01-01 (UTC)
  date2day(dt, start, tz="GMT")
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{date time}
