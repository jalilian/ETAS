
decluster <- function(theta, rbwd, revents, rpoly, tperiod, ndiv, cxxcode, mver)
{
  storage.mode(revents) <- storage.mode(rpoly) <- "double"
  if (cxxcode)
  {
    cbkg <- cxxdeclust(theta, revents, rpoly, rbwd, tperiod, 
                       as.integer(ndiv), as.integer(mver))
  } else
  {
    tht <- sqrt(theta)
    cbkg <- .Call("cdeclust", as.double(tht), as.double(rbwd),
                      revents, rpoly, as.double(tperiod), PACKAGE="ETAS")
  }
  list(revents=cbkg[[1]], integ0=cbkg[[2]])
}

