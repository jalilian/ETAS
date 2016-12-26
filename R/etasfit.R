
etasfit <- function(theta, revents, rpoly, tperiod, integ0, ihess,
                    verbose, ndiv, eps, cxxcode, nthreads)
{
  tht <- sqrt(theta)
  storage.mode(revents) <- storage.mode(rpoly) <- storage.mode(ihess) <- "double"
  if (cxxcode)
  {
    cfit <- cxxfit(tht, revents, rpoly, tperiod, integ0, ihess,
                   as.integer(ndiv),  eps, as.logical(verbose),
                   as.integer(nthreads))
  } else
  {
    rdata <- list(revents, rpoly, as.double(tperiod), as.double(integ0))
    cfit <- .Call("cfit", as.double(tht), rdata, ihess,
                  as.integer(verbose), PACKAGE="ETAS")
  }
  if(length(cfit) == 0)
    stop("Maximum Likelihood optimization failed to converge.\n
         Please try a better starting point.")

  if (!cxxcode)
    cfit[[5]] <- matrix(cfit[[5]], nrow=8, ncol=8)

  H <- cfit[[5]]
  tht <- cfit[[1]]

  avcov <- (1/4) * diag(1/tht) %*% H %*% diag(1/tht)
  list(estimate=tht^2, loglik=cfit[[2]], gradient=cfit[[3]],
       aic=cfit[[4]], ihessian=H, avcov=avcov)
}

