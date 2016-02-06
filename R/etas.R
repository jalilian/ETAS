
etas <- function(object, param0, bwd = NULL, nnp = 5, bwm = 0.05,
                 verbose = TRUE, plot.it = FALSE, no.itr = 11)
{
  spatstat::verifyclass(object, "catalog")
  revents <-  object$revents
  rpoly <-    object$rpoly
  rtperiod <- object$rtperiod
  m0 <- object$mag.threshold
  win <- object$region.win

  # initial prameter values
  if (is.null(param0))
  {
    param0 <- c(mu=0.5, A=0.2, c=0.01, alpha=2, p=1.1, D=0.002, q=1.5, gamma=0.6)
    cat("using non-informative initial parameter values:\n")
    print(param0)
    warning("the algorithm is very sensitive to the choice of starting point")
  }
  # bandwidths for smoothness and integration
  if (is.null(bwd))
  {
    rbwd <- spatstat::nndist.default(revents[, 2], revents[, 3], k=nnp)
    rbwd <- pmax(rbwd, bwm)
  }
  else
  {
    stopifnot(is.numeric(bwd), length(bwd) != nrow(revents))
    rbwd <- bwd
  }

  # check initial values for the model parameters
  if (!is.numeric(param0) || length(param0) != 8 || any(param0 < 0))
    stop("param0 must be a numeric vector of length 8 with positive components")

  param1 <- param0
  names(param1) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  thetar <- matrix(NA, nrow=no.itr, ncol=8)

  for (itr in 1:no.itr)
  {
    bkg <- decluster(param1, rbwd, revents, rpoly, rtperiod)
    revents <- bkg$revents
    integ0 <- bkg$integ0
    bk <- revents[,6]
    pb <- revents[,7]
    if (verbose)
    {
      cat("iteration: ", itr, "\n")
      cat("======================================================\n")
      cat("background seismicity rate:\n")
      print(summary(bk))
      cat("probability of being a background event:\n")
      print(summary(pb))
      cat("integral of background seismicity rate: ", integ0, "\n")
      cat("======================================================\n")
    }
    if (plot.it)
    {
      par(mfrow=c(1, 2), mar=c(4, 4, 3, 1))
      cols <- ifelse(pb < 0.5, "red", "blue")
      plot(object$longlat$long, object$longlat$lat,
           cex = 0.05 + 2.5 * revents[, 4]/m0, col=cols,
           main=paste("iteration: ", itr), xlab="long", ylab="lat")
      polygon(object$region.poly$long, object$region.poly$lat, border=3)
      plot(revents[,1], pb, xlab="time",
           ylab="probability of being a background event")
      rates0 <- rates(param1, object, rbwd, plot.it=plot.it)
    }
    opt <- etasfit(param1, revents, rpoly, rtperiod, integ0, verbose)
    thetar[itr,] <- opt$estimate
    param1 <- thetar[itr,]
    if (verbose)
    {
      cat("======================================================\n")
      cat("MLE:\n")
      print(param1)
      cat("======================================================\n")
    }
  }

  rates0 <- rates(param1, object, rbwd, plot.it=plot.it)
  names(param1) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  out <- list(param = param1, bk=bk, pb=pb, opt=opt, rates=rates0)
  class(out) <- "etas"
  return(out)
}


print.etas <- function (x, ...)
{
  cat("ETAS model: fitted using iterative stochastic declustering method\n")
  cat("ML estimates of model parameters:\n")
  print(x$param)
  cat("log-likelihood: ", x$opt$loglik, "\nAIC: ", x$opt$aic, "\n")
}
