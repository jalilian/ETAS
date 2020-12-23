
etas <- function(object, param0=NULL, bwd = NULL, nnp = 5, bwm = 0.05,
                 verbose = TRUE, plot.it = FALSE, ndiv = 1000,
                 no.itr = 11, rel.tol=1e-03, eps = 1e-06,
                 cxxcode = TRUE, nthreads = 1)
{
  ptm <- proc.time()
  spatstat.geom::verifyclass(object, "catalog")
  revents <-  object$revents
  rpoly <-    object$rpoly
  rtperiod <- object$rtperiod
  m0 <- object$mag.threshold
  win <- object$region.win

  if (nthreads > parallel::detectCores())
  {
    stop(paste("nthreads can not be greater than", parallel::detectCores(),
               "on this machine!"))
  }

  # initial prameter values
  if (is.null(param0))
  {
    mu0 <- nrow(revents)/(4 * diff(rtperiod) * spatstat.geom::area.owin(win))
    param0 <- c(mu=mu0, A=0.01, c=0.01, alpha=1, p=1.3, D=0.01, q=2,
                gamma=1)
    if (object$dist.unit == "km")
      param0["D"] <- 111^2 * param0["D"]
    if (verbose)
    {
      cat("using non-informative initial parameter values:\n")
      print(param0)
    }
    warning("the algorithm is very sensitive to the choice of starting point")
  }
  # bandwidths for smoothness and integration
  if (is.null(bwd))
  {
    if (object$dist.unit == "km")
      bwm <-  6371.3 * pi / 180 * bwm
    rbwd <- spatstat.geom::nndist.default(revents[, 2], revents[, 3], k=nnp)
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
  thetar <- asd <- matrix(NA, nrow=no.itr, ncol=8)
  par.names <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  names(param1) <- colnames(thetar) <- colnames(asd) <- par.names
  loglikfv <- numeric(no.itr)
  rownames(thetar) <- rownames(asd) <- names(loglikfv) <- paste("iteration", 1:no.itr)
  ihess <- diag(8)
  bk <- numeric(nrow(revents))

  for (itr in 1:no.itr)
  {
    cat("declustering:\n")
    bkgpbar <- utils::txtProgressBar(min=0, max=no.itr + 1 - itr, style=3)
    for (l in 1:(no.itr + 1 - itr))
    {
      bkg <- decluster(param1, rbwd, revents, rpoly, rtperiod, ndiv, cxxcode)
      revents <- bkg$revents
      utils::setTxtProgressBar(bkgpbar, l)
    }
    close(bkgpbar)
    integ0 <- bkg$integ0
    dbk <- bk - revents[, 6]
    bk <- revents[, 6]
    pb <- revents[, 7]
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
      rates.inter(param1, object, rbwd, plot.it=plot.it)
    }
    cat("estimating:\n")
    opt <- etasfit(param1, revents, rpoly, rtperiod, integ0, ihess,
                   verbose, ndiv, eps, cxxcode, nthreads)
    thetar[itr, ] <- opt$estimate
    loglikfv[itr] <- opt$loglik
    asd[itr, ] <- sqrt(diag(opt$avcov))
    ihess <- opt$ihess
    param1 <- thetar[itr, ]
    if (verbose)
    {
      cat("======================================================\n")
      cat("MLE:\n")
      print(param1)
      cat("======================================================\n")
    }
    if (itr > 1)
    {
      dtht <- max((thetar[itr, ] - thetar[itr - 1, ])/thetar[itr - 1, ])
      dlrv <- abs(loglikfv[itr] / loglikfv[itr - 1] - 1)
      dbkv <- max(abs(dbk / bk))
      print(c(dtht, dlrv, dbkv))
      if (all(c(dtht, dlrv, dbkv) < rel.tol))
        break
    } else{
      if (itr == no.itr)
        warning("Reached maximum number of iterations\n")
    }

  }

  if (verbose)
  {
    cat("Execution time:\n")
    print(proc.time() - ptm)
  }

  names(param1) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  object$revents <- revents
  out <- list(param = param1, bk=bk, pb=pb, opt=opt, object=object,
              bwd=rbwd, thetar=thetar, loglikfv=loglikfv, asd=asd,
              integ0=integ0, ndiv=ndiv, itr=itr,
              exectime=proc.time() - ptm)
  class(out) <- "etas"
  return(out)
}


print.etas <- function (x, ...)
{
  cat("ETAS model: fitted using iterative stochastic declustering method\n")
  cat("converged after", x$itr, "iterations: elapsed execution time",
      round(x$exectime[3]/60, 2), "minutes\n\n")
  mm <- x$object$revents[x$object$revents[, 5] == 1, 4]
  bt <- 1 / mean(mm)
  asd.bt <- bt^2 / length(mm)
  cat("ML estimates of model parameters:\n")
  ests <- cbind("Estimate" = c(beta=bt, x$param),
                "StdErr" = c(asd.bt, x$asd[x$itr, ]))
  print(round(t(ests), 4))
  cat("\nDeclustering probabilities:\n")
  print(round(summary(x$pb), 4))
  cat("\nlog-likelihood: ", x$opt$loglik, "\tAIC: ", x$opt$aic, "\n")
}


plot.etas <- function(x, which="est", dimyx=NULL, ...)
{
  switch(which, loglik={
    plot(x$loglikfv[1:x$itr], xlab="iterations", 
         ylab="log-likelihood", type="b",
         main="log-likelihood function of the model")
  }, est={
    theta.ts <- stats::ts(x$thetar[1:x$itr, ])
    lattice::xyplot(theta.ts, xlab="iteration", type="b", pch=16,
                    scales=list(x=list(tick.number=x$itr),
                                y=list(tick.number=2)),
                    main="estimates of the model parameters")
  }, dots={
    graphics::dotchart(x$thetar[1:x$itr, ],
                       main="estimates of the model parameters")
  }, rates={
    rates.inter(x$param, x$object, x$bwd)
  }, stop("Wrong type"))
}
