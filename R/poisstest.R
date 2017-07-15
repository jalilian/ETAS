# the code for joint spatio-temporal test is obtained and 
# adapted from 
# http://statistics.berkeley.edu/~stark/Code/Quake/permutest.r
# based on
#  Luen, B., & Stark, P. B. (2012). 
#   Poisson tests of declustered catalogues. 
#   Geophysical journal international, 189(1), 691-700.

poiss.test <- function(object, which="joint", r=NULL, lambda=NULL, bwd=NULL,
                       dimyx=NULL, nsim=299, n.perm=1000, verbose=TRUE, 
                       nthreads=1)
{
  ok <- object$revents[, "flag"] == 1
  tt <- object$revents[ok, "tt"]
  xx <- object$revents[ok, "xx"]
  yy <- object$revents[ok, "yy"]
  
  switch(which, temporal={
    # res1 <- ks.test(diff(tt), pexp, rate=1/mean(diff(tt)))
    res1 <- ks.test(tt, punif, min=object$rtperiod[1], 
                    max=object$rtperiod[2])
    res2 <- goftest::ad.test(tt, punif, min=object$rtperiod[1], 
                             max=object$rtperiod[2])
    return(list(KS=res1, AD=res2))
  }, spatial={
    win <- object$region.win
    unitname <- paste(object$dist.unit, c("", "s"), sep="")
    X <- spatstat::ppp(xx, yy, window=win, unitname=unitname)
    if (is.null(lambda))
      lambda <- Smooth.catalog(object, bwd=bwd, dimyx=dimyx)
    X.sim <- spatstat::rpoint(X$n, lambda, win=win, nsim=nsim)
    X.sim <- lapply(X.sim, function(x) { x$window <- win; x })

    if (is.null(r))
    {
      rmax <- spatstat::rmax.rule("K", win) / 3
      r <- seq(0, rmax, length=200)
    }
    stat <- function(Y, r)
    {
      lamY <- Smooth.catalog(object, bwd=bwd, dimyx=dimyx)
      spatstat::Linhom(Y, lambda=lamY, r=r, correction="translate")
    }
    env <- spatstat::envelope(X, stat, r=r, savefuns=TRUE, use.theory=TRUE, 
                              savepatterns=TRUE, simulate=X.sim, nsim=nsim, 
                              nrank=round(0.02 * nsim))
    res1 <- spatstat::dclf.test(env, use.theory=TRUE)
    res2 <- spatstat::mad.test(env, use.theory=TRUE)
    return(list(X=X, lambda=lambda, env=env, DCLF=res1, MAD=res2))
  }, joint={
    # exclude ties 
    dok <- !(duplicated(xx) | duplicated(yy))
    xx <- xx[dok]
    yy <- yy[dok]
    # extract ranks (assume no ties)
    x.rank <- rank(xx)
    y.rank <- rank(yy)
    
    # number of events
    n <- length(x.rank)
    
    # empirical distribution of spatial ranks
    xy.upper <- t(outer(y.rank, y.rank, "<=")) %*% 
      outer(x.rank, x.rank, "<=")
    
    # Distance function
   # teststat <- cxxstpoisstest(x.rank, y.rank, xy.upper)
    teststat <- cxxstpoisstestMP(x.rank, y.rank, xy.upper, nthreads)
    
    # permuting the occurrence times of events 
    permustat <- rep(NA, n.perm)
    for(permu in 1:n.perm)
    {
      o <- sample(n)
      x.perm <- x.rank[o]
      y.perm <- y.rank[o]
      xy.perm <- xy.upper[o, o]
      
      #permustat[permu] <- cxxstpoisstest(x.perm, y.perm, xy.perm)
      permustat[permu] <- cxxstpoisstestMP(x.perm, y.perm, xy.perm, nthreads)
      if (verbose)
        cat(permu, "\t", permustat[permu], "\n")
    }
    
    p.value <- mean(permustat >= teststat)
    names(teststat) <- "teststat"
    names(n.perm) <- "number of permutations"
    out <- list(statistic = teststat, p.value = p.value,
                method="Permutation test for conditional 
                exchangeability of occurrence times of events",
                data.names=deparse(substitute(object)), 
                parameter=n.perm)
    class(out) <- "htest"
    return(out)
  }, stop("wrong which choice."))
}

Smooth.catalog <- function(object, type="spatial", bwd=NULL, bwm=NULL, 
                           nnp=NULL, dimyx=NULL, convert=FALSE)
{
  ok <- object$revents[, "flag"] == 1
  switch(type, temporal={
    tt <- object$revents[ok, "tt"]
    if (is.null(bwd))
      bwd <- "SJ-ste"
    else
      stopifnot(length(bwd) == 1)
    tdens <- density(tt, bw=bwd)
    plot(tdens, axes=FALSE, main="", xaxt = "n",
         ylab=expression(hat(lambda)(t)))
    idx <- round(seq(1, length(tt), length.out=15))
    years <- substr(object$longlat.coord$dt[ok][idx], 1, 4)
    axis(1, at=tt[idx], labels=years)
    axis(2)
    rug(tt, col="grey75")
  }, spatial={
    xx <- object$revents[ok, "xx"]
    yy <- object$revents[ok, "yy"]
    win <- object$region.win
    
    if (is.null(dimyx))
    {
      rv <- diff(win$xrange)/diff(win$yrange)
      npixel <- spatstat::spatstat.options("npixel")
      if (rv > 1)
      {
        dimyx <- round(npixel * c(1, rv))
      } else
      {
        dimyx <- round(npixel * c(1 / rv, 1))
      }
    }
    
    # bandwidths for smoothness and integration
    if (is.null(bwd))
    {
      if (is.null(nnp))
        nnp <- round(log(length(xx)))
      bwd <- spatstat::nndist.default(xx, yy, k=nnp)
      if (is.null(bwm))
        bwm <- quantile(bwd, probs=0.25)
      bwd <- pmax(bwd, bwm)
    }
    else
    {
      stopifnot(is.numeric(bwd), length(bwd) != length(xx))
    }
    
    gx <- seq(win$xrange[1], win$xrange[2], length.out=dimyx[2])
    gy <- seq(win$yrange[1], win$yrange[2], length.out=dimyx[1])
    out <- cxxSmooth(xx, yy, bwd, gx, gy, TRUE)$out
    out <- out / diff(object$rtperiod)
    
    if (convert)
    {
      gcoords <- expand.grid(gx, gy)
      gcoords <- xy2longlat(gcoords[, 1], gcoords[, 2], 
                            object$region.poly, 
                            dist.unit=object$dist.unit)
      gx <- gcoords$long[1:dimyx[2]]
      gy <- gcoords$lat[dimyx[2] * (1:dimyx[1]) - 1]
    }
    lambda <- spatstat::as.im.default(list(x=gx, y=gy, z=out))
    attr(lambda, 'bwd') <- bwd
    attr(lambda, 'nnp') <- nnp
    return(lambda)
  }, stop("only spatial or temporal smmothing is implemented"))
}
