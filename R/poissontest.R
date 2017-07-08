# the code for joint spatio-temporal test is obtained and 
# adapted from 
# http://statistics.berkeley.edu/~stark/Code/Quake/permutest.r
# based on
#  Luen, B., & Stark, P. B. (2012). 
#   Poisson tests of declustered catalogues. 
#   Geophysical journal international, 189(1), 691-700.

poisson.test <- function(object, which="joint", r=NULL, 
                         n.perm=1000, verbose=TRUE)
{
  ok <- object$revents[, "flag"] == 1
  tt <- object$revents[ok, "tt"]
  xx <- object$revents[ok, "xx"]
  yy <- object$revents[ok, "yy"]
  
  switch(which, temporal={
    # res <- ks.test(diff(tt), pexp, rate=1/mean(diff(tt)))
    res <- ks.test(tt, punif, min=object$rtperiod[1], 
                  max=object$rtperiod[2])
    return(res)
  }, spatial={
    win <- object$region.win
    X <- spatstat::ppp(xx, yy, window=win)
    rxy <- diff(win$xrange)/diff(win$yrange)
    if (rxy >= 1)
      dimyx <- c(128, ceiling(128 * rxy))
    else
      dimyx <- c(ceiling(128 / rxy), 128)
    

    Lam <- spatstat::density.ppp(X, dimyx=dimyx, diggle=TRUE, 
                                 sigma=spatstat::bw.diggle)
    Lam[Lam$v < 0] <- 0
    X.sim <- spatstat::rpoint(X$n, Lam, win=win, nsim=99)
    X.sim <- lapply(X.sim, function(x) { x$window <- win; x })
    env <- spatstat::envelope(X, spatstat::Linhom, 
                              correction="translate", r=r, 
                              global = TRUE, savefuns = TRUE, 
                              use.theory=TRUE, savepatterns=TRUE, 
                              simulate=X.sim)
    res <- spatstat::dclf.test(env)
    par(mar=c(4, 4.1, 1.5, 0.5))
    plot(env, legend=FALSE, axes=FALSE, main="")
    axis(1); axis(2)
    mtext(paste("pvalue =", round(res$p.value, 3)), 3, -1)
    return(res)
  }, joint={
    # extract ranks (assume no ties)
    x.rank <- rank(xx)
    y.rank <- rank(yy)
    
    # number of events
    n <- length(x.rank)
    
    # empirical distribution of spatial ranks
    xy.upper <- t(outer(y.rank, y.rank, "<=")) %*% 
      outer(x.rank, x.rank, "<=")
    
    # Distance function
    distfind <- function(x.rank, y.rank, xy.upper)
    {
      dfv <- 0
      xyz.temp <- matrix(0, n, n)
      for(k in 1:n)
      {
        xyz.temp <- xyz.temp + 
          (y.rank >= y.rank[k]) %*% t(x.rank >= x.rank[k])
        dist.matrix <- xyz.temp/n - xy.upper/n * k/n
        dfv <- max(dfv, abs(dist.matrix))
      }
      return(dfv)
    }
    
    teststat <- distfind(x.rank, y.rank, xy.upper)
    
    # permuting the occurrence times of events 
    permustat <- rep(NA, n.perm)
    for(permu in 1:n.perm)
    {
      o <- sample(n)
      x.perm <- x.rank[o]
      y.perm <- y.rank[o]
      xy.perm <- xy.upper[o, o]
      
      permustat[permu] <- distfind(x.perm,y.perm,xy.perm)
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
  })
}
