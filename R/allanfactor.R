
# Telesca, L., Lovallo, M., Golay, J., & Kanevski, M. (2016). 
# Comparing seismicity declustering techniques by means of 
# the joint use of Allan Factor and Morisita index. 
# Stochastic environmental research and risk assessment, 30(1), 77.

allanfactor <- function(object, K=100, nsim=500, cat.name=NULL)
{
  if (is.null(cat.name))
    cat.name <- deparse(substitute(object))

  ok <- object$revents[, "flag"] == 1
  tt <- object$revents[ok, "tt"]
  tmin <- object$rtperiod[1]
  tmax <- object$rtperiod[2]
  tau <- seq(10, (tmax - tmin)/10, length.out=K)
  
  affun <- function(times)
  {
    ac <- numeric(K)
    for (k in 1:K)
    {
      counts <- table(cut(times, seq(tmin, tmax, by=tau[k])))
      ac[k] <- mean(diff(counts)^2)/(2 * mean(counts))
    }
    return(ac)
  }
  obsaf <- log10(affun(tt))
  # simulating under Poisson hypothesis given the number of events
  simfun <- function(idx)
  {
    affun(sort(runif(length(tt), min=tmin, max=tmax)))
  }
  afsim <- lapply(1:nsim, simfun)
  afsim.mat <- matrix(log10(unlist(afsim)), ncol=nsim)
  q025 <- apply(afsim.mat, 1, stats::quantile, p=0.025)
  q975 <- apply(afsim.mat, 1, stats::quantile, p=0.975)
  ylim <- range(c(obsaf, q025[is.finite(q025)], q975), na.rm=TRUE)
  #oldpar <- par(no.readonly = TRUE)
  #par(mar=c(4, 4.2, 1, 0.5))
  plot(log10(tau), obsaf, type="n", ylim=ylim, axes=FALSE, 
       xlab=expression(log[10]~tau),
       ylab=expression(log[10]~AF(tau)), main=cat.name)
  axis(1); axis(2)
  polygon(c(log10(tau), rev(log10(tau))), c(q025, rev(q975)),
          col="grey70", border="grey70")
  abline(h=0, lty=2, col="green")
  graphics::lines(log10(tau), obsaf, lty=1.25)
  #par(oldpar)
}
