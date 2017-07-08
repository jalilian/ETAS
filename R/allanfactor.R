
allanfactor <- function(object, K=200)
{
  ok <- object$revents[, "flag"] == 1
  tt <- object$revents[ok, "tt"]
  tmin <- object$rtperiod[1]
  tmax <- object$rtperiod[2]
  tau <- seq(10, (tmax - tmin)/5, length.out=K)
  
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
  afsim <- lapply(1:100, simfun)
  q025 <- apply(matrix(log10(unlist(afsim)), ncol=100), 1, 
                stats::quantile, p=0.025)
  q975 <- apply(matrix(log10(unlist(afsim)), ncol=100), 1, 
                stats::quantile, p=0.975)
  ylim <- range(c(obsaf, q025[is.finite(q025)], q975), na.rm = TRUE)
  par(mar=c(4, 4.1, 1, 0.5))
  plot(log10(tau), obsaf, type="n", ylim=ylim, axes=FALSE, 
       xlab=expression(log[10]~tau),
       ylab=expression(log[10]~AF(tau)), main="")
  axis(1); axis(2)
  polygon(c(log10(tau), rev(log10(tau))), c(q025, rev(q975)),
          col="grey70", border="grey70")
  lines(log10(tau), obsaf, lty=1.25)
}
