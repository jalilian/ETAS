
# Telesca, L., Lovallo, M., Golay, J., & Kanevski, M. (2016). 
# Comparing seismicity declustering techniques by means of 
# the joint use of Allan Factor and Morisita index. 
# Stochastic environmental research and risk assessment, 30(1), 77.

morisitaindex <- function(object, K=11, bwd=NULL, dimyx=NULL, cat.name=NULL)
{
  if (is.null(cat.name))
    cat.name <- deparse(substitute(object))
  
  ok <- object$revents[, "flag"] == 1
  xx <- object$revents[ok, "xx"]
  yy <- object$revents[ok, "yy"]
  N <- length(xx)
  win <- object$region.win
  areaW <- spatstat.geom::area.owin(win)
  X <- spatstat.geom::ppp(xx, yy, window=win)
  Lam <- Smooth.catalog(object, bwd=bwd, dimyx=dimyx)
  X.sim <- spatstat.core::rpoint(X$n, Lam, win=win, nsim=100)

  k.add <- 1
  delta <- areaW / (1:K + k.add)^2
  mifun <- function(Y)
  {
    mi <- numeric(K)
    for (k in 1:K)
    {
      #print(c(k=k, n=Y$n))
      Q <- spatstat.geom::quadratcount.ppp(Y, nx=k + k.add, ny=k + k.add)
      mi[k] <- (k + k.add)^2 * sum(Q * (Q - 1)) / (N * (N - 1))
    }
    return(mi)
  }
  obsmi <- log10(mifun(X))
  misim <- lapply(X.sim, mifun)
  
  misim.vals <- matrix(log10(unlist(misim)), ncol=100)
  q025 <- apply(misim.vals, 1, stats::quantile, p=0.025)
  q975 <- apply(misim.vals, 1, stats::quantile, p=0.975)
  ylim <- range(c(obsmi, q025[is.finite(q025)], q975), na.rm = TRUE)
  #oldpar <- par(no.readonly = TRUE)
  par(mar=c(4, 4.2, 1, 0.5))
  plot(log10(delta), obsmi, type="n", ylim=ylim, axes=FALSE, 
       xlab=expression(log[10]~delta),
       ylab=expression(log[10]~I(delta)), main=cat.name)
  axis(1); axis(2)
  graphics::polygon(c(log10(delta), rev(log10(delta))), 
                    c(q025, rev(q975)), col="grey70", border="grey70")
  graphics::lines(log10(delta), obsmi, lty=1.25)
  #par(oldpar)
  invisible(list(X=X, Lam=Lam))
}
