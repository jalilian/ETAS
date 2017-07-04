
resid.etas <- function(fit, type="raw", n.temp=1000, dimyx=NULL)
{
  flg <- fit$object$revents[, "flag"]
  tt <- fit$object$revents[flg == 1, "tt"]
  xx <- fit$object$revents[flg == 1, "xx"]
  yy <- fit$object$revents[flg == 1, "yy"]

  tau <- timetransform(fit)[flg == 1] - lambdatemporal(fit$object$rtperiod[1], fit)
  tg <- seq(min(tt), max(tt), length.out=n.temp)
  tlam <- lambdatemporal(tg, fit)
  dfun <- function(i){ sum((tt <= tg[i]) & (tt > tg[i - 1])) }
  tres <- switch(type, raw = unlist(lapply(2:n.temp, dfun))- tlam[-1] * diff(tg),
                 reciprocal = 1/tlam[-1] - diff(tg),
                 pearson = 1/sqrt(tlam[-1]) - sqrt(tlam[-1]) * diff(tg))

  W <- fit$object$region.win
  Xs <- spatstat::ppp(xx, yy, window=W, check=FALSE)
  qd <- spatstat::quadscheme(Xs)
  xg <- spatstat::x.quad(qd)
  yg <- spatstat::y.quad(qd)
  wg <- spatstat::w.quad(qd)
  slam <- lambdaspatial(xg, yg, fit)
  zg <- spatstat::is.data(qd)

  sres <- switch(type, raw = zg - slam * wg,
                 reciprocal = zg/slam - wg,
                 pearson = zg/sqrt(slam) - sqrt(slam) * wg)

  if (is.null(dimyx))
  {
    rv <- diff(range(xg)) / diff(range(yg))
    if (rv > 1)
    {
      dimyx <- c(128, floor(rv * 128))
    } else
    {
      dimyx <- c(floor(128 / rv), 128)
    }
  }

  Xg <- spatstat::ppp(xg, yg, window=W, check=FALSE)
  spatstat::marks(Xg) <- sres
  sres <- spatstat::Smooth(Xg, dimyx=dimyx, sigma=mean(fit$bwd))
  gr <- expand.grid(x=sres$xcol, y=sres$yrow)
  proj <- xy2longlat(gr$x, gr$y, region.poly=fit$object$region.poly,
                     dist.unit=fit$object$dist.unit)
  sres <- data.frame(x=proj$long, y=proj$lat, z=c(t(sres$v)))
  #sres <- stats::na.omit(sres)

  oldpar <- par(no.readonly = TRUE)

  par(mfrow=c(2, 2), mar=c(3.1, 3.1, 1.5, 1.6))

  plot(tg[-1], tres, type="l", main=paste(type, "temporal residuals"),
       xlab="", ylab="", axes=FALSE)
  abline(h=0, lty=2, col=2)
  axis(1); axis(2)
  mtext("time", 1, 1.95, cex=0.85)
  mtext("residuals", 2, 1.95, cex=0.85)

  zmax <- max(abs(sres$z), na.rm=TRUE)
  fields::quilt.plot(sres$x, sres$y, sres$z, zlim=c(-zmax, zmax),
                     nx=dimyx[2], ny=dimyx[1], asp=TRUE,
                     main=paste(type, "spatial residuals"))
  maps::map('world', add=TRUE, col="grey50")
 # polygon(fit$object$region.poly$long, fit$object$region.poly$lat, border=2)

  plot(tau, type="l", main="transformed times", xlab="",  ylab="",
       asp=TRUE, axes=FALSE)
  graphics::abline(a=0, b=1, col=2)
  graphics::grid(); graphics::axis(1); graphics::axis(2)
  mtext("i", 1, 1.95, cex=0.85)
  mtext(expression(tau[i]), 2, 1.95, cex=0.85)

  U <- 1 - exp(-diff(tau))
  stats::qqplot(U, runif(max(1000, length(U))), axes=FALSE, xlab="",
                ylab="", asp=1, main=expression(Q-Q~plot~of~U[i]))
  mtext(expression(observed~quantiles), 1, 1.95, cex=0.85)
  mtext(expression(quantiles~of~italic(U)(0,1)), 2, 1.95, cex=0.85)
  abline(c(0, 1), col="red")
  graphics::grid(); graphics::axis(1); graphics::axis(2)
  par(oldpar)

  out <- list(tau=tau, U=U, tres=tres, sres=sres, type=type)
  invisible(out)
}
