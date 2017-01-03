\name{resid.etas}
\alias{resid.etas}
\title{Residuals Analysis and Diagnostics Plots}
\description{
  A function to compute and plot spatial and temporal residuals as well as
  transformed times for a fitted ETAS model.
}
\usage{
   resid.etas(fit, type="raw", n.temp=1000, dimyx=NULL)
}
\arguments{
  \item{fit}{A fitted ETAS model. An object of class \code{"etas"}.}
  \item{type}{A character string specifying the type residuals to be
  computed. Options are \code{"raw"} (the default case), \code{"reciprocal"}
  and \code{"pearson"}.}
  \item{n.temp}{An integer specifying the number of partition points for
  temporal residuals.}
  \item{dimyx}{Dimensions of the discretization for the smoothed
  spatial residuals. A numeric vector of length 2.}
}
\value{
  The function produces plots of temporal and smoothed spatial residuals,
  transformed times \eqn{\tau_j} against \eqn{j} and Q-Q plot of \eqn{U_j}.

  It also returns a list with components
  \itemize{
  \item\bold{tau} the transformed times
  \item\bold{U} related quantities with the transformed times
  \item\bold{tres} the temporal residuals
  \item\bold{sres} the smoothed spatial residuals
  }
}
\details{
  The function computes the temporal residuals
  \deqn{R^{temp}(I_j, h) = \sum_{i=1}^{N} \delta_i 1[t_i \in I_j] h(t_i)
    \lambda^{temp}(t_i|H_{t_i}) - \int_{I_j} h(t)\lambda^{temp}(t|H_t) d t}
  for \eqn{I_j=((j-1)T/n.temp, jT/n.temp]}, \eqn{j=1,...,n.temp},
  and the (smoothed version of) spatial residuals
  \deqn{R^{spat}(B_j, h) = h(\tilde{x}_i, \tilde{y}_i)
   \lambda^{spat}(\tilde{x}_i, \tilde{y}_i)(\tilde{\delta}_i - \tilde{w}_i)}
  for a Berman-Turner quadrature scheme with quadrature points
  \eqn{(\tilde{x}_i, \tilde{y}_i)} and quadrature weights \eqn{\tilde{w}_i},
  \eqn{i=1,...,n.spat}. Raw, reciprocal and Pearson residuals obtain
  with \eqn{h=1}, \eqn{h=1/\lambda} and \eqn{h=1/\sqrt{\lambda}},
  respectively.

  In addition, the function computes  transformed times
  \deqn{\tau_j=\int_{0}^{t_j} \lambda^{temp}(t|H_t) d t}
  and
  \deqn{U_j = 1 - \exp(-(t_j - t_{j-1}))}

}
\references{
  Baddeley A, Rubak E, Turner R (2015).
  \emph{ Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press, London. \url{http://www.crcpress.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/9781482210200/}.

  Baddeley A, Turner R (2000).
  Practical Maximum Pseudolikelihood for Spatial Point Patterns.
  \emph{Australian & New Zealand Journal of Statistics},
  \bold{42}(3), 283--322. \href{http://dx.doi.org/10.1111/1467-842X.00128}{doi:10.1111/1467-842X.00128}.

  Baddeley A, Turner R, Moller J, Hazelton M (2005).
  Residual Analysis for Spatial Point Processes.
  \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
  \bold{67}(5), 617--666. \href{http://dx.doi.org/10.1111/j.1467-9868.2005.00519.x}{doi:10.1111/j.1467-9868.2005.00519.x}.

  Ogata Y (1988).
  Statistical Models for Earthquake Occurrences and Residual Analysis for
  Point Processes.
  \emph{Journal of the American Statistical Association},
  \bold{83}(401), 9--27. \href{http://dx.doi.org/10.2307/2288914}{doi:10.2307/2288914}.

  Zhuang J (2006).
  Second-order Residual Analysis of Spatiotemporal Point Processes and
  Applications in Model Evaluation
    \emph{Journal of the Royal Statistical Society: Series B (Statistical
Methodology)},
  \bold{68}(4), 635--653. \href{http://dx.doi.org/10.1111/j.1467-9868.2006.00559.x}{doi:10.1111/j.1467-9868.2006.00559.x}.

}

\seealso{
  \code{\link{etas}}
}
\examples{

  iran.cat <- catalog(iran.quakes, time.begin="1973/01/01",
     study.start="1986/01/01", study.end="2016/01/01",
     lat.range=c(26, 40), long.range=c(44, 63), mag.threshold=5)
  print(iran.cat)
  \dontrun{
  plot(iran.cat)}

  # setting initial parameter values
  param0 <- c(0.46, 0.23, 0.022, 2.8, 1.12, 0.012, 2.4, 0.35)

  # fitting the model
  \dontrun{
  iran.fit <- etas(iran.cat, param0=param0)

  # diagnostic plots
  resid.etas(iran.fit)}

}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
