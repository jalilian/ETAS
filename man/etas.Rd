\name{etas}
\alias{etas}
\title{Fit the space-time ETAS model to data}
\description{
  A function to fit the space-time version of the Epidemic Type
  Aftershock Sequence (ETAS) model to a catalog of earthquakes
  (a spatio-temporal point pattern) and perform a stochastic
  declustering method.
}
\usage{
 etas(object, param0 = NULL, bwd = NULL, nnp = 5, bwm = 0.05,
      verbose = TRUE, plot.it = FALSE, ndiv = 1000, no.itr = 11,
      rel.tol=1e-03, eps = 1e-06, cxxcode = TRUE, nthreads = 1)
}
\arguments{
  \item{object}{An object of class \code{"catalog"} containing an
  earthquake catalog dataset.}
  \item{param0}{Initial guess for model parameters. A numeric vector
                of appropriate length (currently 8). See details.}
  \item{bwd}{Optional. Bandwidths for smoothness and integration
             on the geographical region \code{win}. A numeric vector
             which has the length of the number of events.
             If not supplied, the following arguments \code{nnp} and
             \code{bwm} determine bandwidths.}
  \item{nnp}{Number of nearest neighbors for bandwidth calculations. An integer.}
  \item{bwm}{Minimum bandwidth. A positive numeric value.}
  \item{verbose}{Logical flag indicating whether to print progress reports.}
  \item{plot.it}{Logical flag indicating whether plot probabilities of
    each event being a background event on a map.}
  \item{ndiv}{An integer indicating the number of knots on each
    side of the geographical region for integral approximation.}
  \item{no.itr}{An integer indicating the number of iterations for
   convergence of the iterative approach of simultaneous estimation
   and declustering algorithm. See details.}
  \item{rel.tol}{Relative iteration convergence tolerance of the
  iterative estimation approach.}
  \item{eps}{Optimization  convergence tolerance in the
  Davidon-Fletch-Powell algorithm}
  \item{cxxcode}{Logical flag indicating whether to use the C++ code.
  The C++ code is slightly faster and allows parallel computing.}
  \item{nthreads}{An integer indicating number of threads in
  the parallel region of the C++ code}
}
\value{
  A list with components
  \describe{
    \item{param:}{The ML estimates of model parameters.}
    \item{bk:}{An estimate of the \eqn{u(x, y)}{u(x, y)}.}
    \item{pb:}{The probabilities of being background event.}
    \item{opt:}{The results of optimization: the value of the log-likelihood
      function at the optimum point, its gradient at the optimum point and AIC of the model.}
    \item{rates:}{Pixel images of the estimated total intensity,
          background intensity, clustering intensity and conditional intensity.}
  }
}
\details{
  Ogata (1988) introduced the epidemic type aftershock sequence
  (ETAS) model based on Gutenberg-Richter law and modified Omori law.
  In its space-time representation (Ogata, 1998), the ETAS model is a
  temporal marked point process model, and a special case of marked
  Hawkes process, with conditional intensity function
    \deqn{
      \lambda(t, x, y | H_t) = \mu(x,y) + \sum_{t_i < t} k(m_i)g(t - t_i)f(x - x_i, y - y_i|m_i)
    }{
      lambda(t, x, y | H_t ) = mu(x, y)
       + sum[t[i] < t] k(m[i]) g(t - t[i]) f(x - x[i], y - y[i]|m[i])
    }
  where
  \describe{
      \item{\eqn{H_t}{H_t}: }{is the observational history up to time t, but not
            including t; that is
            \deqn{H_t=\{(t_i, x_i, y_i, m_i): t_i < t\}}{H_t = {(t[i], x[i], y[i], m[i]): t[i] < t}}
           }
      \item{\eqn{\mu(x,y)}{mu(x, y)}: }{is the background intensity. Currently it
            is assumed to take the semi-parametric form
            \deqn{\mu(x,y)=\mu u(x,y)}{mu(x, y) = mu u(x, y)}
            where \eqn{\mu}{mu} is an unknown constant and \eqn{u(x,y)}{u(x, y)}
            is an unknown function.
           }
      \item{\eqn{k(m)}{k(m)}: }{is the expected number of events triggered
            from an event of magnitude \eqn{m}{m} given by
           \deqn{k(m) = A\exp(\alpha(m - m_0))}{k(m[i]) =  A exp(alpha(m - m0))}
           }
      \item{\eqn{g(t)}{g(t)}: }{is the p.d.f of the occurrence
            times of the triggered events, taking the form
           \deqn{g(t) =  \frac{p-1}{c}(1 + \frac{t}{c})^{-p}}{g(t) = ((p - 1)/c)(1 + t/c)^(-p)}
           }
      \item{\eqn{f(x,y|m)}{f(x, y|m)}: }{is the p.d.f of
            the locations of the triggered events, considered to be
            either the long tail inverse power density
           \deqn{ f(x, y|m) = \frac{q-1}{\pi \sigma(m))}
             (1 + \frac{x^2 + y^2}{\sigma(m)})^{-q} }{ f(x, y|m) = (q - 1)/(pi sigma(m))
                 (1 + (x^2 + y^2)/(sigma(m)))^(-q)}
           or the light tail Gaussian density (currently not implemented)
           \deqn{ f(x,y|m)= \frac{1}{2\pi \sigma(m)}\exp(-\frac{x^2 + y^2}{2\sigma(m)}) }{f(x,y|m)= exp(-(x^2 + y^2)/(2 sigma(m)))/(2 pi sigma(m))}
           with
           \deqn{ \sigma(m) = D\exp(\gamma(m - m_0)) }{ sigma(m) = D exp( gamma (m - m0) )}
       }
  }

  The ETAS models classify seismicity into two components, background
  seismicity \eqn{\mu(x, y)}{mu(x,y)} and clustering seismicity
  \eqn{\lambda(t, x, y|H_t) - \mu(x, y)}{lambda(t, x, y|H_t) - mu(x, y)}, where
  each earthquake event, whether it is a background event or generated by
  another event, produces its own offspring according to the branching rules
  controlled by \eqn{k(m)}{k(m)}, \eqn{g(m)}{g(m)} and \eqn{f(x, y|m)}{f(x, y|m)}.

  Background seismicity rate \eqn{u(x, y)}{u(x, y)} and the model parameters
  \deqn{\theta=(\mu, A, c, \alpha, p, D, q, \gamma)}{theta = (mu, A, c, alpha, p, D, q, gamma)}
  are estimated simultaneously using an iterative approach proposed in Zhuang et al. (2002).
  First, for an initial \eqn{u_0(x, y)}{u0(x, y)}, the parameter vector
  \eqn{\theta}{theta}  is estimated by maximizing the log-likelihood function
  \deqn{l(\theta)=\sum_{i} \lambda(t_i, x_i, y_i|H_{t_i}) - \int \lambda(t, x, y|H_t) dx dy dt.}{l(theta) = sum_[i] lambda(t[i], x[i], y[i]|H_t)
         - int lambda(t, x, y|H_t) dx dy dt.}
  Then the procedure calculates the probability of being a background
  event for each event in the catalog by
  \deqn{ \phi_i = \frac{\mu(x_i, y_i)}{\lambda(t_i, x_i, y_i|H_{t_i})}. }{mu(x[i], y[i])/lambda(t[i], x[i], y[i]|H_t[i]).}
  Using these probabilities and kernel smoothing method with Gaussian kernel
  and appropriate choice of bandwidth (determined by \code{bwd} or \code{nnp}
  and \code{bwm} arguments), the background rate \eqn{u_0(x, y)}{u0(x, y)}
  is updated. These steps are repeated until the estimates converge
  (stabilize).

  The \code{no.itr} argument specifies the maximum number of iterations
  in the iterative simultaneous estimation and declustering algorithm.
  The estimates often converge in less than ten iterations. The relative
  iteration convergence tolerance and the optimization
  convergence tolerance are, respectively, determined by
  \code{rel.tol} and \code{eps} arguments.
  The progress of the computations can be traced by setting
  the \code{verbose} and \code{plot.it} arguments to be \code{TRUE}.



  If \code{cxxcode = TRUE}, then the internal function \code{etasfit}
  uses the {C++} code implemented using the \pkg{Rcpp} package,
  which allows multi-thread parallel computing on multi-core processors
  with OpenMP.
  The argument \code{nthreads} in this case determines
  the number of threads in the parallel region of the code.
  If \code{nthreads = 1} (the default case), then a serial version of the
  {C++} code carries out the computations.


  This version of the ETAS model assumes that the earthquake catalog
  is complete and the data are stationary in time. If the catalog
  is incomplete or there is non-stationarity (e.g. increasing or cyclic
  trend)  in the time of events, then the results of this function are
  not reliable.
}
\note{
  This function is based on a \code{C} port of the original
  \code{Fortran} code by Jiancang Zhuang, Yosihiko Ogata and
  their colleagues. The \code{etas} function is intended to be
  used for small and medium-size earthquake catalogs.
  For large earthquake catalogs, due to time-consuming
  computations, it is highly recommended to
  use the parallel \code{Fortran} code on a server machine.
  The \code{Fortran} code (implemented for
  parallel/non-parallel computing) can be obtained from
  \url{http://bemlar.ism.ac.jp/zhuang/software.html}.
}
\references{
  Ogata Y (1988).
  Statistical Models for Earthquake Occurrences and Residual Analysis for
  Point Processes.
  \emph{Journal of the American Statistical Association},
  \bold{83}(401), 9--27. \href{http://dx.doi.org/10.2307/2288914}{doi:10.2307/2288914}.

  Ogata Y (1998).
  Space-time Point-process Models for Earthquake Occurrences.
  \emph{Annals of the Institute of Statistical Mathematics},
  \bold{50}(2), 379--402.
  \href{http://dx.doi.org/10.1023/a:1003403601725}{doi:10.1023/a:1003403601725}.

  Zhuang J, Ogata Y, Vere-Jones D (2002).
  Stochastic Declustering of Space-Time Earthquake Occurrences.
  \emph{Journal of the American Statistical Association},
  \bold{97}(458), 369--380.
  \href{http://dx.doi.org/10.1198/016214502760046925}{doi:10.1198/016214502760046925}.

  Zhuang J, Ogata Y, Vere-Jones D (2006).
  Diagnostic Analysis of Space-Time Branching Processes for Earthquakes.
  In \emph{Case Studies in Spatial Point Process Modeling},
  pp. 275--292. Springer Nature.
  \href{http://dx.doi.org/10.1007/0-387-31144-0_15}{doi:10.1007/0-387-31144-0_15}.

  Zhuang J (2011).
  Next-day Earthquake Forecasts for the Japan Region Generated by
  the ETAS Model.
  \emph{Earth, Planets and Space},
  \bold{63}(3), 207--216.
  \href{http://dx.doi.org/10.5047/eps.2010.12.010}{doi:10.5047/eps.2010.12.010}.

}
\seealso{

  \code{\link{catalog}} for constructing data.
  \code{\link{probs}} for estimated declustering probabilities.
  \code{\link{resid.etas}} for diagnostic plots.
}
\examples{
  # fitting the ETAS model to an Iranian catalog
  # preparing the catalog
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
  iran.fit <- etas(iran.cat, param0=param0)}


  # fitting the ETAS model to an Italian catalog
  # preparing the catalog
  italy.cat <- catalog(italy.quakes, dist.unit="km")
  \dontrun{
  plot(italy.cat)}

  # setting initial parameter values
  mu <- 1
  k0 <- 0.005
  c <- 0.005
  alpha <- 1.05
  p <- 1.01
  D <- 1.1
  q <- 1.52
  gamma <- 0.6
  # reparametrization: transform k0 to A
  A <- pi * k0 / ((p - 1) * c^(p - 1) * (q - 1) * D^(q - 1))
  param0 <- c(mu, A, c, alpha, p, D, q, gamma)

  # fitting the model
  \dontrun{
  nthreads <- parallel::detectCores()
  italy.fit <- etas(italy.cat, param0, nthreads=nthreads)}


  # fitting the ETAS model to a Japanese catalog
  # setting the target polygonal study region
  jpoly <- list(long=c(134.0, 137.9, 143.1, 144.9, 147.8,
      137.8, 137.4, 135.1, 130.6), lat=c(31.9, 33.0, 33.2,
      35.2, 41.3, 44.2, 40.2, 38.0, 35.4))
  # preparing the catalog
  japan.cat <- catalog(japan.quakes, study.start="1953-05-26",
      study.end="1990-01-08", region.poly=jpoly, mag.threshold=4.5)
  \dontrun{
  plot(japan.cat)}

  # setting initial parameter values
  param0 <- c(0.592844590, 0.204288231, 0.022692883, 1.495169224,
  1.109752319, 0.001175925, 1.860044210, 1.041549634)

  # fitting the model
  \dontrun{
  nthreads <- parallel::detectCores()
  japan.fit <- etas(japan.cat, param0, nthreads=nthreads)}
}

\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}

