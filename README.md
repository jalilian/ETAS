<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ETAS)](https://CRAN.R-project.org/package=ETAS)
[![CRAN_Download_Count](http://cranlogs.r-pkg.org/badges/ETAS)](https://CRAN.R-project.org/package=ETAS)
[![R Build Status](https://github.com/jalilian/ETAS/workflows/R-CMD-check/badge.svg)](https://github.com/jalilian/ETAS/actions)
[![Dependencies](https://tinyverse.netlify.com/badge/ETAS)](https://cran.r-project.org/package=ETAS)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](http://www.gnu.org/licenses/gpl-2.0.html)
[![Last Commit](https://img.shields.io/github/last-commit/jalilian/ETAS)](https://github.com/jalilian/ETAS)
<!-- badges: end -->

# ETAS model

An earthquake catalog is a chronologically ordered list of time, coordinates of epicenter, magnitude and focal depth of all recorded earthquakes with magnitudes greater than or equal to a certain threshold that occurred inside or in the vicinity of a geographical region during a specified time period. Among different proposed models, the epidemic type aftershock sequence (ETAS) model is the most widely used statistical model to describe the underlying process that generates an earthquake catalogs. 

The space-time version of the ETAS model is a spatio-temporal marked point process model. It is a semi-parametric model that describes the background and triggering seismic activities in a geographical region and can be used for earthquake declustering. However, estimation of the ETAS model parameters is computationally challenging. The 'ETAS' package fits the ETAS model to an earthquake catalog. The `etas` function in the package is based on a C port of a [Fortran program](http://bemlar.ism.ac.jp/zhuang/software.html) by [Jiancang Zhuang](http://bemlar.ism.ac.jp/zhuang/), [Yosihiko Ogata](https://www.ism.ac.jp/~ogata/) and their colleagues.

## Installation

To install the package from [CRAN](https://CRAN.R-project.org/package=ETAS), run the following in R:
```R
install.packages('ETAS')
```

You can also install the current version of the package on GitHub by running:
```R
require(remotes)
install_github('jalilian/ETAS')
```

If [remotes](https://github.com/r-lib/remotes) is not installed, you should first run:

```R
install.packages('remotes')
```

## Vignette

Jalilian, A. (2019). [ETAS: An R Package for Fitting the Space-Time ETAS Model to Earthquake Data](https://www.jstatsoft.org/htaccess.php?volume=088&type=c&issue=01&filename=paper). Journal of Statistical Software, 88(CN 1), 1-39. [DOI:10.18637/jss.v088.c01](http://dx.doi.org/10.18637/jss.v088.c01)
 
## Parallel computing

Computations of the conditional intensity function, the log-likelihood function, declustering probabilities and the Davidon-Fletcher-Powell algorithm for optimization are all written in C code. As of version 0.3, a new C++ code is implemented using the [Rcpp](https://www.rcpp.org/) package which allows multi-thread parallel computing on multi-core processors with OpenMP and [suported platforms](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#OpenMP-support). The argument `nthreads` in `etas` function determines the number of threads to be used in the parallel region of the code. If `nthreads = 1` (the default), then a serial version of the C++ code carries out the computations. The `detectCores` function in [parallel](http://stat.ethz.ch/R-manual/R-devel/library/parallel/html/parallel-package.html) package can be consulted to find out the overall number of available threads on a given machine:
```R
parallel::detectCores()
```
Parallel computing (`nthreads > 1`) reduces the computation time for large earthquake catalogs. However, resource usage and limitations should be considered when setting `nthreads`.
