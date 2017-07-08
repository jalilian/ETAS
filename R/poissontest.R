# obtained and adapted from 
# http://statistics.berkeley.edu/~stark/Code/Quake/permutest.r
# based on
#  Luen, B., & Stark, P. B. (2012). 
#   Poisson tests of declustered catalogues. 
#   Geophysical journal international, 189(1), 691-700.

poisson.test <- function(object, n.perm=1000, verbose=TRUE)
{
  # Extract ranks (assume no ties)
  long <- object$longlat.coord$long
  lat <- object$longlat.coord$lat
  x.rank <- rank(long)
  y.rank <- rank(lat)
  
  # number of events
  n <- length(x.rank)
  
  # Find empirical distribution of spatial ranks
  xy.upper <- t(outer(y.rank, y.rank, "<=")) %*% outer(x.rank, x.rank, "<=")
  # xy.upper[I,J] is the number of points
  # with y <= y[i], x <= x[j]
  # y is row, x is column
  
  ### Distance function
  
  distfind <- function(x.rank, y.rank, xy.upper)
  {
    dfv <- 0
    xyz.temp <- matrix(0, n, n)
    # xyz.temp is the number of points
    # with y <= y[i], x <= x[j], z <= Z
    # i.e. empirical distribution at time Z
    # Now go through search space chronologically
    # update xyz.temp
    # find the max; check the min isn't close; if it is, look around
    #uu <- lapply(1:n, function(i)
    #{ (y.rank >= y.rank[i]) %*% t(x.rank >= x.rank[i]) })
    
    for(Z in 1:n)
    {
      xyz.temp <- xyz.temp + (y.rank >= y.rank[Z]) %*% t(x.rank >= x.rank[Z])
      dist.matrix <- xyz.temp/n - xy.upper/n * Z/n
      dfv <- max(dfv, abs(dist.matrix))
    }
    return(dfv)
  }
  
  teststat <- distfind(x.rank, y.rank, xy.upper)
  
  permustat <- rep(NA, n.perm)
  
  # permuting the occurrence times of events 
  
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
  
  # P-value
  p.value <- mean(permustat >= teststat)
  return(list(statistic=teststat, p.value=p.value))
}
