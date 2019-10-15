
# number of decimal places
decimalplaces <- function(x) 
{
  out <- integer(length(x))
  idx <- which((x %% 1) != 0)
  for (i in idx)
  {
    tmp <-  strsplit(sub('0+$', '', as.character(x[i])), ".", fixed=TRUE)
    out[i] <- nchar(tmp[[1]][[2]])
  }
  return(out)
}

# 
roundoffErr <- function(x)
{
  x + runif(length(x), -1/2, 1/2) * 10^(-decimalplaces(x))
}
