#
# Swaine Chen; swainechen@gmail.com
# https://github.com/swainechen/hfufs
# MIT License
#
hstrobecks <- function(n, k, theta) {
  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0, then Strobeck's S is 0
  # if k == 1, then Strobeck's S is 0 if theta > 0, 1 if theta == 0
  # Strobeck's S should be 1 if theta is 0, 0 if k == 1 and theta > 0
  if(k <= 1) {
    if (k == 1 & theta == 0) {
      return(1)
    } else {
      return(0)
    }
  }
  # if n is small enough, then just calculate directly
  if(n <= 30) {
    s_n <- stirmat(n,n)
    Sn <- 1;
    for(i in 0:(n-1)) {
      Sn <- Sn * (theta + i)
    }
    if(!is.infinite(Sn)) {
      Ss <- 0
      for(i in 1:k) {	# this is 1:k for Strobeck's S, k:n for Fu's Fs
        Ss <- Ss + abs(s_n[n,i]) * theta**i
      }
      Ss <- Ss / Sn
      if(!is.nan(Ss) & !is.infinite(Ss)) {
        return(Ss)
      }
    }
  }
  # use log approximations to calculate Strobeck's S
  # this is a fallback in case previous calculation overflows
  # make sure to add values in ascending order for accuracy
  # i.e. explicitly don't just use sum(vals)
  lSn <- sum(log(theta + (0:(n-1))))
  vals <- rep(0, k)
  for (i in 1:k) {
    vals[i] <- exp(lstirling(n,i) + i*log(theta) - lSn)
  }
  vals <- sort(vals)
  Ss <- 0
  for (i in 1:k) {
    Ss <- Ss + vals[i]
  }
  return(Ss)
}
