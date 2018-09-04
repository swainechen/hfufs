#
# Swaine Chen; swainechen@gmail.com
# https://github.com/swainechen/hfufs
# MIT License
#
lstirling <- function(n,m) {
  # use approximation from http://oai.cwi.nl/oai/asset/2304/2304A.pdf
  # this is equation (3.5)
  # special cases first
  if (m < 0 | n < 0 | n < m) {
    # we really need everything to be positive to get a value
    return(NULL)
  } else if (m == 0 & n == 0) {
    return(0)
  } else if (m == 0 & n > 0) {
    # the stirling number is 0 if m = 0 and n > 0
    return(NULL)
  } else if (m==1 & n >= 1) {
    # this ends up being (n-1)!, approx with log(x!) ~= x(log x) - x + 1
    # R has something built in - just use lgamma(n)
    return(lgamma(n))
  } else if (n == m) {
    return(0)
  } else if (n > 1 & m > 1 & n > m) {
    phiprime <- function(x) { digamma(x+n+1) - digamma(x+1) - m/x }
    n <- n-1
    m <- m-1
    x0 <- uniroot(phiprime, c(0.1,n*m))$root
    t0 <- m/(n-m)
    B <- lgamma(x0+n+1) - lgamma(x0+1) - m* log(x0) - n*log(t0+1) + m*log(t0)
    gt <- 1/x0 * sqrt(m*(n-m)/n/(trigamma(x0+n+1) - trigamma(x0+1) + m/x0/x0))
    return(B + log(gt) + lchoose(n,m))
  }
}

stirmat <- function(n, m) {
  # from James Cai PGEToolbox
  if (n <= 0 | m <= 0) {
    return(NULL)
  }
  s_n <- matrix(nrow=n, ncol=m)
  s_n[1,1] <- 1
  for(j in 2:m) {
    s_n[1,j] <- 0
  }
  for(i in 2:n) {
    s_n[i,1] <- -(i-1) * s_n[i-1,1]
    for(j in 2:m) {
      s_n[i,j] <- s_n[i-1,j-1] - (i-1)*s_n[i-1,j]
    }
  }
  return(s_n)
}

lstrobecks <- function(n, k, theta) {
  # use log approximations to calculate Strobeck's S
  # this is the "complement" of Fu's Fs - use indices below k instead of above
  lSn <- sum(log(theta + (0:(n-1))))
  S <- 0
  for (i in 1:k) {
    S <- S + exp(lstirling(n,i) + i*log(theta) - lSn)
  }
  return(S)
}

hfufs <- function(n, k, theta) {
  too_small <- 0.1
  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0, then Strobeck's S is 0
  # if k == 1, then Strobeck's S is 0 if theta > 0, 1 if theta == 0
  # if k == 0 or 1, then Fu's Sp is 1, logit (F_s) is infinity
  # Strobeck's S should be 1 if theta is 0, 0 if theta > 0
  if(k <= 1) {
    Fu_Fs <- Inf
    if (k == 1 & theta == 0) {
      Strobeck_S <- 1
    } else {
      Strobeck_S <- 0
    }
    return(Fu_Fs)
  }
  # if n is small enough, then just calculate directly
  if(n <= 30) {
    s_n <- stirmat(n,n)
    Sn <- 1;
    for(i in 0:(n-1)) {
      Sn <- Sn * (theta + i)
    }
    if(!is.infinite(Sn)) {
      Sp <- 0
      for(i in k:n) {
        Sp <- Sp + abs(s_n[n,i]) * theta**i
      }
      Sp <- Sp / Sn
      if ((1-Sp) < too_small) {
        S <- 0
        for(i in 1:(k-1)) {
          S <- S + abs(s_n[n,i]) * theta**i
        }
        S <- S / Sn
        if(!is.nan(S) & S > 0 & 1-S > 0) {
          return(log(1-S) - log(S))
        }
      } else { 
        if (!is.nan(Sp) & Sp > 0 & 1-Sp > 0) {
          return(log(Sp) - log(1-Sp))
        }
      }
    }
  }
  # use log approximations to calculate Fu's Fs
  # this is a fallback in case previous calculation hit infinity
  lSn <- sum(log(theta + (0:(n-1))))
  Sp <- 0
  for (i in k:n) {
    Sp <- Sp + exp(lstirling(n,i) + i*log(theta) - lSn)
  }
  if ((1-Sp) < too_small) {
    S <- 0
    for (i in 1:(k-1)) {
      S <- S + exp(lstirling(n,i) + i*log(theta) - lSn)
    }
    return(log(1-S) - log(S))
  } else {
    return(log(Sp) - log(1-Sp))
  }
}
