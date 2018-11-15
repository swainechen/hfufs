#' Calculate Fu's Fs for arbitrary parameters
#'
#' Returns Fu's Fs statistic.
#'
#' Fu's Fs is a population genetics statistic that is useful for detecting
#' loci that are responsible for population expansion, for example. Fu's Fs
#' can be formulated as a calculation that involves Stirling numbers of the
#' first kind. These can get large very quickly and exceed the floating point
#' range for modern genomic data sets.
#'
#' `hfufs` calculates Fu's Fs either directly, if the numbers are small
#' enough, or using a logarithmic estimator for Stirling numbers implemented
#' in the lstirling function.
#'
#' @param n The number of total sequences/individuals
#' @param k The number of unique alleles
#' @param theta The average pairwise nucleotide divergence
#' @return Fu's Fs
#' @seealso `hstrobecks` for a Strobeck's S calculator and `lstirling` for
#' the Stirling number estimator
#' @export
#' @references
#' Fu, Y.X. (1997) Statistical Tests of Neutrality of Mutations Against
#' Population Growth, Hitchhiking and Background Selection. Genetics
#' 147:915-25
#'
#' @examples
#' n <- 100
#' k <- 30
#' theta <- 12.345
#' hfufs(n, k, theta)
#' # -0.7374915
#'
hfufs <- function(n, k, theta) {
  # this value is pretty arbitrary, but it leaves things accurate when
  # calculate logit (Sp) (i.e. Fu's Fs).
  # Another way would be to split the calculation based on whether Fu's Fs is
  # positive or negative
  too_small <- 0.1

  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0 or 1, then Fu's Sp is 1, logit (Sp) is infinity
  # if k > 1 and theta == 0, then Sp is 0, logit (Sp) is -infinity
  if(k <= 1) {
    return(Inf)
  }
  if(k > 1 & theta == 0) {
    return(-Inf)
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
      for(i in k:n) {	# this is k:n for Fu's Fs, 1:k for Strobeck's S
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
