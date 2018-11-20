#' Calculate Strobeck's S for arbitrary parameters
#'
#' Returns Strobeck's S statistic.
#'
#' Strobeck's S is a population genetics statistic that is used to detect
#' population structure. It is related conceptually (and computationally) to
#' Fu's Fs, which has been argued to be a more appropriate statistic due to
#' it's higher dynamic range, owing to the use of a logit transformation for
#' Fu's Fs. In cases where Fu's Fs becomes more negative, Strobeck's S
#' approaches 1.
#'
#' Computationally, Strobeck's S and Fu's Fs can be defined in terms of Stirling
#' numbers of the first kind. These numbers can grow large very quickly.
#'
#' `hstrobecks` calculates Strobeck's S directly using Stirling numbers,
#' switching to a logarithmic esimator for Stirling numbers (implemented in the
#' lstirling function) when needed to avoid overflow errors.
#'
#' @param n The number of total sequences/individuals
#' @param k The number of unique alleles
#' @param theta The average pairwise nucleotide divergence
#' @return Strobeck's S
#' @seealso `hfufs` for a Fu's Fs calculator and `lstirling` for the Stirling
#' number estimator
#' @export
#' @references
#' Fu, Y.X. (1996) New statistical tests of neutrality for DNA samples from
#' a population. Genetics 143:557-70
#' Strobeck, C. (1987) Average number of nucleotide differences in a sample
#' from a single subpopulation: a test for population subdivision. Genetics
#' 117:149-53
#'
#' @examples
#' n <- 100
#' k <- 30
#' theta <- 12.345
#' hstrobecks(n, k, theta)
#' # 0.7576966
#'
hstrobecks <- function(n, k, theta) {
  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0, then Strobeck's S is 0
  # if k == 1, then Strobeck's S is 0 if theta > 0, 1 if theta == 0
  # Strobeck's S should be 1 if theta is 0 and k >= 1
  if(k <= 1) {
    if (k == 1 & theta == 0) {
      return(1)
    } else {
      return(0)
    }
  }
  if(theta == 0) {
    return(1)
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
