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
  if (base::isTRUE(!base::is.numeric(n) || base::length(n) != 1 || !base::is.finite(n) || n <= 0 || n %% 1 != 0 ||
      !base::is.numeric(k) || base::length(k) != 1 || !base::is.finite(k) || k < 0 || k %% 1 != 0 ||
      !base::is.numeric(theta) || base::length(theta) != 1 || !base::is.finite(theta) || theta < 0 || theta > 2000000000 ||
      k > n)) {
    base::stop("n, k, and theta must be single finite numeric values; n and k must be whole numbers; n > 0, k >= 0, 0 <= theta <= 2,000,000,000, and k <= n")
  }

  if (base::isTRUE(n > 1000000 || k > 1000000)) {
    base::stop("n and k must be <= 1,000,000 to prevent resource exhaustion")
  }

  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0, then Strobeck's S is 0
  # if k == 1, then Strobeck's S is 0 if theta > 0, 1 if theta == 0
  # Strobeck's S should be 1 if theta is 0 and k >= 1
  if (base::isTRUE(k <= 1)) {
    if (base::isTRUE(k == 1 && theta == 0)) {
      return(1)
    } else {
      return(0)
    }
  }
  if (base::isTRUE(theta == 0)) {
    return(1)
  }
  # if n is small enough, then just calculate directly
  if (base::isTRUE(n <= 30)) {
    s_n <- stirmat(n, n)
    Sn <- base::exp(base::lgamma(theta + n) - base::lgamma(theta))
    if (base::isTRUE(!base::is.infinite(Sn))) {
      # Choose the shorter range to minimize iterations (DoS protection).
      if (base::isTRUE(k < (n / 2) + 1)) {
        indices <- 1:k
        vals <- base::vapply(indices, function(i) {
          base::abs(s_n[n, i]) * theta**i
        }, FUN.VALUE = 0.0)
        Ss <- base::sum(base::sort(vals)) / Sn
      } else if (base::isTRUE(k < n)) {
        indices <- (k + 1):n
        vals <- base::vapply(indices, function(i) {
          base::abs(s_n[n, i]) * theta**i
        }, FUN.VALUE = 0.0)
        Ss <- 1 - (base::sum(base::sort(vals)) / Sn)
      } else {
        # k == n, so Ss must be 1. This also avoids (n+1):n crash.
        Ss <- 1
      }
      if (base::isTRUE(!base::is.nan(Ss) && !base::is.infinite(Ss))) {
        return(Ss)
      }
    }
  }

  # use log approximations to calculate Strobeck's S
  # this is a fallback in case previous calculation overflows or for larger n
  # We use lgamma to avoid large vector allocation for 0:(n-1)
  lSn <- base::lgamma(theta + n) - base::lgamma(theta)

  # Choose the shorter range to minimize iterations and calls to lstirling (DoS protection).
  if (base::isTRUE(k < (n / 2) + 1)) {
    indices <- 1:k
    vals <- base::vapply(indices, function(i) {
      base::exp(lstirling(n, i) + i * base::log(theta) - lSn)
    }, FUN.VALUE = 0.0)
    Ss <- base::sum(base::sort(vals))
  } else if (base::isTRUE(k < n)) {
    indices <- (k + 1):n
    vals <- base::vapply(indices, function(i) {
      base::exp(lstirling(n, i) + i * base::log(theta) - lSn)
    }, FUN.VALUE = 0.0)
    Ss <- 1 - base::sum(base::sort(vals))
  } else {
    # k == n
    Ss <- 1
  }

  return(Ss)
}
