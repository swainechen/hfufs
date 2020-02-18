#' Calculate Fu's Fs for arbitrary parameters - asymptotic approximation
#'
#' Returns Fu's Fs statistic.
#'
#' Fu's Fs is a population genetics statistic that is useful for detecting
#' loci that are responsible for population expansion, for example. Fu's Fs
#' can be formulated as a calculation that involves Stirling numbers of the
#' first kind. These can get large very quickly and exceed the floating point
#' range for modern genomic data sets.
#'
#' `afufs` calculates Fu's Fs using a single term asymptotic approximation
#' This contrasts with `hfufs` which uses a Stirling number approximator for
#' each term. `afufs` is therefore simultaneously more accurate and faster.
#'
#' @param n The number of total sequences/individuals
#' @param k The number of unique alleles
#' @param theta The average pairwise nucleotide divergence
#' @return Fu's Fs
#' @export
#' @references
#' Chen, S.L. and Temme, N. (2020) Computing sums of Stirling numbers of the
#' first kind: Application to population genetics statistics. In preparation.
#'
#' @examples
#' n <- 100
#' k <- 30
#' theta <- 12.345
#' afufs(n, k, theta)
#' # -0.7368616
#'
afufs <- function(n, k, theta) {
  # Strobeck's S is prob of k alleles or fewer, Fu's Sp is k alleles or greater
  # if k == 0 or 1, then Fu's Sp is 1, logit (Sp) is infinity
  # if k > 1 and theta == 0, then Sp is 0, logit (Sp) is -infinity
  if(k <= 1) {
    return(Inf)
  }
  if(k > 1 & theta == 0) {
    return(-Inf)
  }
  if(k==n) {
    return(hfufs(n, k, theta))
  }

  # we derive this as S'_(n+1),(k+1)(theta) - so decrement n, k first
  n <- n-1
  k <- k-1
  phi <- function(x) { lgamma(x+n+1) - lgamma(x+1) - k*log(x) }
  phiprime <- function(x) { digamma(x+n+1) - digamma(x+1) - k/x }
  chi <- function(t) { n*log(1+t) - k*log(t) }
  chiprimeprime <- function(t) { -n/(1+t)/(1+t) + k/t/t }
  z0 <- tryCatch (
    { uniroot(phiprime, c(0.1, n*k), tol=.Machine$double.eps, check.conv=TRUE)$root },
    warning = function(w) { return(theta) },
    error = function(e) { return(theta) }
  )
  if (z0 == theta) { return(hfufs(n+1, k+1, theta)) }
  t0 <- k/(n-k)
  B <- phi(z0) - chi(t0)
  chitau <- phi(theta) - B
  revchi <- function(x) { chi(x) - chitau }
  if (z0 > theta) {
    # a check to make sure delta is small enough - chi(t) should approach +inf as t approaches 0 from the positive side
    delta <- 0.00001
    while(revchi(delta) < 0) { delta <- delta/10 }
    tau <- uniroot(revchi, c(delta, t0), tol=.Machine$double.eps, check.conv=TRUE)$root
  } else {
    # otherwise make sure the range is big enough on the right
    tempfactor <- 10
    while (revchi(tempfactor*t0) < 0) { tempfactor <- tempfactor*10 }
    tau <- uniroot(revchi, c(t0, tempfactor*t0), tol=.Machine$double.eps, check.conv=TRUE)$root
  }
  f_t0 <- 1/(z0 - theta) * sqrt(chiprimeprime(t0) / (trigamma(z0+n+1) - trigamma(z0+1) + k/z0/z0))
  G0 <- f_t0 - 1/(t0 - tau)
  temp <- exp(lchoose(n, k-1) - chitau)
  Sprime <- pbeta(tau/(1+tau), k, n-k+1) + temp * G0
  Tprime <- pbeta(1/(1+tau), n-k+1, k) - temp * G0
  if (Sprime < 0.5) {
    return(log(Sprime) - log(1-Sprime))
  } else {
    return(log(1-Tprime) - log(Tprime))
  }
}
