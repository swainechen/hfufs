#' Calculate the logarithm of a Stirling number of the first kind
#'
#' Returns an approximation of the logarithm of a Stirling number of the first
#' kind. This is useful because for larger parameter values, these Stirling
#' numbers can overflow the floating point range.
#'
#' @param n The first parameter (subscript in the notation used in Temme 1993)
#' @param m The second parameter (superscript in parentheses in the notation
#' used in Temme 1993)
#' @return The logarithm of the Stirling number of the first kind.
#' @seealso `stirmat` for a direct Stirling number calculator
#' @export
#' @references
#' Temme, N.M. (1993) Studies in Applied Mathematics 89:233-243 (Equation 3.5)
#'
#' @examples
#' n <- 10
#' m <- 5
#' exp(lstirling(n, m))
#' # 267854.5
#'
lstirling <- function(n,m) {
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
