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
lstirling <- function(n, m) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
      !is.numeric(m) || length(m) != 1 || !is.finite(m)) {
    stop("n and m must be single finite numeric values")
  }

  if (n > 1000000 || m > 1000000) {
    stop("n and m must be <= 1,000,000 to prevent resource exhaustion")
  }

  # special cases first
  if (m < 0 || n < 0 || n < m) {
    # we really need everything to be positive to get a value
    return(-Inf)
  } else if (m == 0 && n == 0) {
    return(0)
  } else if (m == 0 && n > 0) {
    # the stirling number is 0 if m = 0 and n > 0
    return(-Inf)
  } else if (m == 1 && n >= 1) {
    # this ends up being (n-1)!, approx with log(x!) ~= x(log x) - x + 1
    # R has something built in - just use lgamma(n)
    return(base::lgamma(n))
  } else if (n == m) {
    return(0)
  } else if (n > 1 && m > 1 && n > m) {
    phiprime <- function(x) {
      base::digamma(x + n + 1) - base::digamma(x + 1) - m / x
    }
    n <- n - 1
    m <- m - 1
    x0_res <- base::tryCatch(
      {
        stats::uniroot(phiprime, c(0.1, n * m), tol = .Machine$double.eps, check.conv = TRUE)$root
      },
      warning = function(w) {
        return(NULL)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (base::isTRUE(base::is.null(x0_res))) {
      # Fallback or error if uniroot fails
      base::warning("uniroot failed to converge in lstirling")
      return(NaN)
    }
    x0 <- x0_res
    t0 <- m / (n - m)
    B <- base::lgamma(x0 + n + 1) - base::lgamma(x0 + 1) - m * base::log(x0) - n * base::log(t0 + 1) + m * base::log(t0)
    gt <- 1 / x0 * base::sqrt(m * (n - m) / n / (base::trigamma(x0 + n + 1) - base::trigamma(x0 + 1) + m / x0 / x0))
    return(B + base::log(gt) + base::lchoose(n, m))
  }
}
