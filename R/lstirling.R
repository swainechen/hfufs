#
# Swaine Chen; swainechen@gmail.com
# https://github.com/swainechen/hfufs
# MIT License
#
lstirling <- function(n,m) {
  # use approximation from http://oai.cwi.nl/oai/asset/2304/2304A.pdf
  # this is Temme 1993 Studies in Applied Mathematics 89:233-243, equation (3.5)
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
