#' Calculate a matrix of Stirling numbers of the first kind
#'
#' Returns a matrix of Stirling numbers of the first kind
#'
#' Stirling numbers of the first (and second) kind are useful numbers for
#' combinatorics. In the context of sequence analysis, Stirling numbers of the
#' first kind have entered into the calculation of statistics that measure
#' expected numbers of alleles in a population.
#'
#' Stirling numbers are characterized by two parameters, `n` and `m` (or `n`
#' and `k` in some literature). `stirmat` returns an n by m matrix of all the
#' respective Stirling numbers using a recursion relationship.
#'
#' `stirmat` is a direct port of the code from James Cai's PGEToolbox, which
#' was originally written for Matlab.
#'
#' @param n The first parameter for Stirling numbers
#' @param m The second parameter for Stirling numbers
#' @return A matrix with n rows and m columns, the elements of which are
#' Stirling numbers of the first kind, for all combinations of n and m up to
#' the supplied n and m
#' and column, respectively
#' @seealso `lstirling` for a logarithmic estimator of large Stirling numbers
#' @export
#' @references
#' Cai, J.J. (2008) PGEToolsbox: A Matlab toolbox for population genetics and
#' evolution. J Hered 99(4):438-40
#'
#' @examples
#' n <- 20
#' m <- 10
#' stirling_matrix <- stirmat(n, m)
#' stirling_matrix[n,m]
#' # 3.819221e+14
#'
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
