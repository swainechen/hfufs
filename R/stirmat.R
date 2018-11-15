#
# Swaine Chen; swainechen@gmail.com
# https://github.com/swainechen/hfufs
# MIT License
#
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
