#' Calculate population genetics statistics from a PopGenome GENOME object
#'
#' Takes in a PopGenome GENOME object (from PopGenome::readData) and
#' calculates basic population genetics statistics. Uses afufs if needed.
#'
#' @export
#'
#'
#'
hf.alignment.stats <- function(x, slide=F, width=1000, step=500) {
  # takes in a GENOME object as from PopGenome readData
  numindividuals <- length(get.individuals(d)[[1]])
  single_value <- 1/numindividuals
  if(summary(x)[2]=="GENOME") {
    d <- x
    if (slide) {
      window_width <- 1000
      window_step <- 500
      slide_d <- sliding.window.transform(d, width=window_width, jump=window_step, type=2, whole.data=T)
      slide_d <- diversity.stats(slide_d, pi=T)
      slide_d <- neutrality.stats(slide_d, detail=T, do.R2=T)
      n <- data.frame(get.neutrality(slide_d)[[1]])
      n$x.start <- as.numeric(sapply(strsplit(slide_d@region.names, split=" - "), "[")[1,])
      n$x.end <- as.numeric(sapply(strsplit(sub(" :", "", slide_d@region.names), split=" - "), "[")[2,])
      n$numindividuals <- numindividuals
      n$pi <- get.diversity(slide_d)[[1]][,3]
      h <- slide_d@region.stats@haplotype.counts
      n$n.sequences <- as.numeric(lapply(h, sum))
      templist <- lapply(h, ncol)
      templist[sapply(templist, is.null)] <- NA
      n$n.haplotypes <- as.numeric(templist)
      n$n.singleton.haplotypes <- as.numeric(lapply(h, function(x) length(which(x==1))))
      for(i in which(is.nan(n$Fu.F_S))) {
        n$Fu.F_S[i] <- afufs(n$n.sequences[i], n$n.haplotypes[i], n$pi[i])
      }
    } else {
      d <- diversity.stats(d, pi=T)
      d <- neutrality.stats(d, detail=T, do.R2=T)
      n <- data.frame(get.neutrality(d)[[1]])
      n$numindividuals <- numindividuals
      n$pi <- get.diversity(d)[[1]][,3]
      h <- d@region.stats@haplotype.counts[[1]]
      n$n.sequences <- sum(h)
      n$n.haplotypes <- ncol(h)
      n$n.singleton.haplotypes <- length(which(h == 1))
      n$n.consensus.haplotypes <- max(h)
      if(is.nan(n$Fu.F_S)) n$Fu.F_S <- afufs(n$n.sequences, n$n.haplotypes, n$pi)
    }
    return(n)
  }
}
