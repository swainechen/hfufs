#' Calculate population genetics statistics from a PopGenome GENOME object
#'
#' Takes in a PopGenome GENOME object (from PopGenome::readData) and
#' calculates basic population genetics statistics. Uses afufs if needed.
#' Returns a data frame with most of the statistics. The original GENOME
#' object is modified as normal by the PopGenome functions.
#'
#' Runs these functions from PopGenome:
#'   diversity.stats
#'   neutrality.stats
#'
#' @export
#'
#' @examples
#' \dontrun {
#' fasta_file <- "name_of_aligned.fasta"
#' pg.object <- hf.readData(fasta_file)
#' pg.dataframe <- hf.alignment.stats(pg.object)
#' }
#'
hf.alignment.stats <- function(go, slide=F, window=1000, step=500) {
  # takes in a GENOME object as from PopGenome readData
  if(summary(go)[2]=="GENOME") {
    numindividuals <- length(get.individuals(go)[[1]])
    single_value <- 1/numindividuals
    if (slide) {
      slide_go <- sliding.window.transform(go, width=window, jump=step, type=2, whole.data=T)
      slide_go <- diversity.stats(slide_go, pi=T)
      slide_go <- neutrality.stats(slide_go, detail=T, do.R2=T)
      n <- data.frame(get.neutrality(slide_go)[[1]])
      n$x.start <- as.numeric(sapply(strsplit(slide_go@region.names, split=" - "), "[")[1,])
      n$x.end <- as.numeric(sapply(strsplit(sub(" :", "", slide_go@region.names), split=" - "), "[")[2,])
      n$numindividuals <- numindividuals
      n$pi <- get.diversity(slide_go)[[1]][,3]
      h <- slide_go@region.stats@haplotype.counts
      n$n.sequences <- as.numeric(lapply(h, sum))
      templist <- lapply(h, ncol)
      templist[sapply(templist, is.null)] <- NA
      n$n.haplotypes <- as.numeric(templist)
      n$n.singleton.haplotypes <- as.numeric(lapply(h, function(x) length(which(x==1))))
      for(i in which(is.nan(n$Fu.F_S))) {
        n$Fu.F_S[i] <- afufs(n$n.sequences[i], n$n.haplotypes[i], n$pi[i])
      }
    } else {
      go <- diversity.stats(go, pi=T)
      go <- neutrality.stats(go, detail=T, do.R2=T)
      n <- data.frame(get.neutrality(go)[[1]])
      n$numindividuals <- numindividuals
      n$pi <- get.diversity(go)[[1]][,3]
      h <- go@region.stats@haplotype.counts[[1]]
      n$n.sequences <- sum(h)
      n$n.haplotypes <- ncol(h)
      n$n.singleton.haplotypes <- length(which(h == 1))
      n$n.consensus.haplotypes <- max(h)
      if(is.nan(n$Fu.F_S)) n$Fu.F_S <- afufs(n$n.sequences, n$n.haplotypes, n$pi)
    }
    return(n)
  }
}
