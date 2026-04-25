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
#' \dontrun{
#' fasta_file <- "name_of_aligned.fasta"
#' pg.object <- hf.readData(fasta_file)
#' pg.dataframe <- hf.alignment.stats(pg.object)
#' }
#'
hf.alignment.stats <- function(go, slide=FALSE, window=1000, step=500) {
  # Dependency check
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("The 'PopGenome' package is required but not installed.")
  }

  # Input validation
  if (missing(go)) {
    stop("go must be provided")
  }

  # Check if go is a GENOME object safely
  is_genome <- tryCatch({
    summ <- summary(go)
    is.character(summ) && length(summ) >= 2 && isTRUE(summ[2] == "GENOME")
  }, error = function(e) FALSE)

  if (!is_genome) {
    stop("go must be a PopGenome GENOME object")
  }

  if (!is.logical(slide) || length(slide) != 1 || is.na(slide)) {
    stop("slide must be a single non-NA logical value")
  }

  if (!is.numeric(window) || length(window) != 1 || !is.finite(window) || window <= 0) {
    stop("window must be a single finite positive numeric value")
  }

  if (!is.numeric(step) || length(step) != 1 || !is.finite(step) || step <= 0) {
    stop("step must be a single finite positive numeric value")
  }

  # takes in a GENOME object as from PopGenome readData
  if (TRUE) {
    individuals_list <- PopGenome::get.individuals(go)
    if (length(individuals_list) == 0) {
      stop("No individuals found in the GENOME object")
    }
    numindividuals <- length(individuals_list[[1]])
    single_value <- 1/numindividuals
    if (slide) {
      slide_go <- PopGenome::sliding.window.transform(go, width=window, jump=step, type=2, whole.data=T)
      slide_go <- PopGenome::diversity.stats(slide_go, pi=T)
      slide_go <- PopGenome::neutrality.stats(slide_go, detail=T, do.R2=T)

      neutrality_list <- PopGenome::get.neutrality(slide_go)
      if (length(neutrality_list) == 0) {
        stop("PopGenome::get.neutrality(slide_go) returned an empty list")
      }
      n <- data.frame(neutrality_list[[1]])

      n$x.start <- as.numeric(sapply(strsplit(slide_go@region.names, split=" - "), "[")[1,])
      n$x.end <- as.numeric(sapply(strsplit(sub(" :", "", slide_go@region.names), split=" - "), "[")[2,])
      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(slide_go)
      if (length(diversity_list) == 0) {
        stop("PopGenome::get.diversity(slide_go) returned an empty list")
      }
      n$pi <- diversity_list[[1]][,3]

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
      go <- PopGenome::diversity.stats(go, pi=T)
      go <- PopGenome::neutrality.stats(go, detail=T, do.R2=T)

      neutrality_list <- PopGenome::get.neutrality(go)
      if (length(neutrality_list) == 0) {
        stop("PopGenome::get.neutrality(go) returned an empty list")
      }
      n <- data.frame(neutrality_list[[1]])

      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(go)
      if (length(diversity_list) == 0) {
        stop("PopGenome::get.diversity(go) returned an empty list")
      }
      n$pi <- diversity_list[[1]][,3]

      haplotype_counts <- go@region.stats@haplotype.counts
      if (length(haplotype_counts) == 0) {
        stop("go@region.stats@haplotype.counts is an empty list")
      }
      h <- haplotype_counts[[1]]

      n$n.sequences <- sum(h)
      n$n.haplotypes <- ncol(h)
      n$n.singleton.haplotypes <- length(which(h == 1))
      n$n.consensus.haplotypes <- max(h)
      if(isTRUE(is.nan(n$Fu.F_S))) n$Fu.F_S <- afufs(n$n.sequences, n$n.haplotypes, n$pi)
    }
    return(n)
  }
}
