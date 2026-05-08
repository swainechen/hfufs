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
  if (!base::requireNamespace("PopGenome", quietly = TRUE)) {
    base::stop("The 'PopGenome' package is required but not installed.")
  }

  # Input validation
  if (base::missing(go)) {
    base::stop("go must be provided")
  }

  # Check if go is a GENOME object safely
  is_genome <- base::tryCatch({
    summ <- base::summary(go)
    base::is.character(summ) && base::length(summ) >= 2 && base::isTRUE(summ[2] == "GENOME")
  }, error = function(e) FALSE)

  if (!is_genome) {
    base::stop("go must be a PopGenome GENOME object")
  }

  if (!base::is.logical(slide) || base::length(slide) != 1 || base::is.na(slide)) {
    base::stop("slide must be a single non-NA logical value")
  }

  if (!base::is.numeric(window) || base::length(window) != 1 || !base::is.finite(window) || window < 1 || window > 2000000000) {
    base::stop("window must be a single finite numeric value between 1 and 2,000,000,000")
  }

  if (!base::is.numeric(step) || base::length(step) != 1 || !base::is.finite(step) || step < 1 || step > 2000000000) {
    base::stop("step must be a single finite numeric value between 1 and 2,000,000,000")
  }

  # DoS Protection: Limit the number of sliding windows to 1,000,000
  if (base::isTRUE(slide)) {
    # PopGenome treats n.sites as a list when multiple regions are loaded.
    n_sites <- base::sum(base::as.numeric(go@n.sites))
    num_windows <- base::ceiling((n_sites - window + 1) / step)
    if (base::isTRUE(num_windows > 1000000)) {
      base::stop("The requested sliding window parameters would generate > 1,000,000 windows (DoS protection)")
    }
  }

  # takes in a GENOME object as from PopGenome readData
  if (TRUE) {
    individuals_list <- PopGenome::get.individuals(go)
    if (base::length(individuals_list) == 0) {
      base::stop("No individuals found in the GENOME object")
    }
    numindividuals <- base::length(individuals_list[[1]])
    if (numindividuals == 0) {
      base::stop("The GENOME object contains no individuals in the first population")
    }

    if (slide) {
      slide_go <- PopGenome::sliding.window.transform(go, width=window, jump=step, type=2, whole.data=TRUE)
      slide_go <- PopGenome::diversity.stats(slide_go, pi=TRUE)
      slide_go <- PopGenome::neutrality.stats(slide_go, detail=TRUE, do.R2=TRUE)

      neutrality_list <- PopGenome::get.neutrality(slide_go)
      if (base::length(neutrality_list) == 0) {
        base::stop("PopGenome::get.neutrality(slide_go) returned an empty list")
      }
      n <- base::data.frame(neutrality_list[[1]])

      # Robustly extract start/end coordinates. We avoid [1,] and [2,] indexing which can be
      # fragile depending on how sapply simplifies the list of vectors.
      region_splits <- base::strsplit(base::sub(" :$", "", slide_go@region.names), split=" - ")
      n$x.start <- base::as.numeric(base::sapply(region_splits, "[", 1))
      n$x.end <- base::as.numeric(base::sapply(region_splits, "[", 2))
      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(slide_go)
      if (base::length(diversity_list) == 0) {
        base::stop("PopGenome::get.diversity(slide_go) returned an empty list")
      }
      n$pi <- diversity_list[[1]][,3]

      h <- slide_go@region.stats@haplotype.counts
      if (base::length(h) == 0) {
        base::stop("slide_go@region.stats@haplotype.counts is an empty list")
      }
      n$n.sequences <- base::as.numeric(base::lapply(h, base::sum))
      templist <- base::lapply(h, base::ncol)
      templist[base::sapply(templist, base::is.null)] <- NA
      n$n.haplotypes <- base::as.numeric(templist)
      n$n.singleton.haplotypes <- base::as.numeric(base::lapply(h, function(x) base::length(base::which(x==1))))
      for(i in base::which(base::is.na(n$Fu.F_S))) {
        n$Fu.F_S[i] <- base::tryCatch(
          { afufs(n$n.sequences[i], n$n.haplotypes[i], n$pi[i]) },
          error = function(e) { return(NaN) }
        )
      }
    } else {
      go <- PopGenome::diversity.stats(go, pi=TRUE)
      go <- PopGenome::neutrality.stats(go, detail=TRUE, do.R2=TRUE)

      neutrality_list <- PopGenome::get.neutrality(go)
      if (base::length(neutrality_list) == 0) {
        base::stop("PopGenome::get.neutrality(go) returned an empty list")
      }
      n <- base::data.frame(neutrality_list[[1]])

      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(go)
      if (base::length(diversity_list) == 0) {
        base::stop("PopGenome::get.diversity(go) returned an empty list")
      }
      n$pi <- diversity_list[[1]][,3]

      haplotype_counts <- go@region.stats@haplotype.counts
      if (base::length(haplotype_counts) == 0) {
        base::stop("go@region.stats@haplotype.counts is an empty list")
      }
      h <- haplotype_counts[[1]]

      n$n.sequences <- base::sum(h)
      n$n.haplotypes <- base::ncol(h)
      n$n.singleton.haplotypes <- base::length(base::which(h == 1))
      n$n.consensus.haplotypes <- base::max(h)
      for(i in base::which(base::is.na(n$Fu.F_S))) {
        n$Fu.F_S[i] <- base::tryCatch(
          { afufs(n$n.sequences[i], n$n.haplotypes[i], n$pi[i]) },
          error = function(e) { return(NaN) }
        )
      }
    }
    return(n)
  }
}
