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
#' @usage hf.alignment.stats(go, slide = FALSE, window = 1000, step = 500)
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
  if (base::isTRUE(base::missing(go))) {
    base::stop("go must be provided")
  }

  # Check if go is a GENOME object safely
  if (base::isTRUE(!base::inherits(go, "GENOME"))) {
    base::stop("go must be a PopGenome GENOME object")
  }

  if (base::isTRUE(!base::is.logical(slide) || base::length(slide) != 1 || base::is.na(slide))) {
    base::stop("slide must be a single non-NA logical value")
  }

  if (base::isTRUE(!base::is.numeric(window) || base::length(window) != 1 || !base::is.finite(window) || window < 1 || window > 2000000000 || window %% 1 != 0)) {
    base::stop("window must be a single finite numeric whole number between 1 and 2,000,000,000")
  }

  if (base::isTRUE(!base::is.numeric(step) || base::length(step) != 1 || !base::is.finite(step) || step < 1 || step > 2000000000 || step %% 1 != 0)) {
    base::stop("step must be a single finite numeric whole number between 1 and 2,000,000,000")
  }

  # DoS Protection: Limit the number of regions to 1,000,000
  if (base::isTRUE(base::length(go@n.sites) > 1000000)) {
    base::stop("The GENOME object contains > 1,000,000 regions (DoS protection)")
  }

  # DoS Protection: Limit the number of sliding windows to 1,000,000
  if (base::isTRUE(slide)) {
    # PopGenome treats n.sites as a list when multiple regions are loaded.
    # We use unlist() before as.numeric() to safely handle list inputs and prevent coercion crashes.
    n_sites <- base::sum(base::as.numeric(base::unlist(go@n.sites)))
    if (base::isTRUE(!base::is.finite(n_sites) || n_sites < 1)) {
      base::stop("The GENOME object contains no sites or invalid site counts")
    }
    if (base::isTRUE(window > n_sites)) {
      base::stop("The requested window size exceeds the total number of sites in the GENOME object")
    }
    num_windows <- base::ceiling((n_sites - window + 1) / step)
    if (base::isTRUE(!base::is.finite(num_windows) || num_windows > 1000000 || num_windows < 1)) {
      base::stop("The requested sliding window parameters would generate > 1,000,000 windows or an invalid number of windows (DoS protection)")
    }
  }

  # takes in a GENOME object as from PopGenome readData
  if (base::isTRUE(base::inherits(go, "GENOME"))) {
    individuals_list <- PopGenome::get.individuals(go)
    if (base::isTRUE(base::length(individuals_list) == 0)) {
      base::stop("No individuals found in the GENOME object")
    }

    # DoS Protection: Limit the number of populations to 1,000 to prevent resource exhaustion.
    num_populations <- base::length(individuals_list)
    if (base::isTRUE(num_populations > 1000)) {
      base::stop("The GENOME object contains > 1,000 populations (DoS protection)")
    }

    pop_sizes <- base::numeric(base::length(individuals_list))
    # Security: Iteratively validate all populations to prevent hidden large populations
    # from bypassing resource limits and causing a Denial of Service (DoS).
    for (i in base::seq_along(individuals_list)) {
      pop_sizes[i] <- base::length(individuals_list[[i]])
      if (pop_sizes[i] == 0) {
        base::stop(base::paste0("The GENOME object contains no individuals in population ", i))
      }
      # DoS Protection: Limit the number of individuals to 1,000,000 per population.
      if (base::isTRUE(pop_sizes[i] > 1000000)) {
        base::stop(base::paste0("The GENOME object contains > 1,000,000 individuals in population ", i, " (DoS protection)"))
      }
    }

    # DoS Protection: Enforce a cumulative limit of 10,000,000 individuals across all
    # populations to prevent overall memory exhaustion.
    if (base::isTRUE(base::sum(pop_sizes) > 10000000)) {
      base::stop("The GENOME object contains > 10,000,000 total individuals (DoS protection)")
    }

    # Security: Correctly assign from internal state and ensure we do not use undefined
    # external variables (e.g., 'num_ind') which could lead to application crashes.
    numindividuals <- pop_sizes[1]

    if (base::isTRUE(slide)) {
      slide_go <- PopGenome::sliding.window.transform(go, width=window, jump=step, type=2, whole.data=TRUE)
      slide_go <- PopGenome::diversity.stats(slide_go, pi=TRUE)
      slide_go <- PopGenome::neutrality.stats(slide_go, detail=TRUE, do.R2=TRUE)

      neutrality_list <- PopGenome::get.neutrality(slide_go)
      if (base::isTRUE(base::length(neutrality_list) == 0)) {
        base::stop("PopGenome::get.neutrality(slide_go) returned an empty list")
      }
      n <- base::data.frame(neutrality_list[[1]], stringsAsFactors = FALSE)

      # Robustly extract start/end coordinates. We use vapply to ensure numeric return types
      # and prevent unpredictable simplification, providing better type safety.
      # We check the length of splits to prevent out-of-bounds access.
      region_splits <- base::strsplit(base::sub(" :$", "", slide_go@region.names), split=" - ")
      n$x.start <- base::vapply(region_splits, function(x) {
        if (base::length(x) < 1) return(base::as.numeric(NA))
        base::as.numeric(x[[1]])
      }, FUN.VALUE = 0.0)
      n$x.end <- base::vapply(region_splits, function(x) {
        if (base::length(x) < 2) return(base::as.numeric(NA))
        base::as.numeric(x[[2]])
      }, FUN.VALUE = 0.0)
      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(slide_go)
      if (base::isTRUE(base::length(diversity_list) == 0)) {
        base::stop("PopGenome::get.diversity(slide_go) returned an empty list")
      }
      # Security: Validate diversity matrix dimensions before indexing to prevent out-of-bounds errors.
      if (base::isTRUE(base::is.null(base::ncol(diversity_list[[1]])) || base::ncol(diversity_list[[1]]) < 3)) {
        base::stop("PopGenome::get.diversity(slide_go) returned a matrix with insufficient columns (expected at least 3 for pi)")
      }

      # Security: Verify that neutrality and diversity stats have the same number of regions.
      if (base::isTRUE(base::nrow(n) != base::nrow(diversity_list[[1]]))) {
        base::stop("Mismatch in number of regions between neutrality and diversity statistics")
      }
      n$pi <- diversity_list[[1]][,3]

      h_list <- slide_go@region.stats@haplotype.counts
      if (base::isTRUE(base::length(h_list) == 0)) {
        base::stop("slide_go@region.stats@haplotype.counts is an empty list")
      }
      # Security: Verify that haplotype counts have the same number of regions.
      if (base::isTRUE(base::nrow(n) != base::length(h_list))) {
        base::stop("Mismatch in number of regions between statistics and haplotype counts")
      }
      # Use vapply for type safety. We index [[1]] to get the first population, ensuring
      # consistency with neutrality and diversity stats. We handle both vector and matrix
      # return structures from PopGenome.
      n$n.sequences <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        base::sum(reg[[1]])
      }, FUN.VALUE = 0.0)
      n$n.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        x <- reg[[1]]
        nc <- base::ncol(x)
        if (base::is.null(nc)) return(base::as.numeric(base::length(base::which(x > 0))))
        return(base::as.numeric(nc))
      }, FUN.VALUE = 0.0)
      n$n.singleton.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        base::as.numeric(base::length(base::which(reg[[1]] == 1)))
      }, FUN.VALUE = 0.0)
      n$n.consensus.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]]) || base::length(reg[[1]]) == 0) return(0.0)
        base::as.numeric(base::max(reg[[1]]))
      }, FUN.VALUE = 0.0)
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
      if (base::isTRUE(base::length(neutrality_list) == 0)) {
        base::stop("PopGenome::get.neutrality(go) returned an empty list")
      }
      n <- base::data.frame(neutrality_list[[1]], stringsAsFactors = FALSE)

      n$numindividuals <- numindividuals

      diversity_list <- PopGenome::get.diversity(go)
      if (base::isTRUE(base::length(diversity_list) == 0)) {
        base::stop("PopGenome::get.diversity(go) returned an empty list")
      }
      # Security: Validate diversity matrix dimensions before indexing to prevent out-of-bounds errors.
      if (base::isTRUE(base::is.null(base::ncol(diversity_list[[1]])) || base::ncol(diversity_list[[1]]) < 3)) {
        base::stop("PopGenome::get.diversity(go) returned a matrix with insufficient columns (expected at least 3 for pi)")
      }

      # Security: Verify that neutrality and diversity stats have the same number of regions.
      if (base::isTRUE(base::nrow(n) != base::nrow(diversity_list[[1]]))) {
        base::stop("Mismatch in number of regions between neutrality and diversity statistics")
      }
      n$pi <- diversity_list[[1]][,3]

      h_list <- go@region.stats@haplotype.counts
      if (base::isTRUE(base::length(h_list) == 0)) {
        base::stop("go@region.stats@haplotype.counts is an empty list")
      }
      # Security: Verify that haplotype counts have the same number of regions.
      if (base::isTRUE(base::nrow(n) != base::length(h_list))) {
        base::stop("Mismatch in number of regions between statistics and haplotype counts")
      }

      # Use vapply for type safety. We index [[1]] to get the first population, ensuring
      # consistency with neutrality and diversity stats. We handle both vector and matrix
      # return structures from PopGenome.
      n$n.sequences <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        base::sum(reg[[1]])
      }, FUN.VALUE = 0.0)
      n$n.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        x <- reg[[1]]
        nc <- base::ncol(x)
        if (base::is.null(nc)) return(base::as.numeric(base::length(base::which(x > 0))))
        return(base::as.numeric(nc))
      }, FUN.VALUE = 0.0)
      n$n.singleton.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]])) return(0.0)
        base::as.numeric(base::length(base::which(reg[[1]] == 1)))
      }, FUN.VALUE = 0.0)
      n$n.consensus.haplotypes <- base::vapply(h_list, function(reg) {
        if (base::is.null(reg) || base::length(reg) < 1 || base::is.null(reg[[1]]) || base::length(reg[[1]]) == 0) return(0.0)
        base::as.numeric(base::max(reg[[1]]))
      }, FUN.VALUE = 0.0)
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
