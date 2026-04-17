#' Simple interface to read an aligned fasta file into PopGenome
#'
#' Reads a given fasta file into a PopGenome GENOME object.
#' We assume that the fasta file is aligned. This will not do anything
#' related to associated gff files.
#' This uses the defaults of the PopGenome::readData function:
#'   populations = FALSE
#'   format = "fasta"
#' For other information, see the PopGenome documentation.
#' If you need something more complex than a simple fasta file loaded, it is
#' recommended to use the PopGenome functions directly.
#'
#' Note this will load the PopGenome package if needed.
#'
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
hf.readData <- function(fasta_file) {
  if (!is.character(fasta_file) || length(fasta_file) != 1) {
    stop("fasta_file must be a single character string")
  }

  if (file.exists(fasta_file) && !dir.exists(fasta_file)) {
    fasta_file <- normalizePath(fasta_file)

    # readData wants a clean directory with sequences
    # get a new subdir in case there are other temp files already
    hf.tempdir <- tempfile(tmpdir = tempdir(check = TRUE))
    while (file.exists(hf.tempdir) || dir.exists(hf.tempdir)) {
      hf.tempdir <- tempfile(tmpdir = tempdir(check = TRUE))
    }
    dir.create(hf.tempdir, mode = "0700")
    # Ensure temporary directory is cleaned up on exit to prevent resource leaks
    on.exit(unlink(hf.tempdir, recursive = TRUE))

    hf.tempfile <- file.path(hf.tempdir, basename(fasta_file))

    if (grepl("\\.gz$", fasta_file, ignore.case = TRUE)) {
      hf.tempfile <- sub("\\.gz$", "", hf.tempfile, ignore.case = TRUE)
      # Use system2 for more secure command execution, avoiding shell interpretation
      # and protecting against option injection with the '--' flag.
      zcat_res <- system2("zcat", args = c("--", fasta_file), stdout = hf.tempfile)
      if (zcat_res != 0) {
        stop("Failed to decompress fasta file using zcat")
      }
    } else {
      file.symlink(fasta_file, hf.tempfile)
    }

    res <- tryCatch(
      {
        PopGenome::readData(hf.tempdir)
      },
      error = function(e) {
        warning("Error running PopGenome::readData - maybe check if PopGenome is installed")
        return(NULL)
      }
    )

    return(res)
  } else {
    stop(paste0("fasta_file '", fasta_file, "' does not exist or is a directory"))
  }
}
