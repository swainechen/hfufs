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
  fasta_file <- normalizePath(fasta_file)
  if (file.exists(fasta_file)) {
    filehandle <- file(fasta_file)
    file_checks <- summary(filehandle)

    # readData wants a clean directory with sequences
    # get a new subdir in case there are other temp files already
    hf.tempdir <- tempfile(tmpdir=tempdir(check=T))
    while (file.exists(hf.tempdir) | dir.exists(hf.tempdir)) {
      hf.tempdir <- tempfile(tmpdir=tempdir(check=T))
    }
    dir.create(hf.tempdir)
    hf.tempfile <- file.path(hf.tempdir, basename(fasta_file))

    if (file_checks$class == "gzfile") {
      hf.tempfile <- sub(".gz$", "", hf.tempfile)
      system(paste("zcat ", fasta_file, " > ", hf.tempfile))
    } else {
      file.symlink(fasta_file, hf.tempfile)
    }
    close(filehandle)

    return <- tryCatch (
      { PopGenome::readData(hf.tempdir) },
      error = function(e) { warning("Error running PopGenome::readData - maybe check if PopGenome is installed") }
    )

    unlink(hf.tempfile)
    # leave this here, should be empty though
    # unlink(hf.tempdir, recursive = T)

    return(return)
  }
}
