#' Simple interface to read an aligned fasta file into PopGenome
#'
#' Reads a given fasta file into a PopGenome GENOME object.
#'
#'
#' @export
#'
#' @examples
#' \dontrun {
#' fasta_file <- "name_of_aligned.fasta"
#' hf.readData(fasta_file)
#' }
#'
hf.readData <- function(fasta_file) {
  fasta_file <- normalizePath(fasta_file)
  hf.tempdir <- tempdir(check=T)
  hf.tempfile <- tempfile(tmpdir=hf.tempdir)	# fallback temp file
  if (file.exists(fasta_file)) {
    file_checks <- summary(file(fasta_file))
    hf.tempfile <- file.path(hf.tempdir, basename(fasta_file))
    if (file_checks$class == "gzfile") {
      hf.tempfile <- sub(".gz$", "", hf.tempfile)
      system(paste("zcat ", fasta_file, " > ", hf.tempfile))
    } else {
      file.symlink(fasta_file, hf.tempfile)
    }
    return <- PopGenome::readData(hf.tempdir)
    unlink(hf.tempfile)
    return(return)
  }
}
