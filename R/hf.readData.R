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

  # Use utils::file_test("-f", ...) to ensure it is a regular file.
  # This is more robust than file.exists() as it excludes special files like FIFOs
  # that could cause the process to block, while remaining portable across platforms.
  if (utils::file_test("-f", fasta_file)) {
    fasta_file <- normalizePath(fasta_file)

    # readData wants a clean directory with sequences
    # get a new subdir in case there are other temp files already.
    # We check the R version as tempdir(check = TRUE) was introduced in R 4.0.0.
    hf.tempdir_root <- if (getRversion() >= "4.0.0") tempdir(check = TRUE) else tempdir()
    hf.tempdir <- tempfile(tmpdir = hf.tempdir_root)
    iter <- 0
    while ((file.exists(hf.tempdir) || dir.exists(hf.tempdir)) && iter < 100) {
      hf.tempdir <- tempfile(tmpdir = hf.tempdir_root)
      iter <- iter + 1
    }
    if (!dir.create(hf.tempdir, mode = "0700")) {
      stop(paste0("Failed to create temporary directory: ", hf.tempdir))
    }
    # Ensure temporary directory is cleaned up on exit to prevent resource leaks.
    # We use add = TRUE to avoid overwriting any existing exit handlers.
    # We only register this if dir.create succeeded.
    on.exit(unlink(hf.tempdir, recursive = TRUE), add = TRUE)

    # Preserve the original filename to maintain sample identifiers in PopGenome.
    # basename() is used to prevent path traversal when constructing the temporary path.
    orig_basename <- basename(fasta_file)
    max_size <- 2 * 1024 * 1024 * 1024 # 2GB limit for DoS protection (disk exhaustion)

    if (grepl("\\.gz$", fasta_file, ignore.case = TRUE)) {
      # Decompress while stripping .gz extension
      hf.tempfile <- file.path(hf.tempdir, sub("\\.gz$", "", orig_basename, ignore.case = TRUE))
      # Native R decompression with size limit to prevent DoS.
      # This avoids dependency on external 'zcat' and is more portable.
      con_in <- gzfile(fasta_file, "rb")
      con_out <- file(hf.tempfile, "wb")
      total_bytes <- 0
      chunk_size <- 10 * 1024 * 1024 # 10MB chunks
      tryCatch({
        while (TRUE) {
          chunk <- readBin(con_in, "raw", n = chunk_size)
          if (length(chunk) == 0) break
          total_bytes <- total_bytes + length(chunk)
          if (total_bytes > max_size) {
            stop("Decompressed file exceeds 2GB limit (DoS protection)")
          }
          writeBin(chunk, con_out)
        }
      }, finally = {
        close(con_out)
        close(con_in)
      })
    } else {
      hf.tempfile <- file.path(hf.tempdir, orig_basename)
      # Use file.symlink for performance with large genomic files, as per bioinformatics standards.
      file_info <- file.info(fasta_file)
      if (isTRUE(file_info$size > max_size)) {
        stop("File exceeds 2GB limit (DoS protection)")
      }
      if (!file.symlink(fasta_file, hf.tempfile)) {
        # Fallback to file.copy if symlink fails (e.g., on some Windows configurations)
        if (!file.copy(fasta_file, hf.tempfile, overwrite = TRUE)) {
          stop("Failed to link or copy fasta file to temporary directory")
        }
      }
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
