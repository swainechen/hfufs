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
  # Dependency check
  if (!base::requireNamespace("PopGenome", quietly = TRUE)) {
    base::stop("The 'PopGenome' package is required but not installed.")
  }

  if (!base::is.character(fasta_file) || base::length(fasta_file) != 1) {
    base::stop("fasta_file must be a single character string")
  }

  # Use utils::file_test("-f", ...) to ensure it is a regular file.
  # This is more robust than file.exists() as it excludes special files like FIFOs
  # that could cause the process to block, while remaining portable across platforms.
  if (utils::file_test("-f", fasta_file)) {
    fasta_file <- base::normalizePath(fasta_file)

    # readData wants a clean directory with sequences
    # get a new subdir in case there are other temp files already.
    # We check the R version as tempdir(check = TRUE) was introduced in R 4.0.0.
    hf.tempdir_root <- if (base::getRversion() >= "4.0.0") base::tempdir(check = TRUE) else base::tempdir()
    hf.tempdir <- base::tempfile(tmpdir = hf.tempdir_root)
    iter <- 0
    while ((base::file.exists(hf.tempdir) || base::dir.exists(hf.tempdir)) && iter < 100) {
      hf.tempdir <- base::tempfile(tmpdir = hf.tempdir_root)
      iter <- iter + 1
    }
    if (!base::dir.create(hf.tempdir, mode = "0700")) {
      base::stop(base::paste0("Failed to create temporary directory: ", hf.tempdir))
    }
    # Ensure temporary directory is cleaned up on exit to prevent resource leaks.
    # We use add = TRUE to avoid overwriting any existing exit handlers.
    # We only register this if dir.create succeeded.
    base::on.exit(base::unlink(hf.tempdir, recursive = TRUE), add = TRUE)

    # Preserve the original filename to maintain sample identifiers in PopGenome.
    # basename() is used to prevent path traversal when constructing the temporary path.
    orig_basename <- base::basename(fasta_file)
    if (orig_basename == "" || orig_basename == "." || orig_basename == "..") {
      orig_basename <- "input.fasta"
    }
    max_size <- 2 * 1024 * 1024 * 1024 # 2GB limit for DoS protection (disk exhaustion)

    if (base::grepl("\\.gz$", fasta_file, ignore.case = TRUE)) {
      # Decompress while stripping .gz extension
      decompressed_basename <- base::sub("\\.gz$", "", orig_basename, ignore.case = TRUE)
      if (decompressed_basename == "" || decompressed_basename == "." || decompressed_basename == "..") {
        decompressed_basename <- "input.fasta"
      }
      hf.tempfile <- base::file.path(hf.tempdir, decompressed_basename)

      # Native R decompression with size limit to prevent DoS.
      # This avoids dependency on external 'zcat' and is more portable.
      con_in <- base::gzfile(fasta_file, "rb")
      base::tryCatch({
        con_out <- base::file(hf.tempfile, "wb")
        base::tryCatch({
          total_bytes <- 0
          chunk_size <- 10 * 1024 * 1024 # 10MB chunks
          while (TRUE) {
            chunk <- base::readBin(con_in, "raw", n = chunk_size)
            if (base::length(chunk) == 0) break
            total_bytes <- total_bytes + base::length(chunk)
            if (total_bytes > max_size) {
              base::stop("Decompressed file exceeds 2GB limit (DoS protection)")
            }
            base::writeBin(chunk, con_out)
          }
        }, finally = {
          base::close(con_out)
        })
      }, finally = {
        base::close(con_in)
      })
    } else {
      hf.tempfile <- base::file.path(hf.tempdir, orig_basename)
      # Use file.symlink for performance with large genomic files, as per bioinformatics standards.
      file_info <- base::file.info(fasta_file)
      if (base::isTRUE(file_info$size > max_size)) {
        base::stop("File exceeds 2GB limit (DoS protection)")
      }
      if (!base::file.symlink(fasta_file, hf.tempfile)) {
        # Fallback to file.copy if symlink fails (e.g., on some Windows configurations)
        if (!base::file.copy(fasta_file, hf.tempfile, overwrite = TRUE)) {
          base::stop("Failed to link or copy fasta file to temporary directory")
        }
      }
    }

    res <- base::tryCatch(
      {
        PopGenome::readData(hf.tempdir)
      },
      error = function(e) {
        base::warning("Error running PopGenome::readData - maybe check if PopGenome is installed")
        return(NULL)
      }
    )

    return(res)
  } else {
    base::stop(base::paste0("fasta_file '", fasta_file, "' does not exist or is a directory"))
  }
}
