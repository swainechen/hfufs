SENTINEL'S JOURNAL - CRITICAL LEARNINGS ONLY

## 2026-04-15 - Command Injection and Insecure Temp Permissions in R
**Vulnerability:** Shell command injection via unsanitized file paths in `system()` calls, and world-readable temporary directories.
**Learning:** Using `system()` with `paste()` to construct shell commands is highly risky. Even with `shQuote()`, shell redirection and argument handling can be problematic. Temporary directories created with default permissions can leak sensitive data.
**Prevention:** Always use `system2()` with arguments passed as a character vector. Use the `--` flag to separate command options from positional arguments (like file paths) to prevent option injection. Set `mode = "0700"` when creating temporary directories with `dir.create()` to ensure they are private.

## 2024-05-24 - Resource Exhaustion and Input Validation in R
**Vulnerability:** Denial of Service via temporary directory leakage and potential R process blocking on special files.
**Learning:** Failing to clean up temporary directories in R (especially those with restricted permissions) can lead to disk space and inode exhaustion. Using `file()` on untrusted paths can block the R process indefinitely if the path points to a named pipe (FIFO).
**Prevention:** Always use `on.exit(unlink(hf.tempdir, recursive = TRUE), add = TRUE)` immediately after directory creation to guarantee cleanup without overwriting existing handlers. Validate inputs to ensure they are single character strings and regular files (not directories or special files) before opening connections or performing operations.

## 2025-05-15 - Portability and Security Theater in R System Calls
**Vulnerability:** Use of non-portable system calls for file validation when native, safer alternatives exist.
**Learning:** Attempting to use `system2("test", args = c("-f", path))` for file validation introduces a platform dependency (POSIX) that breaks R package portability (e.g., on Windows) and adds unnecessary complexity.
**Prevention:** Use `utils::file_test("-f", path)` for robust, portable validation that a path is a regular file. This avoids blocking on FIFOs while maintaining cross-platform compatibility and avoiding the risks associated with external system calls.

## 2024-05-25 - Resource Exhaustion (DoS) via Vector Allocation and Missing Validation
**Vulnerability:** Large vector allocations using the '0:(n-1)' pattern and 'rep(0, k)' can lead to out-of-memory (OOM) crashes if 'n' or 'k' are controlled by untrusted input. Missing validation on complex objects like PopGenome GENOME objects can lead to unexpected failures.
**Learning:** Even if 'n' is not technically "user input" in a web sense, in a library context, unbounded parameters that drive memory allocation are a DoS risk. 'sum(log(theta + (0:(n-1))))' is a common but dangerous pattern for large 'n'.
**Prevention:** Use mathematically equivalent but resource-efficient functions like 'lgamma' to avoid large intermediate vector allocations. Implement reasonable upper bounds on dimensions and iteration counts. Use robust type checking and 'tryCatch' when validating complex objects that rely on specific internal structures or summary outputs.

## 2024-05-26 - Defensive R Programming and Namespace Security
**Vulnerability:** Application crashes due to 'missing value where TRUE/FALSE needed' when comparisons result in NA/NaN, and potential namespace masking of external library calls.
**Learning:** R's default behavior for logical comparisons with NA is to return NA, which causes if() and other control flow to fail. Unqualified function calls can be 'masked' by functions in the global environment or other packages.
**Prevention:** Wrap logical conditions in if() statements with isTRUE() to handle NA/NaN results safely. Always prefix external library calls with their namespace (e.g., PopGenome::function) to ensure the intended code is executed and provide defense in depth against namespace pollution.

## 2024-05-27 - Functional Regression via Security Over-Hardening
**Vulnerability:** Potential path traversal or special character handling issues in temporary file staging.
**Learning:** Attempting to mitigate risks by renaming all input files to a generic "input.fasta" broke downstream analysis in `PopGenome`, which relies on filenames as sample identifiers. Security hardening must not compromise the core functional requirements of the application.
**Prevention:** Use `basename()` to safely extract filenames for staging in temporary directories, preserving metadata required by external parsers. Balance security with performance by using `file.symlink` with a fallback to `file.copy`, and use native R decompression (`gzfile`) to avoid external dependencies while maintaining portability.
