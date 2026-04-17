SENTINEL'S JOURNAL - CRITICAL LEARNINGS ONLY

## 2026-04-15 - Command Injection and Insecure Temp Permissions in R
**Vulnerability:** Shell command injection via unsanitized file paths in `system()` calls, and world-readable temporary directories.
**Learning:** Using `system()` with `paste()` to construct shell commands is highly risky. Even with `shQuote()`, shell redirection and argument handling can be problematic. Temporary directories created with default permissions can leak sensitive data.
**Prevention:** Always use `system2()` with arguments passed as a character vector. Use the `--` flag to separate command options from positional arguments (like file paths) to prevent option injection. Set `mode = "0700"` when creating temporary directories with `dir.create()` to ensure they are private.

## 2024-05-24 - Resource Exhaustion and Input Validation in R
**Vulnerability:** Denial of Service via temporary directory leakage and potential R process blocking on special files.
**Learning:** Failing to clean up temporary directories in R (especially those with restricted permissions) can lead to disk space and inode exhaustion. Using `file()` on untrusted paths can block the R process indefinitely if the path points to a named pipe (FIFO).
**Prevention:** Always use `on.exit(unlink(hf.tempdir, recursive = TRUE))` immediately after directory creation to guarantee cleanup. Validate inputs to ensure they are single character strings and regular files (not directories or special files) before opening connections or performing operations.
