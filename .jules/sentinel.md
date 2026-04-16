SENTINEL'S JOURNAL - CRITICAL LEARNINGS ONLY

## 2026-04-15 - Command Injection and Insecure Temp Permissions in R
**Vulnerability:** Shell command injection via unsanitized file paths in `system()` calls, and world-readable temporary directories.
**Learning:** Using `system()` with `paste()` to construct shell commands is highly risky. Even with `shQuote()`, shell redirection and argument handling can be problematic. Temporary directories created with default permissions can leak sensitive data.
**Prevention:** Always use `system2()` with arguments passed as a character vector. Use the `--` flag to separate command options from positional arguments (like file paths) to prevent option injection. Set `mode = "0700"` when creating temporary directories with `dir.create()` to ensure they are private.
