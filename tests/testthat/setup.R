# Pre-create the R6 reference worktree before parallel workers start.
# In testthat 3e parallel mode, setup.R runs in the main process while
# helper_*.R runs in each worker. Using a fixed path (not tempdir())
# ensures all workers find the same pre-created worktree.
repo_dir <- tryCatch({
  gcd <- system2("git", c("rev-parse", "--git-common-dir"),
                  stdout = TRUE, stderr = FALSE)
  normalizePath(file.path(gcd, ".."))
}, error = function(e) NULL)

if (!is.null(repo_dir)) {
  ref_source <- file.path(repo_dir, ".git", "mvsusieR_r6_ref")
  if (!dir.exists(ref_source)) {
    system2("git", c("-C", repo_dir, "worktree", "remove", "--force",
                     ref_source), stdout = FALSE, stderr = FALSE)
    system2("git", c("-C", repo_dir, "worktree", "add", "--detach",
                     ref_source, "master"), stdout = FALSE, stderr = FALSE)
  }
}
