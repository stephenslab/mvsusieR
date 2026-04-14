# Diagnose why mvSuSiE gives all-zero PIPs on the vignette data.
# Exact same data as vignette (N=100, R=128, L=1).

library(fsusieR)
library(mvsusieR)

set.seed(1)
data(N3finemapping)
X_full <- N3finemapping$X[1:100, ]

effect <- 1.2 * cos((1:128) / 128 * 3 * pi)
effect[effect < 0] <- 0
effect[1:40] <- 0

true_pos <- 700
y <- matrix(0, 100, 128)
for (i in 1:100)
  y[i, ] <- X_full[i, true_pos] * effect + rnorm(128)

col_var <- apply(X_full, 2, var)
keep <- which(col_var > 0)
X <- X_full[, keep]
true_pos_filt <- which(keep == true_pos)

cat("N =", nrow(X), ", J =", ncol(X), ", R =", ncol(y), "\n")
cat("True SNP (filtered index):", true_pos_filt, "\n\n")

# Priors
U_true <- outer(effect, effect)
prior_true <- create_mixture_prior(
  mixture_prior = list(matrices = list(U_true), weights = 1)
)
prior_canon <- create_mixture_prior(R = 128)

# Residual variances
V_true <- diag(128)

print_results <- function(fit, label) {
  cat(sprintf("\n=== %s ===\n", label))
  cat("PIP at true SNP:", fit$pip[true_pos_filt], "\n")
  cat("Max PIP:", max(fit$pip), "at filtered SNP", which.max(fit$pip),
      "(original:", keep[which.max(fit$pip)], ")\n")
  cat("Converged:", fit$converged, "\n")
  cat("Number of CS:", length(fit$sets$cs), "\n")
  if (length(fit$sets$cs) > 0)
    cat("CS SNPs (original):", paste(keep[fit$sets$cs[[1]]], collapse = ", "), "\n")
  cat("ELBO:", paste(round(tail(fit$elbo, 3), 2), collapse = " -> "), "\n")
}

# Test 1: true prior + true resid (fixed)
cat("Running Test 1...\n")
fit1 <- mvsusie(X, y, L = 1,
                prior_variance = prior_true,
                residual_variance = V_true,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit1, "Test 1: true prior + true resid (fixed)")

# Test 2: true prior + estimated resid (true init)
cat("\nRunning Test 2...\n")
fit2 <- mvsusie(X, y, L = 1,
                prior_variance = prior_true,
                residual_variance = V_true,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit2, "Test 2: true prior + estimated resid (true init)")

# Test 3: canon prior + true resid (fixed)
cat("\nRunning Test 3...\n")
fit3 <- mvsusie(X, y, L = 1,
                prior_variance = prior_canon,
                residual_variance = V_true,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit3, "Test 3: canon prior + true resid (fixed)")

# Test 4: true prior + estimated resid (default init)
cat("\nRunning Test 4...\n")
fit4 <- mvsusie(X, y, L = 1,
                prior_variance = prior_true,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit4, "Test 4: true prior + estimated resid")

# Test 5: true prior + fixed resid (default init)
cat("\nRunning Test 5...\n")
fit5 <- mvsusie(X, y, L = 1,
                prior_variance = prior_true,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit5, "Test 5: true prior + fixed resid (default init)")

# Test 6: canon prior + estimated resid (default init)
cat("\nRunning Test 6...\n")
fit6 <- mvsusie(X, y, L = 1,
                prior_variance = prior_canon,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit6, "Test 6: canon prior + estimated resid")

# Test 7: canon prior + fixed resid (default init)
cat("\nRunning Test 7...\n")
fit7 <- mvsusie(X, y, L = 1,
                prior_variance = prior_canon,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = TRUE,
                max_iter = 100, verbose = TRUE)
print_results(fit7, "Test 7: canon prior + fixed resid (default init)")

# Summary table
cat("\n\n=== SUMMARY ===\n")
labels <- c(
  "1: true prior   + true resid (fixed)",
  "2: true prior   + estimated resid (true init)",
  "3: canon prior  + true resid (fixed)",
  "4: true prior   + estimated resid",
  "5: true prior   + fixed resid (default init)",
  "6: canon prior  + estimated resid",
  "7: canon prior  + fixed resid (default init)"
)
cat(sprintf("%-55s  PIP@700  MaxPIP  #CS  Conv\n", "Test"))
for (i in 1:7) {
  f <- get(paste0("fit", i))
  cat(sprintf("%-55s  %7.4f  %7.4f  %3d  %s\n",
              labels[i], f$pip[true_pos_filt], max(f$pip),
              length(f$sets$cs), f$converged))
}
