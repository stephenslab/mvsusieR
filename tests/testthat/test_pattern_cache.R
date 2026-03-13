# Tests for pattern-based missingness caching and imputation correctness.
#
# Verifies: extract_missing_patterns, precompute_pattern_cache,
# impute_missing_Y (R version), impute_missing_Y_rcpp (C++ version),
# and consistency between the two implementations.


# ============================================================================
# Helper: naive per-observation imputation (no caching)
# ============================================================================
naive_impute_missing_Y <- function(Y, mu, Vinv, Y_missing) {
  N <- nrow(Y)
  R <- ncol(Y)
  Y_cov <- matrix(0, R, R)
  sum_neg_ent <- 0

  for (i in seq_len(N)) {
    miss <- which(Y_missing[i, ])
    obs  <- which(!Y_missing[i, ])
    n_miss <- length(miss)
    if (n_miss == 0) next

    Vinv_mm <- Vinv[miss, miss, drop = FALSE]
    Vinv_mm_chol <- chol(Vinv_mm)
    cov_mm <- chol2inv(Vinv_mm_chol)

    Y_cov[miss, miss] <- Y_cov[miss, miss] + cov_mm
    sum_neg_ent <- sum_neg_ent + 0.5 * (n_miss * log(1 / (2 * pi * exp(1))) +
                                         mvsusieR:::chol2ldet(Vinv_mm_chol))

    if (length(obs) > 0) {
      Vinv_mo <- Vinv[miss, obs, drop = FALSE]
      resid_obs <- Y[i, obs] - mu[i, obs]
      correction <- cov_mm %*% Vinv_mo %*% resid_obs
      Y[i, miss] <- mu[i, miss] - as.vector(correction)
    } else {
      Y[i, miss] <- mu[i, miss]
    }
  }

  list(Y = Y, Y_cov = Y_cov, sum_neg_ent_Y_miss = sum_neg_ent)
}


# ============================================================================
# 1. Pattern extraction basic
# ============================================================================
test_that("Pattern extraction identifies correct patterns from MCAR data", {
  set.seed(42)
  N <- 50; R <- 4
  Y_missing <- matrix(runif(N * R) < 0.3, N, R)

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)

  # Should have at least 1 pattern (since some rows have missing)

  expect_true(miss_info$n_patterns >= 1)
  expect_equal(length(miss_info$patterns), miss_info$n_patterns)
  expect_equal(length(miss_info$pattern_obs), miss_info$n_patterns)

  # All observation indices should be valid
  all_obs <- unlist(miss_info$pattern_obs)
  expect_true(all(all_obs >= 1 & all_obs <= N))

  # Each pattern should have at least one missing entry
  for (k in seq_len(miss_info$n_patterns)) {
    expect_true(any(miss_info$patterns[[k]]))
  }

  # No observation should appear in multiple patterns
  expect_equal(length(all_obs), length(unique(all_obs)))
})


# ============================================================================
# 2. Pattern extraction: all same pattern
# ============================================================================
test_that("Pattern extraction: single pattern when all rows have same missingness", {
  N <- 100; R <- 3
  Y_missing <- matrix(FALSE, N, R)
  Y_missing[, 2] <- TRUE  # outcome 2 always missing

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)

  expect_equal(miss_info$n_patterns, 1)
  expect_equal(length(miss_info$pattern_obs[[1]]), N)
  expect_equal(miss_info$patterns[[1]], c(FALSE, TRUE, FALSE))
})


# ============================================================================
# 3. Pattern extraction: all unique patterns
# ============================================================================
test_that("Pattern extraction: unique pattern per observation", {
  N <- 5; R <- 5
  Y_missing <- matrix(FALSE, N, R)
  for (i in seq_len(N)) Y_missing[i, i] <- TRUE

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)

  expect_equal(miss_info$n_patterns, N)
  for (k in seq_len(N)) {
    expect_equal(length(miss_info$pattern_obs[[k]]), 1)
  }
})


# ============================================================================
# 4. Pattern extraction: no missing data
# ============================================================================
test_that("Pattern extraction: empty when no missing data", {
  N <- 30; R <- 3
  Y_missing <- matrix(FALSE, N, R)

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)

  expect_equal(miss_info$n_patterns, 0)
  expect_equal(length(miss_info$patterns), 0)
  expect_equal(length(miss_info$pattern_obs), 0)
})


# ============================================================================
# 5. Cache vs naive per-observation loop
# ============================================================================
test_that("Pattern-cached imputation matches naive per-obs loop", {
  set.seed(123)
  N <- 80; R <- 4
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.5), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  # Create MCAR missingness (30%)
  Y_missing <- matrix(runif(N * R) < 0.3, N, R)
  # Ensure at least some rows have all obs (for variety)
  Y_missing[1:5, ] <- FALSE
  Y[Y_missing] <- 0

  # Naive per-observation
  naive_res <- naive_impute_missing_Y(Y, mu, Vinv, Y_missing)

  # Pattern-cached
  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  cached_res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                             pattern_cache, version = "R")

  expect_equal(cached_res$Y, naive_res$Y, tolerance = 1e-14)
  expect_equal(cached_res$Y_cov, naive_res$Y_cov, tolerance = 1e-14)
  expect_equal(cached_res$sum_neg_ent_Y_miss, naive_res$sum_neg_ent_Y_miss,
               tolerance = 1e-14)
})


# ============================================================================
# 6. Cache vs non-cache reference (formula verification)
# ============================================================================
test_that("Pattern-cached imputation matches manual conditional normal formula", {
  set.seed(456)
  N <- 20; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.3), N, R)
  V <- matrix(c(2, 0.5, 0.3,  0.5, 1.5, 0.2,  0.3, 0.2, 1.8), R, R)
  Vinv <- solve(V)

  Y_missing <- matrix(FALSE, N, R)
  Y_missing[3, 2] <- TRUE    # obs 3: outcome 2 missing
  Y_missing[7, c(1,3)] <- TRUE  # obs 7: outcomes 1,3 missing
  Y[Y_missing] <- 0

  # Manual imputation for obs 3 (outcome 2 missing)
  miss3 <- 2; obs3 <- c(1, 3)
  cov_mm_3 <- solve(matrix(Vinv[miss3, miss3], 1, 1))
  Vinv_mo_3 <- matrix(Vinv[miss3, obs3], 1, 2)
  resid3 <- Y[3, obs3] - mu[3, obs3]
  y3_imp <- mu[3, miss3] - as.numeric(cov_mm_3 %*% Vinv_mo_3 %*% resid3)

  # Manual imputation for obs 7 (outcomes 1,3 missing)
  miss7 <- c(1, 3); obs7 <- 2
  cov_mm_7 <- solve(Vinv[miss7, miss7])
  Vinv_mo_7 <- Vinv[miss7, obs7, drop = FALSE]
  resid7 <- Y[7, obs7] - mu[7, obs7]
  y7_imp <- mu[7, miss7] - as.vector(cov_mm_7 %*% Vinv_mo_7 %*% resid7)

  # Pattern-cached
  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  cached_res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                             pattern_cache, version = "R")

  expect_equal(cached_res$Y[3, 2], y3_imp, tolerance = 1e-14)
  expect_equal(cached_res$Y[7, c(1,3)], y7_imp, tolerance = 1e-14)
})


# ============================================================================
# 7. R vs Rcpp implementation
# ============================================================================
test_that("R and Rcpp imputation produce identical results", {
  set.seed(789)
  N <- 60; R <- 5
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.5), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.25, N, R)
  Y_missing[1:3, ] <- FALSE  # some fully observed
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)

  r_res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                        pattern_cache, version = "R")
  rcpp_res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                           pattern_cache, version = "Rcpp")

  expect_equal(rcpp_res$Y, r_res$Y, tolerance = 1e-12)
  expect_equal(rcpp_res$Y_cov, r_res$Y_cov, tolerance = 1e-12)
  expect_equal(rcpp_res$sum_neg_ent_Y_miss, r_res$sum_neg_ent_Y_miss,
               tolerance = 1e-12)
})


# ============================================================================
# 8. Rcpp vs naive per-obs (end-to-end)
# ============================================================================
test_that("Rcpp imputation matches naive per-observation loop", {
  set.seed(321)
  N <- 40; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.5), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.35, N, R)
  Y[Y_missing] <- 0

  naive_res <- naive_impute_missing_Y(Y, mu, Vinv, Y_missing)

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  rcpp_res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                           pattern_cache, version = "Rcpp")

  expect_equal(rcpp_res$Y, naive_res$Y, tolerance = 1e-12)
  expect_equal(rcpp_res$Y_cov, naive_res$Y_cov, tolerance = 1e-12)
  expect_equal(rcpp_res$sum_neg_ent_Y_miss, naive_res$sum_neg_ent_Y_miss,
               tolerance = 1e-12)
})


# ============================================================================
# 9. Y_cov symmetry
# ============================================================================
test_that("Y_cov is symmetric after imputation", {
  set.seed(111)
  N <- 50; R <- 4
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.3, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  expect_equal(res$Y_cov, t(res$Y_cov), tolerance = 1e-15)
})


# ============================================================================
# 10. Y_cov positive semidefinite
# ============================================================================
test_that("Y_cov is positive semidefinite", {
  set.seed(222)
  N <- 50; R <- 4
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.3, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  evals <- eigen(res$Y_cov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals >= -1e-14))
})


# ============================================================================
# 11. Imputed Y: observed entries unchanged
# ============================================================================
test_that("Observed entries of Y are unchanged after imputation", {
  set.seed(333)
  N <- 40; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.3, N, R)
  Y_orig <- Y
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  # Observed entries should be exactly the same as the original (zero-filled) Y
  expect_equal(res$Y[!Y_missing], Y[!Y_missing])
})


# ============================================================================
# 12. Single missing entry: scalar conditional normal
# ============================================================================
test_that("Single missing entry matches scalar conditional normal formula", {
  set.seed(444)
  N <- 10; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- matrix(c(1, 0.3, 0.2,  0.3, 1.5, 0.4,  0.2, 0.4, 2), R, R)
  Vinv <- solve(V)

  Y_missing <- matrix(FALSE, N, R)
  Y_missing[5, 2] <- TRUE  # only obs 5, outcome 2 missing
  Y[Y_missing] <- 0

  # Manual scalar conditional normal:
  # E[Y2 | Y1, Y3] = mu2 - Vinv[2,2]^-1 * (Vinv[2,1]*(Y1-mu1) + Vinv[2,3]*(Y3-mu3))
  cov_22 <- 1 / Vinv[2, 2]
  resid <- c(Y[5, 1] - mu[5, 1], Y[5, 3] - mu[5, 3])
  vinv_mo <- Vinv[2, c(1, 3)]
  expected <- mu[5, 2] - cov_22 * sum(vinv_mo * resid)

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  expect_equal(res$Y[5, 2], expected, tolerance = 1e-14)
})


# ============================================================================
# 13. All entries missing: imputed value = mu
# ============================================================================
test_that("All entries missing: imputed value equals mu", {
  set.seed(555)
  N <- 10; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(FALSE, N, R)
  Y_missing[4, ] <- TRUE  # obs 4: all outcomes missing
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  # When all entries are missing, no conditioning data -> imputed = mu
  expect_equal(res$Y[4, ], mu[4, ], tolerance = 1e-14)
})


# ============================================================================
# 14. Precomputed cache reuse with different mu/Y
# ============================================================================
test_that("Same cache works correctly with different mu and Y", {
  set.seed(666)
  N <- 30; R <- 3
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.3, N, R)

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)

  # First call
  Y1 <- matrix(rnorm(N * R), N, R); Y1[Y_missing] <- 0
  mu1 <- matrix(rnorm(N * R), N, R)
  res1 <- mvsusieR:::impute_missing_Y(Y1, mu1, Vinv, miss_info,
                                       pattern_cache, version = "R")
  naive1 <- naive_impute_missing_Y(Y1, mu1, Vinv, Y_missing)

  # Second call with different data
  Y2 <- matrix(rnorm(N * R), N, R); Y2[Y_missing] <- 0
  mu2 <- matrix(rnorm(N * R), N, R)
  res2 <- mvsusieR:::impute_missing_Y(Y2, mu2, Vinv, miss_info,
                                       pattern_cache, version = "R")
  naive2 <- naive_impute_missing_Y(Y2, mu2, Vinv, Y_missing)

  expect_equal(res1$Y, naive1$Y, tolerance = 1e-14)
  expect_equal(res2$Y, naive2$Y, tolerance = 1e-14)
})


# ============================================================================
# 15. Large R: numerical stability
# ============================================================================
test_that("Pattern-cached imputation works for large R without numerical issues", {
  set.seed(777)
  N <- 100; R <- 20
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.5), N, R)
  # Create well-conditioned PD matrix
  V <- crossprod(matrix(rnorm(R * R), R, R)) / R + diag(R) * 2
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.15, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pattern_cache <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info,
                                      pattern_cache, version = "R")

  # No NaN or Inf values
  expect_true(all(is.finite(res$Y)))
  expect_true(all(is.finite(res$Y_cov)))
  expect_true(is.finite(res$sum_neg_ent_Y_miss))

  # Y_cov should be PSD
  evals <- eigen(res$Y_cov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals >= -1e-10))
})
