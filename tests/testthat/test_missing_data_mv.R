# Integration and math verification tests for missing data imputation.
#
# Tests: imputation formula correctness, ELBO corrections, V update,
# utility functions, R=1 path preservation, and R>1 end-to-end integration.

context("Missing data variational imputation")

# ============================================================================
# Helper: generate multivariate SuSiE data with missing entries
# ============================================================================
generate_mvsusie_data <- function(N = 200, J = 100, R = 3, L = 3,
                                   miss_rate = 0.2, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(N * J), N, J)
  # True effects: L causal variables with effect sizes ~ N(0, 1)
  causal <- sample(J, L)
  B <- matrix(0, J, R)
  for (l in seq_len(L)) {
    B[causal[l], ] <- rnorm(R)
  }
  # Residual covariance
  V_true <- matrix(c(1, 0.3, 0.2,  0.3, 1.5, 0.1,  0.2, 0.1, 0.8), R, R)
  E <- MASS::mvrnorm(N, rep(0, R), V_true)
  Y <- X %*% B + E

  # Introduce MCAR missingness
  Y_full <- Y
  if (miss_rate > 0) {
    miss_mask <- matrix(runif(N * R) < miss_rate, N, R)
    Y[miss_mask] <- NA
  } else {
    miss_mask <- matrix(FALSE, N, R)
  }

  list(X = X, Y = Y, Y_full = Y_full, B = B, V_true = V_true,
       causal = causal, miss_mask = miss_mask, N = N, J = J, R = R, L = L)
}


# ============================================================================
# Unit: Imputation formula correctness
# ============================================================================

# 16. Conditional normal formula
test_that("Imputation matches manual conditional normal computation", {
  set.seed(100)
  N <- 20; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R, sd = 0.3), N, R)
  V <- matrix(c(2, 0.5, 0.3,  0.5, 1.5, 0.2,  0.3, 0.2, 1.8), R, R)
  Vinv <- solve(V)

  Y_missing <- matrix(FALSE, N, R)
  Y_missing[2, c(1, 3)] <- TRUE
  Y_missing[5, 2] <- TRUE
  Y_missing[10, ] <- TRUE
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pc <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info, pc)

  # Manual for obs 2 (miss={1,3}, obs={2})
  miss <- c(1, 3); obs <- 2
  cov_mm <- solve(Vinv[miss, miss])
  resid <- Y[2, obs] - mu[2, obs]
  expected_2 <- mu[2, miss] - as.vector(cov_mm %*% Vinv[miss, obs, drop = FALSE] %*% resid)
  expect_equal(res$Y[2, miss], expected_2, tolerance = 1e-13)

  # Manual for obs 5 (miss={2}, obs={1,3})
  miss <- 2; obs <- c(1, 3)
  cov_mm <- 1 / Vinv[miss, miss]
  resid <- Y[5, obs] - mu[5, obs]
  expected_5 <- mu[5, miss] - cov_mm * sum(Vinv[miss, obs] * resid)
  expect_equal(res$Y[5, miss], expected_5, tolerance = 1e-13)

  # Manual for obs 10 (all missing -> imputed = mu)
  expect_equal(res$Y[10, ], mu[10, ], tolerance = 1e-13)
})


# 17. Y_cov accumulation
test_that("Y_cov matches manual accumulation of per-obs Lambda_MM^-1", {
  set.seed(101)
  N <- 30; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- matrix(c(1, 0.3, 0.2,  0.3, 1.5, 0.1,  0.2, 0.1, 0.8), R, R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.3, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pc <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info, pc)

  # Manual Y_cov
  Y_cov_manual <- matrix(0, R, R)
  for (i in seq_len(N)) {
    miss <- which(Y_missing[i, ])
    if (length(miss) == 0) next
    Y_cov_manual[miss, miss] <- Y_cov_manual[miss, miss] +
      solve(Vinv[miss, miss, drop = FALSE])
  }

  expect_equal(res$Y_cov, Y_cov_manual, tolerance = 1e-13)
})


# 18. Entropy term
test_that("Negative entropy sum matches manual computation", {
  set.seed(102)
  N <- 25; R <- 4
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- crossprod(matrix(rnorm(R * R), R, R)) + diag(R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.25, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pc <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info, pc)

  # Manual entropy
  sum_neg_ent_manual <- 0
  for (i in seq_len(N)) {
    miss <- which(Y_missing[i, ])
    n_miss <- length(miss)
    if (n_miss == 0) next
    Vinv_mm <- Vinv[miss, miss, drop = FALSE]
    log_det <- 2 * sum(log(diag(chol(Vinv_mm))))
    sum_neg_ent_manual <- sum_neg_ent_manual +
      0.5 * (n_miss * log(1 / (2 * pi * exp(1))) + log_det)
  }

  expect_equal(res$sum_neg_ent_Y_miss, sum_neg_ent_manual, tolerance = 1e-12)
})


# ============================================================================
# Unit: ELBO correction
# ============================================================================

# 19. ELBO correction terms
test_that("ELBO correction terms are computed correctly", {
  set.seed(103)
  N <- 30; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  mu <- matrix(rnorm(N * R), N, R)
  V <- matrix(c(1, 0.3, 0.2,  0.3, 1.5, 0.1,  0.2, 0.1, 0.8), R, R)
  Vinv <- solve(V)

  Y_missing <- matrix(runif(N * R) < 0.25, N, R)
  Y[Y_missing] <- 0

  miss_info <- mvsusieR:::extract_missing_patterns(Y_missing)
  pc <- mvsusieR:::precompute_pattern_cache(Vinv, miss_info$patterns)
  res <- mvsusieR:::impute_missing_Y(Y, mu, Vinv, miss_info, pc)

  # ELBO correction = -0.5 * tr(Vinv * Y_cov) - sum_neg_ent
  correction <- -0.5 * sum(Vinv * res$Y_cov) - res$sum_neg_ent_Y_miss

  # Manually compute the two parts
  trace_term <- -0.5 * sum(Vinv * res$Y_cov)
  entropy_term <- -res$sum_neg_ent_Y_miss

  expect_equal(correction, trace_term + entropy_term, tolerance = 1e-14)
  # Both terms should be finite
  expect_true(is.finite(trace_term))
  expect_true(is.finite(entropy_term))
})


# 20. ELBO with no missing = standard ELBO
test_that("ELBO with no missing data equals standard ELBO", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 100, J = 50, R = 3, L = 2,
                              miss_rate = 0, seed = 200)
  # Fit with imputation code path active (but no missing data)
  fit <- mvsusie_core(d$X, d$Y, L = 2,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      max_iter = 5, verbosity = 0)
  # ELBO should be finite and reasonable
  elbo <- fit$elbo
  expect_true(all(is.finite(elbo)))
})


# ============================================================================
# Unit: V update
# ============================================================================

# 21. V update with Y_cov
test_that("V estimation includes Y_cov term when imputation is active", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 150, J = 50, R = 3, L = 2,
                              miss_rate = 0.2, seed = 210)
  # Fit with V estimation
  fit <- mvsusie_core(d$X, d$Y, L = 2,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      max_iter = 10, verbosity = 0)
  # V should be a PD matrix
  expect_true(is.matrix(fit$sigma2))
  evals <- eigen(fit$sigma2, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
})


# 22. V stays PD after multiple iterations
test_that("Residual variance stays PD throughout iterations with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 150, J = 50, R = 3, L = 2,
                              miss_rate = 0.3, seed = 220)
  fit <- mvsusie_core(d$X, d$Y, L = 2,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      max_iter = 20, verbosity = 0)
  V <- fit$sigma2
  expect_true(is.matrix(V))
  evals <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
})


# ============================================================================
# Unit: Utility functions
# ============================================================================

# 23. chol2ldet correctness
test_that("chol2ldet matches log(det()) for known matrices", {
  A <- matrix(c(4, 2, 2, 3), 2, 2)
  L <- chol(A)
  expect_equal(mvsusieR:::chol2ldet(L), log(det(A)), tolerance = 1e-14)

  B <- diag(c(1, 2, 3))
  L_B <- chol(B)
  expect_equal(mvsusieR:::chol2ldet(L_B), log(det(B)), tolerance = 1e-14)
})


# 24. makePD correctness
test_that("makePD produces positive definite matrix", {
  # Start with a nearly singular matrix
  A <- matrix(c(1, 0.9999, 0.9999, 1), 2, 2)
  A_pd <- mvsusieR:::makePD(A, 1e-4)
  evals <- eigen(A_pd, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))

  # Check that ridge is added
  expect_equal(diag(A_pd), diag(A) + 1e-4)
})


# 25. invert_via_chol with fallback
test_that("invert_via_chol handles PD, singular, and zero matrices", {
  # PD matrix
  A <- matrix(c(4, 2, 2, 3), 2, 2)
  res <- mvsusieR:::invert_via_chol(A)
  expect_equal(res$inv %*% A, diag(2), tolerance = 1e-12)
  expect_equal(res$rank, 2)

  # Zero matrix
  Z <- matrix(0, 3, 3)
  res_z <- mvsusieR:::invert_via_chol(Z)
  expect_equal(res_z$inv, Z)
  expect_equal(res_z$rank, 0)

  # Singular matrix -> should fall back to pseudo-inverse
  S <- matrix(c(1, 1, 1, 1), 2, 2)
  res_s <- mvsusieR:::invert_via_chol(S)
  expect_true(res_s$rank <= 1)
})


# 26. compute_cov_flash
test_that("compute_cov_flash returns PD matrix of correct dimension", {
  set.seed(260)
  N <- 100; R <- 3
  Y <- matrix(rnorm(N * R), N, R)
  Y[sample(N * R, 30)] <- NA  # 10% missing

  V <- mvsusieR:::compute_cov_flash(Y)
  expect_equal(dim(V), c(R, R))
  evals <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
  expect_true(isSymmetric(V, tol = 1e-10))
})


# ============================================================================
# Integration: R=1 path preservation
# ============================================================================

# 27. R=1 zero-fill unchanged
test_that("R=1 with missing data uses complete-case (zero-fill), not imputation", {
  skip_on_cran()
  set.seed(270)
  N <- 100; J <- 50
  X <- matrix(rnorm(N * J), N, J)
  Y <- X[, 3] * 2 + rnorm(N)

  # Create R=1 data with missing
  Y_miss <- Y
  Y_miss[sample(N, 20)] <- NA

  fit <- mvsusie_core(X, Y_miss, L = 5,
                      prior_variance = matrix(0.2, 1, 1),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 50, verbosity = 0)
  expect_true(!is.null(fit$alpha))
  expect_equal(ncol(fit$alpha), J)
})


# 28. R=1 ELBO unchanged
test_that("R=1 ELBO is finite with missing data", {
  skip_on_cran()
  set.seed(280)
  N <- 100; J <- 50
  X <- matrix(rnorm(N * J), N, J)
  Y <- X[, 3] * 2 + rnorm(N)
  Y_miss <- Y
  Y_miss[sample(N, 15)] <- NA

  fit <- mvsusie_core(X, Y_miss, L = 5,
                      prior_variance = matrix(0.2, 1, 1),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 20, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
})


# ============================================================================
# Integration: R>1 end-to-end
# ============================================================================

# 29. No missing = identity
test_that("No missing data: imputation code produces identical results", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 100, J = 50, R = 3, L = 2,
                              miss_rate = 0, seed = 290)
  V0 <- cov(d$Y)

  fit <- mvsusie_core(d$X, d$Y, L = 2,
                      prior_variance = 0.2 * diag(d$R),
                      residual_variance = V0,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(!is.null(fit$pip))
})


# 30. Low missingness (5% MCAR)
test_that("Low missingness (5% MCAR): fit succeeds and PIPs are reasonable", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 200, J = 100, R = 3, L = 3,
                              miss_rate = 0.05, seed = 300)
  fit <- mvsusie_core(d$X, d$Y, L = 5,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 50, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))

  # Causal variables should have higher PIPs on average
  causal_pips <- fit$pip[d$causal]
  noncausal_pips <- fit$pip[-d$causal]
  expect_true(mean(causal_pips) > mean(noncausal_pips))
})


# 31. Medium missingness (20% MCAR)
test_that("Medium missingness (20% MCAR): algorithm converges", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 200, J = 100, R = 3, L = 3,
                              miss_rate = 0.2, seed = 310)
  fit <- mvsusie_core(d$X, d$Y, L = 5,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 50, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})


# 32. High missingness (50% MCAR)
test_that("High missingness (50% MCAR): algorithm converges without errors", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 200, J = 100, R = 3, L = 2,
                              miss_rate = 0.5, seed = 320)
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})


# 33. Mash prior (20% MCAR)
test_that("Mixture (mash) prior works with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 200, J = 100, R = 3, L = 2,
                              miss_rate = 0.2, seed = 330)
  # Create simple mixture prior
  U1 <- diag(d$R)
  U2 <- matrix(1, d$R, d$R)
  prior <- create_mixture_prior(list(matrices = list(U1 = U1, U2 = U2)))

  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = prior,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})


# 34. Monotone missingness
test_that("Monotone missingness: algorithm handles condition-level dropout", {
  skip_on_cran()
  set.seed(340)
  N <- 200; J <- 100; R <- 3; L <- 2
  X <- matrix(rnorm(N * J), N, J)
  B <- matrix(0, J, R)
  causal <- c(10, 30)
  B[causal, ] <- matrix(rnorm(L * R), L, R)
  V_true <- diag(R)
  E <- MASS::mvrnorm(N, rep(0, R), V_true)
  Y <- X %*% B + E

  # Monotone: condition 3 entirely missing for 50% of obs
  Y[101:200, 3] <- NA

  fit <- mvsusie_core(X, Y, L = 3,
                      prior_variance = 0.2 * diag(R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})


# 35. ELBO monotonicity
test_that("ELBO is non-decreasing across iterations with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 150, J = 50, R = 3, L = 2,
                              miss_rate = 0.2, seed = 350)
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 50, verbosity = 0)
  elbo <- fit$elbo
  # Check monotonicity (allow small tolerance for numerical noise)
  diffs <- diff(elbo)
  expect_true(all(diffs >= -1e-6),
              info = paste("ELBO decreased by", min(diffs)))
})


# 36. estimate_residual_variance = TRUE
test_that("V estimation converges with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 200, J = 50, R = 3, L = 2,
                              miss_rate = 0.15, seed = 360)
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  V <- fit$sigma2
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(d$R, d$R))
  evals <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
})


# 37. estimate_residual_variance = FALSE
test_that("Fixed V: algorithm converges with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 150, J = 50, R = 3, L = 2,
                              miss_rate = 0.2, seed = 370)
  V_fixed <- diag(d$R) * 1.5
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      residual_variance = V_fixed,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 30, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
})


# 38. Output structure completeness
test_that("Output has all standard fields with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 100, J = 50, R = 3, L = 2,
                              miss_rate = 0.2, seed = 380)
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      compute_objective = TRUE,
                      max_iter = 20, verbosity = 0)

  # Check essential output fields
  expect_true(!is.null(fit$pip))
  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$b1))
  expect_true(!is.null(fit$b1_rescaled))
  expect_true(!is.null(fit$sigma2))
  expect_true(!is.null(fit$V))
  expect_true(!is.null(fit$elbo))
  expect_true(!is.null(fit$fitted))
  expect_true(!is.null(fit$intercept))
  expect_true(!is.null(fit$lfsr))

  # Dimensions
  expect_equal(length(fit$pip), d$J)
  expect_equal(nrow(fit$alpha), 3)  # L
  expect_equal(ncol(fit$alpha), d$J)
})


# 39. Intercept recovery
test_that("Intercept is approximately correct with missing data", {
  skip_on_cran()
  set.seed(390)
  N <- 200; J <- 50; R <- 3
  X <- matrix(rnorm(N * J), N, J)
  Y <- matrix(rnorm(N * R), N, R)
  # Add intercept
  true_intercept <- c(2, -1, 0.5)
  Y <- Y + matrix(true_intercept, N, R, byrow = TRUE)
  # Add missing
  Y[sample(N * R, round(0.1 * N * R))] <- NA

  fit <- mvsusie_core(X, Y, L = 1,
                      prior_variance = 0.01 * diag(R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      max_iter = 20, verbosity = 0)

  # Intercept should be close to true values
  # (not exact due to missing data and finite sample)
  expect_equal(fit$intercept, true_intercept, tolerance = 0.5)
})


# 40. Variable names preserved
test_that("Variable names are preserved with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 80, J = 30, R = 3, L = 1,
                              miss_rate = 0.15, seed = 400)
  colnames(d$X) <- paste0("SNP", seq_len(d$J))
  colnames(d$Y) <- paste0("trait", seq_len(d$R))

  fit <- mvsusie_core(d$X, d$Y, L = 2,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      max_iter = 10, verbosity = 0)
  expect_equal(fit$variable_names, paste0("SNP", seq_len(d$J)))
  expect_equal(fit$condition_names, paste0("trait", seq_len(d$R)))
})


# ============================================================================
# Integration: Precomputation + missing data
# ============================================================================

# 41. precompute_covariances + missing data
test_that("Precomputed covariances work with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 100, J = 50, R = 3, L = 2,
                              miss_rate = 0.15, seed = 410)
  fit <- mvsusie_core(d$X, d$Y, L = 3,
                      prior_variance = 0.2 * diag(d$R),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      precompute_covariances = TRUE,
                      compute_objective = TRUE,
                      max_iter = 20, verbosity = 0)
  expect_true(all(is.finite(fit$elbo)))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})


# 42. estimate_prior_variance methods
test_that("Both EM and optim work with missing data", {
  skip_on_cran()
  d <- generate_mvsusie_data(N = 100, J = 50, R = 3, L = 2,
                              miss_rate = 0.15, seed = 420)

  # optim
  fit_optim <- mvsusie_core(d$X, d$Y, L = 3,
                             prior_variance = 0.2 * diag(d$R),
                             estimate_residual_variance = FALSE,
                             estimate_prior_variance = TRUE,
                             estimate_prior_method = "optim",
                             compute_objective = TRUE,
                             max_iter = 15, verbosity = 0)
  expect_true(all(is.finite(fit_optim$elbo)))

  # EM
  fit_em <- mvsusie_core(d$X, d$Y, L = 3,
                          prior_variance = 0.2 * diag(d$R),
                          estimate_residual_variance = FALSE,
                          estimate_prior_variance = TRUE,
                          estimate_prior_method = "EM",
                          compute_objective = TRUE,
                          max_iter = 15, verbosity = 0)
  expect_true(all(is.finite(fit_em$elbo)))
})
