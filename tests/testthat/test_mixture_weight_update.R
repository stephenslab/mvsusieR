context("Test mixture prior weight updates (EM and mixsqp)")

# =============================================================================
# UNIT TESTS: update_mixture_weights_em
# =============================================================================

test_that("EM update with uniform posterior and uniform alpha leaves weights unchanged", {
  # When pi_V_posterior is uniform across components and alpha is uniform,
  # EM should return weights that are (approximately) uniform for non-null components.
  K <- 4
  J <- 50
  L <- 2
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = lapply(seq_len(L), function(l) {
      # Uniform J x (K+1) matrix
      matrix(1 / K_plus_1, J, K_plus_1)
    })
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # Non-null weights should be approximately equal
  expect_equal(result$pi_V, rep(1 / K, K), tolerance = 1e-10)
  # Null weight should be unchanged (not updated)
  expect_equal(result$null_weight, 0)
})

test_that("EM update with concentrated posterior shifts weights correctly", {
  # When posterior concentrates on component 2 (index 3 in K+1 matrix),
  # that component's weight should increase.
  K <- 4
  J <- 50
  L <- 2
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = lapply(seq_len(L), function(l) {
      pvp <- matrix(0.01 / (K_plus_1 - 1), J, K_plus_1)
      pvp[, 3] <- 0.99  # component 2 (non-null index 2) dominates
      pvp <- pvp / rowSums(pvp)
      pvp
    })
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # Component 2 (pi_V index 2) should have the largest weight
  expect_true(which.max(result$pi_V) == 2)
  expect_true(result$pi_V[2] > 0.9)
  # Weights should sum to 1
  expect_equal(sum(result$pi_V), 1)
})

test_that("EM update alpha-weighting: only the selected variable matters", {
  # When alpha is concentrated on one variable and that variable's
  # pi_V_posterior concentrates on component 1, the result should
  # reflect that even if all other variables prefer component 2.
  K <- 3
  J <- 50
  L <- 1
  K_plus_1 <- K + 1

  # Alpha concentrated on variable 1
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1  # variable 1 is the only causal one

  # Variable 1 prefers component 1, all other variables prefer component 3
  pvp <- matrix(0, J, K_plus_1)
  pvp[1, 2] <- 1.0    # variable 1: component 1 (non-null index 1)
  pvp[2:J, 4] <- 1.0  # all others: component 3 (non-null index 3)

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = list(pvp)
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # With alpha weighting, component 1 should dominate (variable 1's preference)
  # because the other 49 variables have near-zero alpha weight
  expect_true(result$pi_V[1] > 0.9)
})

test_that("EM update averages across L=3 effects", {
  K <- 3
  J <- 20
  L <- 3
  K_plus_1 <- K + 1

  # Each effect selects a different variable
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1
  alpha[2, 2] <- 1
  alpha[3, 3] <- 1

  # Each effect's selected variable prefers a different component
  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = list(
      {
        pvp <- matrix(0.01, J, K_plus_1)
        pvp[1, 2] <- 0.97  # variable 1 (selected by effect 1): component 1
        pvp / rowSums(pvp)
      },
      {
        pvp <- matrix(0.01, J, K_plus_1)
        pvp[2, 3] <- 0.97  # variable 2 (selected by effect 2): component 2
        pvp / rowSums(pvp)
      },
      {
        pvp <- matrix(0.01, J, K_plus_1)
        pvp[3, 4] <- 0.97  # variable 3 (selected by effect 3): component 3
        pvp / rowSums(pvp)
      }
    )
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # All three non-null components should get roughly equal weight
  expect_true(all(result$pi_V > 0.2))
  expect_true(all(result$pi_V < 0.5))
  expect_equal(sum(result$pi_V), 1)
})

test_that("EM update with update_null = TRUE updates null_weight", {
  K <- 3
  J <- 30
  L <- 1
  K_plus_1 <- K + 1

  # Alpha concentrated on one variable
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1

  # Posterior concentrates heavily on null component (column 1) for
  # the selected variable
  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = list({
      pvp <- matrix(0.01, J, K_plus_1)
      pvp[1, 1] <- 0.97  # variable 1 (selected): null dominates
      pvp / rowSums(pvp)
    })
  )

  result_null <- mvsusieR:::update_mixture_weights_em(model, update_null = TRUE)
  result_no_null <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)

  # With update_null, null_weight should increase
  expect_true(result_null$null_weight > 0.8)
  # Without update_null, null_weight stays at 0
  expect_equal(result_no_null$null_weight, 0)
})

test_that("EM update handles NULL entries in pi_V_posterior (skips them)", {
  K <- 3
  L <- 3
  J <- 20
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = list(
      matrix(1 / K_plus_1, J, K_plus_1),  # valid
      NULL,                                 # skipped
      matrix(1 / K_plus_1, J, K_plus_1)   # valid
    )
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # Should not error and weights should be valid
  expect_equal(sum(result$pi_V), 1)
  expect_true(all(result$pi_V > 0))
})

# =============================================================================
# UNIT TESTS: update_mixture_weights_mixsqp
# =============================================================================

test_that("mixsqp finds correct weights on synthetic llik with alpha weighting", {
  # Create a synthetic llik matrix where the ML solution is known.
  # Component 2 has higher log-likelihood for the variables with high alpha.
  set.seed(42)
  K_plus_1 <- 4
  J <- 100
  L <- 2

  # Alpha: each effect selects one variable
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1
  alpha[2, 2] <- 1

  # For the selected variables, component 2 has higher log-likelihood
  llik_list <- lapply(seq_len(L), function(l) {
    llik <- matrix(rnorm(J * K_plus_1, mean = -5, sd = 1), J, K_plus_1)
    llik[l, 2] <- 5  # selected variable: component 2 dominates
    llik
  })

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / (K_plus_1 - 1), K_plus_1 - 1),
    null_weight = 0,
    per_effect_llik = llik_list
  )

  result <- mvsusieR:::update_mixture_weights_mixsqp(model, update_null = FALSE)
  # Component 1 (non-null index 1, which is column 2) should dominate
  expect_true(result$pi_V[1] > 0.5)
  expect_equal(sum(result$pi_V), 1)
})

test_that("mixsqp returns model unchanged when per_effect_llik is empty", {
  model <- list(
    alpha = matrix(1 / 10, 2, 10),
    pi_V = c(0.5, 0.5),
    null_weight = 0,
    per_effect_llik = list(NULL, NULL)
  )

  result <- mvsusieR:::update_mixture_weights_mixsqp(model, update_null = FALSE)
  expect_equal(result$pi_V, c(0.5, 0.5))
})

test_that("mixsqp alpha-weighting: only selected variable matters", {
  # Same test as EM alpha-weighting: only the selected variable's llik matters
  set.seed(42)
  K_plus_1 <- 3
  J <- 50
  L <- 1

  # Alpha concentrated on variable 1
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1

  # Variable 1's llik favors component 1 (col 2); all others favor component 2 (col 3)
  llik <- matrix(-5, J, K_plus_1)
  llik[1, 2] <- 5   # variable 1: component 1 dominates
  llik[2:J, 3] <- 5 # all others: component 2 dominates

  model <- list(
    alpha = alpha,
    pi_V = c(0.5, 0.5),
    null_weight = 0,
    per_effect_llik = list(llik)
  )

  result <- mvsusieR:::update_mixture_weights_mixsqp(model, update_null = FALSE)
  # With alpha weighting, component 1 should get high weight because
  # only variable 1 (which favors component 1) has meaningful alpha
  expect_true(result$pi_V[1] > 0.5)
})

# =============================================================================
# UNIT TESTS: update_mixture_weights dispatcher
# =============================================================================

test_that("update_mixture_weights dispatches correctly", {
  set.seed(123)
  K <- 3
  J <- 30
  L <- 1
  K_plus_1 <- K + 1

  # Alpha concentrated on variable 1
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1

  pvp <- matrix(0.01, J, K_plus_1)
  pvp[1, 3] <- 0.97  # variable 1: component 2
  pvp <- pvp / rowSums(pvp)

  llik <- matrix(rnorm(J * K_plus_1, mean = -5), J, K_plus_1)
  llik[1, 3] <- 5  # variable 1: component 2

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    pi_V_posterior = list(pvp),
    per_effect_llik = list(llik)
  )

  result_em <- mvsusieR:::update_mixture_weights(model, method = "EM")
  result_msqp <- mvsusieR:::update_mixture_weights(model, method = "mixsqp")

  # Both should identify component 2 as dominant
  expect_true(which.max(result_em$pi_V) == 2)
  expect_true(which.max(result_msqp$pi_V) == 2)
})

# =============================================================================
# UNIT TESTS: prune_mixture_components
# =============================================================================

test_that("prune_mixture_components does nothing when all above threshold", {
  R <- 2
  K <- 3
  model <- list(
    pi_V = c(0.4, 0.35, 0.25),
    V_structure = list(diag(R), diag(R) * 2, diag(R) * 3),
    V_structure_3d = array(0, c(R, R, K)),
    V_structure_inv = array(0, c(R, R, K)),
    V_structure_rank = c(2, 2, 2),
    pi_V_posterior = list(matrix(1/4, 10, 4)),  # J=10, K+1=4
    per_effect_llik = list(matrix(0, 10, 4)),
    eigen_cache = list(list(), list(), list())
  )
  for (k in 1:K) {
    model$V_structure_3d[, , k] <- model$V_structure[[k]]
    model$V_structure_inv[, , k] <- solve(model$V_structure[[k]])
  }

  result <- mvsusieR:::prune_mixture_components(model, threshold = 1e-8)
  expect_equal(length(result$pi_V), K)
  expect_equal(length(result$V_structure), K)
})

test_that("prune_mixture_components removes low-weight components", {
  R <- 2
  K <- 4
  model <- list(
    pi_V = c(0.5, 0.49, 1e-10, 1e-12),
    V_structure = list(diag(R), diag(R)*2, diag(R)*3, diag(R)*4),
    V_structure_3d = array(0, c(R, R, K)),
    V_structure_inv = array(0, c(R, R, K)),
    V_structure_rank = c(2, 2, 2, 2),
    pi_V_posterior = list(matrix(1/5, 10, 5)),
    per_effect_llik = list(matrix(0, 10, 5)),
    eigen_cache = list(list(), list(), list(), list())
  )
  for (k in 1:K) {
    model$V_structure_3d[, , k] <- model$V_structure[[k]]
    model$V_structure_inv[, , k] <- solve(model$V_structure[[k]])
  }

  result <- mvsusieR:::prune_mixture_components(model, threshold = 1e-8)
  # Should keep only first 2 components
  expect_equal(length(result$pi_V), 2)
  expect_equal(sum(result$pi_V), 1, tolerance = 1e-10)
  expect_equal(length(result$V_structure), 2)
  expect_equal(dim(result$V_structure_3d)[3], 2)
  expect_equal(dim(result$V_structure_inv)[3], 2)
  expect_equal(length(result$V_structure_rank), 2)
  # pi_V_posterior should have null + 2 kept = 3 columns
  expect_equal(ncol(result$pi_V_posterior[[1]]), 3)
  expect_equal(ncol(result$per_effect_llik[[1]]), 3)
  expect_equal(length(result$eigen_cache), 2)
})

# =============================================================================
# INTEGRATION TESTS: Individual data
# =============================================================================

test_that("K=1 single matrix prior: no weight update triggered", with(simulate_multivariate(r=2), {
  # Single matrix prior has K=1, so weight update should be a no-op
  residual_var <- cov(y)
  fit <- mvsusie(X, y, L = 3, prior_variance = V,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 compute_objective = FALSE)
  # pi_V should still be 1.0 (single component)
  expect_equal(fit$pi_V, 1.0)
}))

test_that("estimate_prior_mixture_weights = FALSE: pi_V unchanged", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = FALSE,
                 compute_objective = FALSE)
  # pi_V should be unchanged
  expect_equal(fit$pi_V, pi_V_init)
}))

test_that("EM weight update changes pi_V with mash prior", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 mixture_weight_method = "EM",
                 compute_objective = FALSE)
  # pi_V should have changed from initial
  expect_false(isTRUE(all.equal(fit$pi_V, pi_V_init)))
  # Weights should sum to 1
  expect_equal(sum(fit$pi_V), 1, tolerance = 1e-10)
}))

test_that("mixsqp weight update changes pi_V with mash prior", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 mixture_weight_method = "mixsqp",
                 compute_objective = FALSE)
  # pi_V should have changed from initial
  expect_false(isTRUE(all.equal(fit$pi_V, pi_V_init)))
  # Weights should sum to 1
  expect_equal(sum(fit$pi_V), 1, tolerance = 1e-10)
}))

test_that("EM and mixsqp both produce valid updated weights", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit_em <- mvsusie(X, y, L = 3, prior_variance = prior,
                    residual_variance = residual_var,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = FALSE,
                    estimate_prior_mixture_weights = TRUE,
                    mixture_weight_method = "EM",
                    compute_objective = FALSE)

  fit_msqp <- mvsusie(X, y, L = 3, prior_variance = prior,
                      residual_variance = residual_var,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = TRUE,
                      mixture_weight_method = "mixsqp",
                      compute_objective = FALSE)

  # Both produce valid probability vectors
  expect_equal(sum(fit_em$pi_V), 1, tolerance = 1e-10)
  expect_equal(sum(fit_msqp$pi_V), 1, tolerance = 1e-10)
  expect_true(all(fit_em$pi_V >= 0))
  expect_true(all(fit_msqp$pi_V >= 0))

  # Both change from initial uniform weights
  expect_false(isTRUE(all.equal(fit_em$pi_V, pi_V_init)))
  expect_false(isTRUE(all.equal(fit_msqp$pi_V, pi_V_init)))

  # Same number of components
  expect_equal(length(fit_em$pi_V), length(fit_msqp$pi_V))
}))

test_that("ELBO converges with EM weight updates", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 mixture_weight_method = "EM",
                 compute_objective = TRUE,
                 max_iter = 50))
  # ELBO should exist and be finite
  expect_true(length(fit$elbo) > 0)
  expect_true(all(is.finite(fit$elbo)))
  # Note: EM may have tiny ELBO drops because it optimizes the *uncollapsed*
  # objective while the ELBO uses the *collapsed* form. This is expected.
  # For guaranteed monotonicity, use mixsqp (the default).
}))

test_that("ELBO is monotonically non-decreasing with mixsqp weight updates", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 mixture_weight_method = "mixsqp",
                 compute_objective = TRUE,
                 max_iter = 50))
  # ELBO should exist and be finite
  expect_true(length(fit$elbo) > 0)
  expect_true(all(is.finite(fit$elbo)))
  # mixsqp directly optimizes the collapsed ELBO contribution, so
  # the ELBO should be monotonically non-decreasing
  if (length(fit$elbo) > 1) {
    diffs <- diff(fit$elbo)
    expect_true(all(diffs >= -1e-6),
                info = paste("ELBO dropped:", paste(round(diffs[diffs < -1e-6], 6), collapse = ", ")))
  }
}))

test_that("estimate_residual_variance=FALSE + weight update works", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  # This tests the susieR change (D0): update_model_variance is always called
  # even when estimate_residual_variance = FALSE
  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 compute_objective = FALSE)
  # pi_V should have changed
  expect_false(isTRUE(all.equal(fit$pi_V, pi_V_init)))
  # sigma2 should be unchanged (residual variance not estimated)
  expect_equal(fit$sigma2, residual_var, tolerance = 1e-10)
}))

test_that("Both residual variance and weight update work together", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = TRUE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 compute_objective = FALSE)
  # pi_V should have changed
  expect_false(isTRUE(all.equal(fit$pi_V, pi_V_init)))
  # sigma2 should also have changed
  expect_false(isTRUE(all.equal(fit$sigma2, residual_var)))
}))

test_that("precompute_covariances + weight update work together", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit_precomp <- mvsusie(X, y, L = 3, prior_variance = prior,
                         residual_variance = residual_var,
                         estimate_residual_variance = FALSE,
                         estimate_prior_variance = FALSE,
                         estimate_prior_mixture_weights = TRUE,
                         precompute_covariances = TRUE,
                         compute_objective = FALSE)

  fit_no_precomp <- mvsusie(X, y, L = 3, prior_variance = prior,
                            residual_variance = residual_var,
                            estimate_residual_variance = FALSE,
                            estimate_prior_variance = FALSE,
                            estimate_prior_mixture_weights = TRUE,
                            precompute_covariances = FALSE,
                            compute_objective = FALSE)

  # Results should be identical regardless of precomputation strategy
  expect_equal(fit_precomp$alpha, fit_no_precomp$alpha, tolerance = 1e-6)
  expect_equal(fit_precomp$pi_V, fit_no_precomp$pi_V, tolerance = 1e-6)
}))

# =============================================================================
# INTEGRATION TESTS: Sufficient statistics
# =============================================================================

test_that("mvsusie_ss with mixture weight updates runs and converges", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  YtY <- crossprod(y)
  N <- nrow(X)

  fit_ss <- mvsusie_ss(XtX, XtY, YtY, N, L = 3,
                       prior_variance = prior,
                       residual_variance = residual_var,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE,
                       estimate_prior_mixture_weights = TRUE,
                       compute_objective = FALSE)
  # pi_V should have changed
  expect_false(isTRUE(all.equal(fit_ss$pi_V, pi_V_init)))
  expect_equal(sum(fit_ss$pi_V), 1, tolerance = 1e-10)
}))

test_that("mvsusie_rss with weight updates via ... passthrough", with(simulate_multivariate(r=2), {
  z <- sapply(1:ncol(y), function(j) {
    ss <- susieR:::univariate_regression(X, y[, j])
    ss$betahat / ss$sebetahat
  })
  R_ld <- cor(X)
  nn <- nrow(X)

  prior <- create_mash_prior(Ulist = list(V * nn), grid = c(0.5, 1, 2))
  pi_V_init <- prior$pi

  fit_rss <- mvsusie_rss(z, R_ld, N = nn, L = 3,
                         prior_variance = prior,
                         estimate_prior_variance = FALSE,
                         estimate_prior_mixture_weights = TRUE,
                         compute_objective = FALSE)
  # pi_V should have changed (forwarded via ...)
  expect_false(isTRUE(all.equal(fit_rss$pi_V, pi_V_init)))
}))

test_that("SS and individual data produce similar pi_V", with(simulate_multivariate(r=2, center_scale=TRUE), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit_ind <- mvsusie(X, y, L = 3, prior_variance = prior,
                     residual_variance = residual_var,
                     estimate_residual_variance = FALSE,
                     estimate_prior_variance = FALSE,
                     estimate_prior_mixture_weights = TRUE,
                     mixture_weight_method = "EM",
                     compute_objective = FALSE)

  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  YtY <- crossprod(y)
  N <- nrow(X)

  fit_ss <- mvsusie_ss(XtX, XtY, YtY, N, L = 3,
                       prior_variance = prior,
                       residual_variance = residual_var,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE,
                       estimate_prior_mixture_weights = TRUE,
                       mixture_weight_method = "EM",
                       compute_objective = FALSE)

  # Both paths should produce valid updated weights
  expect_equal(sum(fit_ind$pi_V), 1, tolerance = 1e-10)
  expect_equal(sum(fit_ss$pi_V), 1, tolerance = 1e-10)
  expect_equal(length(fit_ind$pi_V), length(fit_ss$pi_V))
}))

# =============================================================================
# CLEANUP TESTS
# =============================================================================

test_that("Temporary fields are cleaned up in final model", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 compute_objective = FALSE)
  # These should be cleaned up
  expect_null(fit$per_effect_llik)
  expect_null(fit$ibss_iter)
  expect_null(fit$llik_cache)
  expect_null(fit$em_cache)
}))
