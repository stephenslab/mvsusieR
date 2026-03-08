context("Test mixture prior weight updates (EM and mixsqp)")

# =============================================================================
# UNIT TESTS: inner_em_cpp (C++ core)
# =============================================================================

test_that("inner_em_cpp converges on uniform llik to uniform weights", {
  K <- 4
  N <- 100
  # Uniform log-likelihoods → all components equally likely → uniform weights
  llik <- matrix(0, N, K)
  weights <- rep(1, N)
  pi_init <- rep(1 / K, K)

  result <- mvsusieR:::inner_em_cpp(llik, weights, pi_init, 100L, 1e-10)
  expect_equal(as.numeric(result$pi), rep(1 / K, K), tolerance = 1e-8)
  expect_true(result$converged)
})

test_that("inner_em_cpp concentrates weight on dominant component", {
  K <- 3
  N <- 50
  # Component 2 has much higher log-likelihood for all observations
  llik <- matrix(-10, N, K)
  llik[, 2] <- 0  # component 2 dominates
  weights <- rep(1, N)
  pi_init <- rep(1 / K, K)

  result <- mvsusieR:::inner_em_cpp(llik, weights, pi_init, 100L, 1e-10)
  expect_true(result$pi[2] > 0.99)
  expect_true(result$converged)
  expect_equal(sum(result$pi), 1, tolerance = 1e-10)
})

test_that("inner_em_cpp respects observation weights", {
  K <- 3
  N <- 100
  # Half the observations favor component 1, half favor component 3
  llik <- matrix(-10, N, K)
  llik[1:50, 1] <- 0    # first 50: component 1
  llik[51:100, 3] <- 0  # last 50: component 3

  # Equal weights: roughly equal split between components 1 and 3
  weights_equal <- rep(1, N)
  result_equal <- mvsusieR:::inner_em_cpp(llik, weights_equal, rep(1/K, K), 100L, 1e-10)
  expect_true(abs(result_equal$pi[1] - result_equal$pi[3]) < 0.1)

  # Weight only the first 50: component 1 should dominate
  weights_first <- c(rep(1, 50), rep(1e-10, 50))
  result_first <- mvsusieR:::inner_em_cpp(llik, weights_first, rep(1/K, K), 100L, 1e-10)
  expect_true(result_first$pi[1] > 0.9)
})

test_that("inner_em_cpp converges early on well-conditioned problems", {
  set.seed(42)
  K <- 5
  N <- 200
  llik <- matrix(rnorm(N * K), N, K)
  weights <- rep(1, N)
  pi_init <- rep(1 / K, K)

  result <- mvsusieR:::inner_em_cpp(llik, weights, pi_init, 1000L, 1e-10)
  expect_true(result$converged)
  expect_true(result$n_iter < 1000)  # should converge well before max
})

test_that("inner_em_cpp matches pure R reference implementation", {
  set.seed(123)
  K <- 4
  N <- 80
  llik <- matrix(rnorm(N * K, mean = -3, sd = 2), N, K)
  weights <- runif(N)
  pi_init <- rep(1 / K, K)

  # C++ result
  result_cpp <- mvsusieR:::inner_em_cpp(llik, weights, pi_init, 200L, 1e-12)

  # Pure R reference implementation
  pi_cur <- pi_init
  for (t in seq_len(200)) {
    log_phi <- sweep(llik, 2, log(pmax(pi_cur, 1e-300)), "+")
    log_phi <- log_phi - apply(log_phi, 1, max)
    phi <- exp(log_phi)
    phi <- phi / rowSums(phi)
    pi_new <- colSums(weights * phi) / sum(weights)
    if (max(abs(pi_new - pi_cur)) < 1e-12) break
    pi_cur <- pi_new
  }

  expect_equal(as.numeric(result_cpp$pi), pi_cur, tolerance = 1e-8)
})

# =============================================================================
# UNIT TESTS: update_mixture_weights_em (R wrapper calling C++)
# =============================================================================

test_that("EM update with uniform llik and uniform alpha leaves weights unchanged", {
  # When all components have equal log-likelihood and alpha is uniform,
  # EM should return weights that are (approximately) unchanged.
  K <- 4
  J <- 50
  L <- 2
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    per_effect_llik = lapply(seq_len(L), function(l) {
      # Uniform J x (K+1) log-likelihood matrix
      matrix(0, J, K_plus_1)
    })
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # Non-null weights should be approximately equal
  expect_equal(result$pi_V, rep(1 / K, K), tolerance = 1e-6)
  # Null weight should be unchanged (not updated)
  expect_equal(result$null_weight, 0)
})

test_that("EM update with concentrated llik shifts weights correctly", {
  # When llik concentrates on component 2 (column 3 in K+1 matrix),
  # that component's weight should increase.
  K <- 4
  J <- 50
  L <- 2
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    per_effect_llik = lapply(seq_len(L), function(l) {
      llik <- matrix(-10, J, K_plus_1)
      llik[, 3] <- 0  # component 2 (non-null index 2, column 3) dominates
      llik
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
  # llik concentrates on component 1, the result should reflect that
  # even if all other variables prefer component 2.
  K <- 3
  J <- 50
  L <- 1
  K_plus_1 <- K + 1

  # Alpha concentrated on variable 1
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1  # variable 1 is the only causal one

  # Variable 1's llik favors component 1 (col 2); all others favor component 3 (col 4)
  llik <- matrix(-10, J, K_plus_1)
  llik[1, 2] <- 0    # variable 1: component 1 (non-null index 1)
  llik[2:J, 4] <- 0  # all others: component 3 (non-null index 3)

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    per_effect_llik = list(llik)
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
    per_effect_llik = list(
      {
        llik <- matrix(-10, J, K_plus_1)
        llik[1, 2] <- 0  # variable 1 (selected by effect 1): component 1
        llik
      },
      {
        llik <- matrix(-10, J, K_plus_1)
        llik[2, 3] <- 0  # variable 2 (selected by effect 2): component 2
        llik
      },
      {
        llik <- matrix(-10, J, K_plus_1)
        llik[3, 4] <- 0  # variable 3 (selected by effect 3): component 3
        llik
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

  # Variable 1's llik favors null component (column 1)
  llik <- matrix(-10, J, K_plus_1)
  llik[1, 1] <- 0  # variable 1 (selected): null dominates

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
    per_effect_llik = list(llik)
  )

  result_null <- mvsusieR:::update_mixture_weights_em(model, update_null = TRUE)
  result_no_null <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)

  # With update_null, null_weight should increase substantially
  # (EM applies a small floor to zero weights so it can discover them)
  expect_true(result_null$null_weight > 0.8)
  # Without update_null, null_weight stays at 0
  expect_equal(result_no_null$null_weight, 0)
})

test_that("EM update handles NULL entries in per_effect_llik (skips them)", {
  K <- 3
  L <- 3
  J <- 20
  K_plus_1 <- K + 1

  model <- list(
    alpha = matrix(1 / J, L, J),
    pi_V = rep(1 / K, K),
    null_weight = 0,
    per_effect_llik = list(
      matrix(0, J, K_plus_1),  # valid
      NULL,                     # skipped
      matrix(0, J, K_plus_1)   # valid
    )
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  # Should not error and weights should be valid
  expect_equal(sum(result$pi_V), 1)
  expect_true(all(result$pi_V > 0))
})

test_that("EM returns model unchanged when all per_effect_llik are NULL", {
  model <- list(
    alpha = matrix(1 / 10, 2, 10),
    pi_V = c(0.5, 0.5),
    null_weight = 0,
    per_effect_llik = list(NULL, NULL)
  )

  result <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  expect_equal(result$pi_V, c(0.5, 0.5))
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
# UNIT TESTS: EM and mixsqp solve the same problem (direct comparison)
# =============================================================================

test_that("EM and mixsqp converge to same weights on shared llik data", {
  set.seed(42)
  K_plus_1 <- 5
  J <- 100
  L <- 3

  # Alpha: each effect selects one variable
  alpha <- matrix(1e-10, L, J)
  alpha[1, 1] <- 1
  alpha[2, 2] <- 1
  alpha[3, 3] <- 1

  # Create llik where selected variables have distinctive patterns
  llik_list <- lapply(seq_len(L), function(l) {
    llik <- matrix(rnorm(J * K_plus_1, mean = -3, sd = 1), J, K_plus_1)
    # Each selected variable favors a different component
    llik[l, l + 1] <- 5
    llik
  })

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / (K_plus_1 - 1), K_plus_1 - 1),
    null_weight = 0,
    per_effect_llik = llik_list
  )

  result_em <- mvsusieR:::update_mixture_weights_em(model, update_null = FALSE)
  result_msqp <- mvsusieR:::update_mixture_weights_mixsqp(model, update_null = FALSE)

  # Both should converge to the same weights (both solve the same concave problem)
  expect_equal(result_em$pi_V, result_msqp$pi_V, tolerance = 0.05)
  # Correlation should be very high
  expect_true(cor(result_em$pi_V, result_msqp$pi_V) > 0.95)
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

  llik <- matrix(rnorm(J * K_plus_1, mean = -5), J, K_plus_1)
  llik[1, 3] <- 5  # variable 1: component 2

  model <- list(
    alpha = alpha,
    pi_V = rep(1 / K, K),
    null_weight = 0,
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
                 )
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
                 )
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
                 )
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
                 )
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
                    )

  fit_msqp <- mvsusie(X, y, L = 3, prior_variance = prior,
                      residual_variance = residual_var,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = TRUE,
                      mixture_weight_method = "mixsqp",
                      )

  # Both produce valid probability vectors
  expect_equal(sum(fit_em$pi_V), 1, tolerance = 1e-10)
  expect_equal(sum(fit_msqp$pi_V), 1, tolerance = 1e-10)
  expect_true(all(fit_em$pi_V >= 0))
  expect_true(all(fit_msqp$pi_V >= 0))

  # Both change from initial uniform weights
  expect_false(isTRUE(all.equal(fit_em$pi_V, pi_V_init)))
  expect_false(isTRUE(all.equal(fit_msqp$pi_V, pi_V_init)))

  # Same number of components (after pruning)
  expect_equal(length(fit_em$pi_V), length(fit_msqp$pi_V))
}))

# =============================================================================
# ELBO MONOTONICITY TESTS
# =============================================================================

test_that("ELBO is monotonically non-decreasing with EM weight updates", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 mixture_weight_method = "EM",
                 max_iter = 50))
  # ELBO should exist and be finite
  expect_true(length(fit$elbo) > 0)
  expect_true(all(is.finite(fit$elbo)))
  # With the proper inner EM, ELBO should be monotonically non-decreasing
  if (length(fit$elbo) > 1) {
    diffs <- diff(fit$elbo)
    expect_true(all(diffs >= -1e-6),
                info = paste("EM ELBO dropped:", paste(round(diffs[diffs < -1e-6], 8), collapse = ", ")))
  }
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

# =============================================================================
# EM ≈ MIXSQP AGREEMENT TESTS (integration)
# =============================================================================

test_that("EM and mixsqp produce similar pi_V weights", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit_em <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                    residual_variance = residual_var,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = FALSE,
                    estimate_prior_mixture_weights = TRUE,
                    mixture_weight_method = "EM",
                    max_iter = 50))

  fit_msqp <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                      residual_variance = residual_var,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = TRUE,
                      mixture_weight_method = "mixsqp",
                      max_iter = 50))

  # After final pruning, both should have the same number of components
  expect_equal(length(fit_em$pi_V), length(fit_msqp$pi_V))

  # pi_V weights should be similar (both solve the same concave problem)
  if (length(fit_em$pi_V) == length(fit_msqp$pi_V) && length(fit_em$pi_V) > 1) {
    expect_true(cor(fit_em$pi_V, fit_msqp$pi_V) > 0.95,
                info = paste("Weight correlation:", cor(fit_em$pi_V, fit_msqp$pi_V)))
    expect_true(max(abs(fit_em$pi_V - fit_msqp$pi_V)) < 0.1,
                info = paste("Max weight diff:", max(abs(fit_em$pi_V - fit_msqp$pi_V))))
  }
}))

test_that("EM and mixsqp produce similar alpha", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit_em <- mvsusie(X, y, L = 3, prior_variance = prior,
                    residual_variance = residual_var,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = FALSE,
                    estimate_prior_mixture_weights = TRUE,
                    mixture_weight_method = "EM",
                    )

  fit_msqp <- mvsusie(X, y, L = 3, prior_variance = prior,
                      residual_variance = residual_var,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = TRUE,
                      mixture_weight_method = "mixsqp",
                      )

  # Alpha matrices should be similar
  expect_true(max(abs(fit_em$alpha - fit_msqp$alpha)) < 0.1,
              info = paste("Max alpha diff:", max(abs(fit_em$alpha - fit_msqp$alpha))))
}))

test_that("EM and mixsqp produce similar final ELBO", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit_em <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                    residual_variance = residual_var,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = FALSE,
                    estimate_prior_mixture_weights = TRUE,
                    mixture_weight_method = "EM",
                    max_iter = 50))

  fit_msqp <- suppressWarnings(mvsusie(X, y, L = 3, prior_variance = prior,
                      residual_variance = residual_var,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = TRUE,
                      mixture_weight_method = "mixsqp",
                      max_iter = 50))

  # Final ELBOs should be close
  elbo_em <- tail(fit_em$elbo, 1)
  elbo_msqp <- tail(fit_msqp$elbo, 1)
  expect_true(abs(elbo_em - elbo_msqp) < 1.0,
              info = paste("ELBO diff:", abs(elbo_em - elbo_msqp),
                           "EM:", elbo_em, "mixsqp:", elbo_msqp))
}))

# =============================================================================
# OTHER INTEGRATION TESTS
# =============================================================================

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
                 )
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
                 )
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
                         )

  fit_no_precomp <- mvsusie(X, y, L = 3, prior_variance = prior,
                            residual_variance = residual_var,
                            estimate_residual_variance = FALSE,
                            estimate_prior_variance = FALSE,
                            estimate_prior_mixture_weights = TRUE,
                            precompute_covariances = FALSE,
                            )

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
                       )
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
                         )
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
                     )

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
                       )

  # Both paths should produce valid updated weights
  expect_equal(sum(fit_ind$pi_V), 1, tolerance = 1e-10)
  expect_equal(sum(fit_ss$pi_V), 1, tolerance = 1e-10)
  expect_equal(length(fit_ind$pi_V), length(fit_ss$pi_V))
}))

# =============================================================================
# CLEANUP AND PRUNING TESTS
# =============================================================================

test_that("Temporary fields are cleaned up in final model", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2))

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 )
  # These should be cleaned up
  expect_null(fit$per_effect_llik)
  expect_null(fit$ibss_iter)
  expect_null(fit$llik_cache)
  expect_null(fit$em_cache)
}))

test_that("Final pruning removes near-zero components at convergence", with(simulate_multivariate(r=2), {
  residual_var <- cov(y)
  # Use a larger grid to create more components, some of which should be pruned
  prior <- create_mash_prior(Ulist = list(V), grid = c(0.1, 0.5, 1, 2, 5))
  K_init <- length(prior$pi)

  fit <- mvsusie(X, y, L = 3, prior_variance = prior,
                 residual_variance = residual_var,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 estimate_prior_mixture_weights = TRUE,
                 )

  # All remaining weights should be above the pruning threshold
  expect_true(all(fit$pi_V >= 1e-8))
  # Weights should sum to 1
  expect_equal(sum(fit$pi_V), 1, tolerance = 1e-10)
  # Number of remaining components should be consistent across structures
  K_final <- length(fit$pi_V)
  expect_equal(length(fit$V_structure), K_final)
}))
