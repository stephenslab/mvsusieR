# =============================================================================
# Integration tests: mvsusie with approximate/exact methods
# =============================================================================

test_that("approximate method runs and converges with 20% MCAR", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  expect_true(!is.null(fit$alpha))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(max(fit$pip) > 0.1)
})

test_that("exact method runs and converges with 20% MCAR", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "exact",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  expect_true(!is.null(fit$alpha))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(max(fit$pip) > 0.1)
})

test_that("approximate method equals standard when no data is missing", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2)
  V <- cov(sim$y)
  fit_standard <- mvsusie(sim$X, sim$y, L = 5,
                          prior_variance = sim$V,
                          missing_y_method = "impute",
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          residual_variance = V,
                          max_iter = 5, verbose = FALSE)
  fit_approx <- mvsusie(sim$X, sim$y, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "approximate",
                         estimate_prior_variance = FALSE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V,
                         max_iter = 5, verbose = FALSE)
  expect_equal(fit_standard$alpha, fit_approx$alpha)
  expect_equal(fit_standard$pip, fit_approx$pip)
})

test_that("approximate and impute PIPs are correlated with 20% MCAR", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2, y_missing = 0.2)
  V <- cov(sim$y)
  fit_impute <- mvsusie(sim$X, sim$y_missing, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "impute",
                         estimate_prior_variance = TRUE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V,
                         max_iter = 10, tol = 1e-3, verbose = FALSE)
  fit_approx <- mvsusie(sim$X, sim$y_missing, L = 5,
                          prior_variance = sim$V,
                          missing_y_method = "approximate",
                          estimate_prior_variance = TRUE,
                          estimate_residual_variance = FALSE,
                          residual_variance = V,
                          max_iter = 10, tol = 1e-3, verbose = FALSE)
  pip_cor <- cor(fit_impute$pip, fit_approx$pip)
  expect_true(pip_cor > 0.5)
})

test_that("exact method has monotonically non-decreasing ELBO", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "exact",
                 estimate_prior_variance = FALSE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  elbo <- fit$elbo[!is.na(fit$elbo)]
  if (length(elbo) > 1) {
    diffs <- diff(elbo)
    expect_true(all(diffs > -1e-3))
  }
})

test_that("approximate correctly identifies null effects in pure noise", {
  set.seed(42)
  n <- 60; p <- 30; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  for (i in 1:n) Y[i, runif(r) <= 0.2] <- NA
  fit <- mvsusie(X, Y, L = 3,
                 prior_variance = diag(r),
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = diag(r),
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  expect_true(all(fit$V < 1.0))
})

test_that("missing_y_method is ignored for R=1", {
  sim <- simulate_univariate()
  fit1 <- mvsusie(sim$X, sim$y, L = 5,
                  missing_y_method = "approximate",
                  estimate_prior_variance = FALSE,
                  max_iter = 3, verbose = FALSE)
  fit2 <- mvsusie(sim$X, sim$y, L = 5,
                  missing_y_method = "impute",
                  estimate_prior_variance = FALSE,
                  max_iter = 3, verbose = FALSE)
  expect_equal(fit1$alpha, fit2$alpha)
})

test_that("approximate/exact allow estimate_residual_variance = TRUE", {
  sim <- simulate_multivariate(n = 40, p = 20, r = 2, y_missing = 0.2)
  # est_rv=TRUE is now allowed for approximate/exact methods
  # (the update formula uses expected sufficient statistics, not ELBO)
  fit <- suppressWarnings(
    mvsusie(sim$X, sim$y_missing, L = 3,
            prior_variance = sim$V,
            missing_y_method = "approximate",
            estimate_residual_variance = TRUE,
            max_iter = 3, verbose = FALSE)
  )
  expect_true(!is.null(fit$alpha))
  expect_true(all(is.finite(fit$pip)))
})

test_that("approximate equals exact when residual variance is diagonal", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 2, y_missing = 0.2)
  V_diag <- diag(diag(cov(sim$y)))
  fit_approx <- mvsusie(sim$X, sim$y_missing, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "approximate",
                         estimate_prior_variance = FALSE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V_diag,
                         max_iter = 5, verbose = FALSE)
  fit_exact <- mvsusie(sim$X, sim$y_missing, L = 5,
                        prior_variance = sim$V,
                        missing_y_method = "exact",
                        estimate_prior_variance = FALSE,
                        estimate_residual_variance = FALSE,
                        residual_variance = V_diag,
                        max_iter = 5, verbose = FALSE)
  expect_equal(fit_approx$alpha, fit_exact$alpha, tolerance = 0.1)
  expect_equal(fit_approx$pip, fit_exact$pip, tolerance = 0.1)
})

test_that("approximate converges for R=3 with mixture prior", {
  sim <- simulate_multivariate(n = 60, p = 30, r = 3, y_missing = 0.15)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = c(0.5, 1, 2))
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = mash_init,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("approximate method handles non-overlapping sample groups", {
  set.seed(42)
  n <- 60; p <- 30; r <- 4
  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, r)
  beta[1:3, ] <- matrix(rnorm(12), 3, 4)
  Y <- X %*% beta + matrix(rnorm(n * r), n, r)
  Y[1:30, 3:4] <- NA
  Y[31:60, 1:2] <- NA
  V <- diag(r)
  fit <- mvsusie(X, Y, L = 5,
                 prior_variance = V,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = V,
                 max_iter = 10, tol = 1e-3, verbose = FALSE)
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# =============================================================================
# Math verification: precomputed svs_inv vs naive per-sample loop
# =============================================================================

test_that("precomputed svs_inv matches naive per-sample computation", {
  set.seed(123)
  n <- 30; p <- 10; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y_missing <- matrix(FALSE, n, r)
  Y_missing[1:5, 1] <- TRUE
  Y_missing[6:10, 2] <- TRUE
  Y[Y_missing] <- NA
  V <- diag(r)

  data <- mvsusieR:::create_mvsusie_data(X, Y, center = TRUE, scale = TRUE,
                                          missing_y_method = "approximate")
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
  my <- data$miss3d

  for (j in 1:p) {
    naive_svs_inv_j <- matrix(0, r, r)
    for (i in 1:n) {
      xij <- my$X_3d[i, j, ]
      Vinv_i <- my$Vinv[[my$pattern_assign[i]]]
      naive_svs_inv_j <- naive_svs_inv_j + t(Vinv_i * xij) * xij
    }
    expect_equal(data$svs_inv[[j]], naive_svs_inv_j, tolerance = 1e-10,
                 info = paste("svs_inv mismatch at j =", j))
  }
})

test_that("compute_XtR_3d matches naive per-variable computation", {
  set.seed(123)
  n <- 30; p <- 10; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y_missing <- matrix(FALSE, n, r)
  Y_missing[1:5, 1] <- TRUE
  Y_missing[6:10, 2] <- TRUE
  Y[Y_missing] <- NA
  V <- diag(r)

  data <- mvsusieR:::create_mvsusie_data(X, Y, center = TRUE, scale = TRUE,
                                          missing_y_method = "approximate")
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
  my <- data$miss3d

  R_mat <- matrix(rnorm(n * r), n, r)
  R_mat[!my$Y_non_missing] <- 0

  XtR_opt <- mvsusieR:::compute_XtR_3d(data, R_mat)

  Vinv <- my$Vinv
  VinvR <- matrix(0, n, r)
  for (i in 1:n) {
    VinvR[i, ] <- as.numeric(Vinv[[my$pattern_assign[i]]] %*% R_mat[i, ])
  }
  XtR_naive <- t(sapply(1:p, function(j) {
    colSums(my$X_3d[, j, ] * VinvR)
  }))

  expect_equal(XtR_opt, XtR_naive, tolerance = 1e-10)
})

# =============================================================================
# C++ vs R: missing data functions (parametric over method and R)
# =============================================================================

# Helper: create test data for a given method and R dimension
make_test_data <- function(n, p, r, miss_frac, method, V = NULL) {
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y[sample(n * r, round(miss_frac * n * r))] <- NA
  data <- mvsusieR:::create_mvsusie_data(X, Y, missing_y_method = method)
  if (is.null(V)) {
    M <- matrix(rnorm(r * r), r, r)
    V <- crossprod(M) / r + diag(r) * 0.5
  }
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
  data
}

test_that("C++ VinvR, XtR, Xb, betahat match R across methods and R dims", {
  for (r in c(2, 5)) {
    for (method in c("approximate", "exact")) {
      set.seed(100 + r)
      data <- make_test_data(40, 20, r, 0.2, method)
      label <- paste(method, "R =", r)

      # VinvR
      mat <- matrix(rnorm(data$n * r), data$n, r)
      expect_equal(mvsusieR:::compute_VinvR_3d(data, mat),
                   mvsusieR:::compute_VinvR_3d_R(data, mat),
                   tolerance = 1e-10, info = paste("VinvR", label))

      # XtR
      R_mat <- matrix(rnorm(data$n * r), data$n, r)
      expect_equal(mvsusieR:::compute_XtR_3d(data, R_mat),
                   mvsusieR:::compute_XtR_3d_R(data, R_mat),
                   tolerance = 1e-10, info = paste("XtR", label))

      # Xb
      b <- matrix(rnorm(data$p * r), data$p, r)
      expect_equal(mvsusieR:::compute_Xb_3d(data, b),
                   mvsusieR:::compute_Xb_3d_R(data, b),
                   tolerance = 1e-10, info = paste("Xb", label))

      # betahat
      XtR <- mvsusieR:::compute_XtR_3d(data, R_mat)
      expect_equal(mvsusieR:::compute_betahat_3d(data, XtR),
                   mvsusieR:::compute_betahat_3d_R(data, XtR),
                   tolerance = 1e-10, info = paste("betahat", label))
    }
  }
})

test_that("C++ VinvR matches R with near-singular V", {
  set.seed(103)
  data <- make_test_data(40, 15, 5, 0.15, "approximate")
  V_base <- crossprod(matrix(rnorm(25), 5, 5)) / 5 + diag(5) * 0.5
  eig <- eigen(V_base, symmetric = TRUE)
  eig$values[5] <- 1e-6
  V_nearsing <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V_nearsing)

  mat <- matrix(rnorm(data$n * 5), data$n, 5)
  expect_equal(mvsusieR:::compute_VinvR_3d(data, mat),
               mvsusieR:::compute_VinvR_3d_R(data, mat),
               tolerance = 1e-8)
})

# =============================================================================
# C++ vs R: svs_inv naive loop verification (approximate + exact)
# =============================================================================

test_that("C++ svs_inv matches naive loop for approximate and exact", {
  for (method in c("approximate", "exact")) {
    set.seed(500)
    n <- 30; p <- 10; r <- 3
    X <- matrix(rnorm(n * p), n, p)
    Y <- matrix(rnorm(n * r), n, r)
    Y[1:8, 1] <- NA; Y[9:15, 2] <- NA; Y[16:20, 3] <- NA
    data <- mvsusieR:::create_mvsusie_data(X, Y, missing_y_method = method)
    V <- diag(r) * 1.5
    data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
    my <- data$miss3d

    if (method == "approximate") {
      for (j in seq_len(p)) {
        naive <- matrix(0, r, r)
        for (i in seq_len(n)) {
          xij <- my$X_3d[i, j, ]
          naive <- naive + t(my$Vinv[[my$pattern_assign[i]]] * xij) * xij
        }
        expect_equal(data$svs_inv[[j]], naive, tolerance = 1e-10,
                     info = paste("approx svs_inv j =", j))
      }
    } else {
      Vinvsum <- Reduce("+", lapply(seq_len(nrow(my$pattern)), function(k) {
        my$Vinv[[k]] * my$n_k[k]
      }))
      for (j in seq_len(p)) {
        A1 <- matrix(0, r, r)
        A2 <- matrix(0, r, r)
        for (i in seq_len(n)) {
          xij <- my$X_3d[i, j, ]
          Vinv_i <- my$Vinv[[my$pattern_assign[i]]]
          A1 <- A1 + t(Vinv_i * xij) * xij
          A2 <- A2 + t(t(Vinv_i) * xij)
        }
        Xbar_j <- my$Xbar[j, , ]
        if (!is.matrix(Xbar_j)) Xbar_j <- matrix(Xbar_j, r, r)
        naive <- A1 - crossprod(Xbar_j, A2) - crossprod(A2, Xbar_j) +
          crossprod(Xbar_j, Vinvsum %*% Xbar_j)
        expect_equal(data$svs_inv[[j]], naive, tolerance = 1e-9,
                     info = paste("exact svs_inv j =", j))
      }
    }
  }
})

# =============================================================================
# C++ vs R: Xbar from precomputed sums
# =============================================================================

test_that("C++ compute_Xbar_from_sums matches naive N-loop", {
  set.seed(42)
  sim <- simulate_multivariate(n = 40, p = 15, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "exact")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  my <- data$miss3d
  J <- data$p; N <- data$n; R_dim <- data$R

  Xbar_naive <- array(0, dim = c(J, R_dim, R_dim))
  for (j in seq_len(J)) {
    inner <- matrix(0, R_dim, R_dim)
    for (i in seq_len(N)) {
      inner <- inner + t(t(my$Vinv[[my$pattern_assign[i]]]) * my$X_3d[i, j, ])
    }
    Xbar_naive[j, , ] <- my$Vinvsuminv %*% inner
  }

  expect_equal(my$Xbar, Xbar_naive, tolerance = 1e-10)
})

# =============================================================================
# C++ vs R: eigendecomposition precomputation
# =============================================================================

# Helper: generate random SPD svs and PSD V_structure
make_eigen_data <- function(J, R, K, near_singular = FALSE) {
  svs <- lapply(seq_len(J), function(j) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) + diag(R) * 0.5
  })
  V_structure <- lapply(seq_len(K), function(k) {
    M <- matrix(rnorm(R * R), R, R)
    U <- crossprod(M) / R
    if (near_singular) {
      eig <- eigen(U, symmetric = TRUE)
      eig$values[R] <- 1e-10
      U <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    }
    U
  })
  list(svs = svs, V_structure = V_structure)
}

compare_eigen_caches <- function(cache_cpp, cache_r, K, J, tol = 1e-10) {
  expect_equal(as.numeric(cache_cpp$log_det_svs),
               as.numeric(cache_r$log_det_svs), tolerance = tol)
  for (k in seq_len(K)) {
    check_j <- unique(c(1, ceiling(J / 2), J))
    for (j in check_j) {
      eig_cpp <- sort(cache_cpp$components[[k]]$eigenvalues[, j], decreasing = TRUE)
      eig_r   <- sort(cache_r$components[[k]]$eigenvalues[, j], decreasing = TRUE)
      expect_equal(eig_cpp, eig_r, tolerance = tol,
                   info = paste("eigenvalue k =", k, "j =", j))
    }
  }
}

compare_posterior <- function(post_cpp, post_r, tol_mean = 1e-8,
                               tol_neg = 1e-6) {
  expect_equal(post_cpp$post_mean, post_r$post_mean, tolerance = tol_mean)
  expect_equal(post_cpp$post_mean2, post_r$post_mean2, tolerance = tol_mean)
  expect_equal(post_cpp$post_neg, post_r$post_neg, tolerance = tol_neg)
  expect_equal(post_cpp$post_zero, post_r$post_zero, tolerance = 1e-10)
  expect_equal(as.numeric(post_cpp$prior_scale_em_update),
               as.numeric(post_r$prior_scale_em_update), tolerance = tol_mean)
}

test_that("C++ eigen cache matches R for non-common-cov", {
  for (params in list(
    list(J = 20, R = 2, K = 3, ns = FALSE),
    list(J = 20, R = 5, K = 6, ns = FALSE),
    list(J = 15, R = 4, K = 3, ns = TRUE)
  )) {
    set.seed(601 + params$R)
    d <- make_eigen_data(params$J, params$R, params$K, params$ns)
    cache_cpp <- mvsusieR:::precompute_eigen_cache(d$svs, d$V_structure, FALSE)
    cache_r   <- mvsusieR:::precompute_eigen_cache_R(d$svs, d$V_structure, FALSE)
    tol <- if (params$ns) 1e-8 else 1e-10
    compare_eigen_caches(cache_cpp, cache_r, params$K, params$J, tol)
  }
})

# =============================================================================
# C++ vs R: loglik non-common-cov
# =============================================================================

test_that("C++ loglik_non_common matches R across parameter regimes", {
  for (params in list(
    list(J = 20, R = 2, K = 3, V = 1.5),
    list(J = 20, R = 5, K = 6, V = 1.0),
    list(J = 15, R = 3, K = 2, V = 1e-6),
    list(J = 15, R = 3, K = 4, V = 100.0)
  )) {
    set.seed(700 + params$R)
    d <- make_eigen_data(params$J, params$R, params$K)
    betahat <- matrix(rnorm(params$J * params$R), params$J, params$R)
    cache <- mvsusieR:::precompute_eigen_cache(d$svs, d$V_structure, FALSE)

    llik_cpp <- mvsusieR:::loglik_non_common_rcpp(
      betahat, params$V, cache$log_det_svs, cache$components)
    llik_r <- mvsusieR:::loglik_precomputed_R(betahat, params$V, cache)

    expect_equal(llik_cpp, llik_r, tolerance = 1e-9,
                 info = paste("R =", params$R, "V =", params$V))
  }
})

# =============================================================================
# C++ vs R: posterior non-common-cov
# =============================================================================

test_that("C++ posterior_non_common matches R across parameter regimes", {
  for (params in list(
    list(J = 20, R = 2, K = 3, V = 1.5, em = FALSE),
    list(J = 20, R = 5, K = 6, V = 1.0, em = FALSE),
    list(J = 15, R = 3, K = 4, V = 2.0, em = TRUE),
    list(J = 15, R = 3, K = 1, V = 1.5, em = FALSE)
  )) {
    set.seed(800 + params$R)
    d <- make_eigen_data(params$J, params$R, params$K)
    betahat <- matrix(rnorm(params$J * params$R), params$J, params$R)
    cache <- mvsusieR:::precompute_eigen_cache(d$svs, d$V_structure, FALSE)

    pi_raw <- matrix(runif(params$J * (params$K + 1)), params$J, params$K + 1)
    pi_V_post <- pi_raw / rowSums(pi_raw)

    em_wt <- matrix(0, 0, 0)
    em_wt_r <- NULL
    if (params$em) {
      em_raw <- matrix(runif((params$K + 1) * params$J), params$K + 1, params$J)
      em_wt <- em_raw / rowSums(em_raw)
      em_wt_r <- em_wt
    }

    post_cpp <- mvsusieR:::posterior_non_common_rcpp(
      betahat, params$V, cache$components, pi_V_post, em_wt)
    post_r <- mvsusieR:::posterior_precomputed_R(
      betahat, params$V, cache, pi_V_post, em_wt_r)

    compare_posterior(post_cpp, post_r)
  }
})

# =============================================================================
# C++ vs R: accumulate_post_mean2_common
# =============================================================================

test_that("C++ accumulate_post_mean2_common matches R loop", {
  for (R in c(2, 4)) {
    set.seed(900 + R)
    J <- 30
    post_mean2 <- array(rnorm(J * R * R), c(J, R, R))
    M_k <- matrix(rnorm(J * R), J, R)
    C_k <- crossprod(matrix(rnorm(R * R), R, R)) / R
    w_k <- runif(J)

    pm2_r <- post_mean2
    for (j in seq_len(J)) {
      pm2_r[j, , ] <- pm2_r[j, , ] + w_k[j] * (C_k + tcrossprod(M_k[j, ]))
    }
    pm2_cpp <- mvsusieR:::accumulate_post_mean2_common_rcpp(
      post_mean2, M_k, C_k, w_k)
    expect_equal(pm2_cpp, pm2_r, tolerance = 1e-10,
                 info = paste("R =", R))
  }
})

# =============================================================================
# Full pipeline: create_data -> set_V -> all ops (approx + exact)
# =============================================================================

test_that("Full pipeline C++ vs R for R=5 approximate and exact", {
  for (method in c("approximate", "exact")) {
    set.seed(1001)
    data <- make_test_data(40, 20, 5, 0.25, method)

    mat <- matrix(rnorm(data$n * 5), data$n, 5)
    expect_equal(mvsusieR:::compute_VinvR_3d(data, mat),
                 mvsusieR:::compute_VinvR_3d_R(data, mat),
                 tolerance = 1e-10, info = method)

    R_mat <- matrix(rnorm(data$n * 5), data$n, 5)
    expect_equal(mvsusieR:::compute_XtR_3d(data, R_mat),
                 mvsusieR:::compute_XtR_3d_R(data, R_mat),
                 tolerance = 1e-10, info = method)

    b <- matrix(rnorm(data$p * 5), data$p, 5)
    expect_equal(mvsusieR:::compute_Xb_3d(data, b),
                 mvsusieR:::compute_Xb_3d_R(data, b),
                 tolerance = 1e-10, info = method)

    XtR <- mvsusieR:::compute_XtR_3d(data, R_mat)
    expect_equal(mvsusieR:::compute_betahat_3d(data, XtR),
                 mvsusieR:::compute_betahat_3d_R(data, XtR),
                 tolerance = 1e-10, info = method)
  }
})

# =============================================================================
# Full eigen pipeline: cache -> loglik -> posterior end-to-end
# =============================================================================

test_that("Full eigen pipeline C++ vs R end-to-end", {
  set.seed(1101)
  J <- 20; R <- 5; K <- 6
  d <- make_eigen_data(J, R, K)
  betahat <- matrix(rnorm(J * R), J, R)
  V_scalar <- 1.5

  cache_cpp <- mvsusieR:::precompute_eigen_cache(d$svs, d$V_structure, FALSE)
  cache_r   <- mvsusieR:::precompute_eigen_cache_R(d$svs, d$V_structure, FALSE)

  llik_cpp <- mvsusieR:::loglik_non_common_rcpp(
    betahat, V_scalar, cache_cpp$log_det_svs, cache_cpp$components)
  llik_r <- mvsusieR:::loglik_precomputed_R(betahat, V_scalar, cache_r)
  expect_equal(llik_cpp, llik_r, tolerance = 1e-9)

  pi_raw <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_raw / rowSums(pi_raw)

  post_cpp <- mvsusieR:::posterior_non_common_rcpp(
    betahat, V_scalar, cache_cpp$components, pi_V_post, matrix(0, 0, 0))
  post_r <- mvsusieR:::posterior_precomputed_R(
    betahat, V_scalar, cache_r, pi_V_post)

  compare_posterior(post_cpp, post_r)
})

# =============================================================================
# Edge case: block missingness patterns
# =============================================================================

test_that("C++ matches R with block missingness patterns R=5", {
  set.seed(1301)
  n <- 60; p <- 20; r <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y[1:15, 1:2] <- NA
  Y[16:30, 3:4] <- NA
  Y[31:45, 5] <- NA

  for (method in c("approximate", "exact")) {
    data <- mvsusieR:::create_mvsusie_data(X, Y, missing_y_method = method)
    M <- matrix(rnorm(r * r), r, r)
    V <- crossprod(M) / r + diag(r) * 0.5
    data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

    mat <- matrix(rnorm(data$n * r), data$n, r)
    expect_equal(mvsusieR:::compute_VinvR_3d(data, mat),
                 mvsusieR:::compute_VinvR_3d_R(data, mat),
                 tolerance = 1e-10, info = method)

    R_mat <- matrix(rnorm(data$n * r), data$n, r)
    expect_equal(mvsusieR:::compute_XtR_3d(data, R_mat),
                 mvsusieR:::compute_XtR_3d_R(data, R_mat),
                 tolerance = 1e-10, info = method)

    b <- matrix(rnorm(data$p * r), data$p, r)
    expect_equal(mvsusieR:::compute_Xb_3d(data, b),
                 mvsusieR:::compute_Xb_3d_R(data, b),
                 tolerance = 1e-10, info = method)
  }
})
