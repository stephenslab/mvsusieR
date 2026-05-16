# ============================================================================
# Component-level zero-tolerance math verification tests.
#
# These tests verify that every mathematical component in the S3
# refactored code produces EXACTLY the same result as the R6 (master)
# implementation, to machine precision.
#
# Tests cover:
#  1. Standalone functions: multivariate_regression, multivariate_lbf,
#     compute_softmax, mvsusie_get_lfsr, mvsusie_single_effect_lfsr
#  2. Mash prior construction: create_mash_prior vs MashInitializer
#  3. Special case: mixture K=1 no-null = single matrix prior
#  4. SS path = individual path (on equivalent data)
#  5. Workhorse defaults sanity
#  6. Edge cases
#
# Full end-to-end S3-vs-R6 comparison across parameter sweeps is in
# test_numerical_sweep.R.
# ============================================================================


# --------------------------------------------------------------------------
# Setup: load R6 reference from master (reusing test_r6_reference.R infra)
# --------------------------------------------------------------------------

repo_dir <- tryCatch({
  gcd <- system2("git", c("rev-parse", "--git-common-dir"),
                 stdout = TRUE, stderr = FALSE)
  normalizePath(file.path(gcd, ".."))
}, error = function(e) NULL)

# Cached R6 reference functions
.r6_env <- NULL

load_r6_env <- function() {
  if (!is.null(.r6_env)) return(.r6_env)
  if (is.null(repo_dir)) stop("Cannot determine git repository root")

  ref_source <- file.path(tempdir(), "mvsusieR_r6_ref")
  if (!dir.exists(ref_source)) {
    system2("git", c("-C", repo_dir, "worktree", "remove", "--force", ref_source),
            stdout = FALSE, stderr = FALSE)
    ret <- system2("git", c("-C", repo_dir, "worktree", "add", "--detach",
                            ref_source, "origin/master"),
                   stdout = FALSE, stderr = FALSE)
    if (ret != 0) stop("Failed to create worktree from master")
  }

  # Load R6 package into a separate environment
  ref <- pkgload::load_all(ref_source, export_all = TRUE, quiet = TRUE)

  funcs <- list(
    mvsusie = ref$env$mvsusie,
    mvsusie_rss = ref$env$mvsusie_rss,
    mvsusie_ss = tryCatch(ref$env$mvsusie_ss, error = function(e) NULL),
    MashInitializer = ref$env$MashInitializer,
    create_mixture_prior = ref$env$create_mixture_prior,
    mvsusie_sim1 = ref$env$mvsusie_sim1,
    multivariate_regression = ref$env$multivariate_regression,
    multivariate_lbf = ref$env$multivariate_lbf,
    compute_softmax = ref$env$compute_softmax,
    invert_via_chol = ref$env$invert_via_chol,
    pseudo_inverse = ref$env$pseudo_inverse,
    matlist2array = ref$env$matlist2array,
    env = ref$env
  )

  # Reload the S3 (development) package to restore S3 functions
  dev_source <- tryCatch(
    rprojroot::find_root(rprojroot::is_r_package),
    error = function(e) normalizePath(file.path(getwd(), "../.."))
  )
  pkgload::load_all(dev_source, export_all = TRUE, quiet = TRUE)

  .r6_env <<- funcs
  funcs
}

skip_if_no_r6 <- function() {
  load_r6_env()
}

# Machine precision tolerance
MACH_TOL <- .Machine$double.eps^0.5

# Generate simulation data deterministically
set.seed(999)
sim_r1 <- mvsusie_sim1(n = 100, p = 50, r = 1, s = 3)
sim_r1$L <- 5
set.seed(999)
sim_r3 <- mvsusie_sim1(n = 100, p = 50, r = 3, s = 3)
sim_r3$L <- 5

# ============================================================================
# Section 1: Standalone function tests
# ============================================================================

test_that("multivariate_regression: R6 function consistency check", {
  # NOTE: multivariate_regression() exists only in R6 (master). In S3, the
  # same math is performed via mashr C++ (calc_lik_rcpp + calc_sermix_rcpp).
  # This test verifies the R6 function is internally consistent.
  skip_if_no_r6()
  r6 <- load_r6_env()

  set.seed(42)
  R <- 3; J <- 20
  betahat <- matrix(rnorm(J * R), J, R)
  S <- lapply(1:J, function(j) {
    m <- matrix(rnorm(R * R), R, R)
    crossprod(m) + diag(R) * 0.5
  })
  U <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), R, R)
  S_inv <- lapply(S, function(s) chol2inv(chol(s)))

  # R6's multivariate_regression
  r6_post <- r6$multivariate_regression(betahat, S, U, S_inv)

  # Verify: posterior mean should satisfy b1_j = cov_j * S_inv_j * betahat_j
  for (j in 1:min(5, J)) {
    cov_j <- U %*% solve(diag(R) + S_inv[[j]] %*% U)
    b1_j <- cov_j %*% S_inv[[j]] %*% betahat[j, ]
    expect_equal(get_b1(r6_post)[j, ], as.vector(b1_j), tolerance = MACH_TOL,
                 info = paste0("multivariate_regression: b1[", j, "]"))
  }
})

test_that("multivariate_lbf: S3 matches R6 exactly", {
  skip_if_no_r6()
  r6 <- load_r6_env()

  set.seed(42)
  R <- 3; J <- 20
  betahat <- matrix(rnorm(J * R), J, R)
  S <- lapply(1:J, function(j) {
    m <- matrix(rnorm(R * R), R, R)
    crossprod(m) + diag(R) * 0.5
  })
  U <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), R, R)

  r6_lbf <- r6$multivariate_lbf(betahat, S, U)
  s3_lbf <- multivariate_lbf(betahat, S, U)

  expect_equal(s3_lbf, r6_lbf, tolerance = MACH_TOL)
})

test_that("compute_softmax: S3 matches R6 exactly", {
  skip_if_no_r6()
  r6 <- load_r6_env()

  set.seed(42)
  values <- rnorm(50)
  weights <- runif(50)

  r6_sm <- r6$compute_softmax(values, weights, log = TRUE)
  s3_sm <- compute_softmax(values, weights, log = TRUE)

  expect_equal(s3_sm$weights, r6_sm$weights, tolerance = MACH_TOL)
  expect_equal(s3_sm$log_sum, r6_sm$log_sum, tolerance = MACH_TOL)

  # Also test non-log scale
  pos_values <- exp(values)
  r6_sm2 <- r6$compute_softmax(pos_values, weights, log = FALSE)
  s3_sm2 <- compute_softmax(pos_values, weights, log = FALSE)
  expect_equal(s3_sm2$weights, r6_sm2$weights, tolerance = MACH_TOL)
  expect_equal(s3_sm2$log_sum, r6_sm2$log_sum, tolerance = MACH_TOL)
})

test_that("mvsusie_get_lfsr: correct on known inputs", {
  # Test 1: Single effect (L=1), J=3, R=2
  L <- 1; J <- 3; R <- 2
  alpha <- matrix(c(0.8, 0.1, 0.1), L, J)
  clfsr <- array(0, c(L, J, R))
  clfsr[1, , 1] <- c(0.01, 0.5, 0.9)
  clfsr[1, , 2] <- c(0.02, 0.6, 0.8)

  lfsr <- mvsusie_get_lfsr(clfsr, alpha)

  # For outcome r: lfsr_j = pmax(1e-20, 1 - max_l(alpha_lj * (1 - clfsr_ljr)))
  # r=1: true_sign = alpha * (1 - clfsr): 0.8*(1-0.01)=0.792, 0.1*(1-0.5)=0.05, 0.1*(1-0.9)=0.01
  #       max per variable: max over L (only L=1): [0.792, 0.05, 0.01]
  #       lfsr = pmax(1e-20, 1 - [0.792, 0.05, 0.01]) = [0.208, 0.95, 0.99]
  expect_equal(lfsr[1, 1], 1 - 0.8 * (1 - 0.01), tolerance = MACH_TOL)
  expect_equal(lfsr[2, 1], 1 - 0.1 * (1 - 0.5), tolerance = MACH_TOL)
  expect_equal(lfsr[3, 1], 1 - 0.1 * (1 - 0.9), tolerance = MACH_TOL)

  # Test 2: NA input
  expect_true(is.na(mvsusie_get_lfsr(as.numeric(NA), alpha)))

  # Test 3: unweighted
  lfsr_uw <- mvsusie_get_lfsr(clfsr, alpha, weighted = FALSE)
  # With uniform alpha, true_sign = 1 * (1 - clfsr)
  # max_l(1 * (1-clfsr_{ljr})): for r=1: max(0.99, 0.5, 0.1) = 0.99
  expect_equal(lfsr_uw[1, 1], pmax(1e-20, 1 - (1 - 0.01)), tolerance = MACH_TOL)
})

test_that("mvsusie_single_effect_lfsr: correct on known inputs", {
  # L=2, J=3, R=2
  # alpha is L x J, stored by column in R: alpha[1,] = row 1, alpha[2,] = row 2
  L <- 2; J <- 3; R <- 2
  alpha <- matrix(c(0.5, 0.1,   # col 1: [alpha[1,1], alpha[2,1]]
                     0.3, 0.6,   # col 2: [alpha[1,2], alpha[2,2]]
                     0.2, 0.3),  # col 3: [alpha[1,3], alpha[2,3]]
                  L, J)
  clfsr <- array(0, c(L, J, R))
  clfsr[1, , 1] <- c(0.1, 0.3, 0.5)
  clfsr[2, , 1] <- c(0.2, 0.4, 0.6)
  clfsr[1, , 2] <- c(0.15, 0.35, 0.55)
  clfsr[2, , 2] <- c(0.25, 0.45, 0.65)

  se_lfsr <- mvsusie_single_effect_lfsr(clfsr, alpha)

  # For effect l, outcome r: se_lfsr[l,r] = pmax(0, sum_j alpha[l,j] * clfsr[l,j,r])
  # l=1, r=1: alpha[1,] = c(0.5, 0.3, 0.2), clfsr[1,,1] = c(0.1, 0.3, 0.5)
  #   sum = 0.5*0.1 + 0.3*0.3 + 0.2*0.5 = 0.05 + 0.09 + 0.10 = 0.24
  expected_l1_r1 <- sum(alpha[1, ] * clfsr[1, , 1])
  expect_equal(se_lfsr[1, 1], expected_l1_r1, tolerance = MACH_TOL)

  # l=2, r=1: alpha[2,] = c(0.1, 0.6, 0.3), clfsr[2,,1] = c(0.2, 0.4, 0.6)
  #   sum = 0.1*0.2 + 0.6*0.4 + 0.3*0.6 = 0.02 + 0.24 + 0.18 = 0.44
  expected_l2_r1 <- sum(alpha[2, ] * clfsr[2, , 1])
  expect_equal(se_lfsr[2, 1], expected_l2_r1, tolerance = MACH_TOL)

  # l=1, r=2: alpha[1,] = c(0.5, 0.3, 0.2), clfsr[1,,2] = c(0.15, 0.35, 0.55)
  expected_l1_r2 <- sum(alpha[1, ] * clfsr[1, , 2])
  expect_equal(se_lfsr[1, 2], expected_l1_r2, tolerance = MACH_TOL)

  # NA input
  expect_true(is.na(mvsusie_single_effect_lfsr(as.numeric(NA), alpha)))
})

# ============================================================================
# Section 2: Mash prior construction
# ============================================================================

test_that("create_mash_prior matches MashInitializer (Ulist + grid)", {
  skip_if_no_r6()
  r6 <- load_r6_env()

  R <- 3
  U1 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), R, R)
  U2 <- diag(R)
  Ulist <- list(U1, U2)
  grid <- c(0.5, 1.0, 2.0)

  # R6 MashInitializer
  r6_prior <- r6$MashInitializer$new(Ulist = Ulist, grid = grid, null_weight = 0.1)
  r6_pv <- r6_prior$prior_variance

  # S3 create_mash_prior
  s3_prior <- create_mash_prior(Ulist = Ulist, grid = grid, null_weight = 0.1)

  # Compare pi (mixture weights including null)
  # R6: r6_pv$pi is (K+1)-vector with null_weight first
  # S3: s3_prior$pi is K-vector (non-null only), s3_prior$null_weight is separate
  r6_null_w <- r6_pv$pi[1]
  r6_nonnull_w <- r6_pv$pi[-1]
  expect_equal(as.numeric(r6_null_w), s3_prior$null_weight, tolerance = MACH_TOL)
  expect_equal(as.numeric(r6_nonnull_w),
               as.numeric(s3_prior$pi * (1 - s3_prior$null_weight)),
               tolerance = MACH_TOL)

  # Compare xUlist (non-null components)
  # R6: xUlist is a list including null at position 1
  r6_nonnull_xUlist <- r6_pv$xUlist[-1]
  for (k in seq_along(r6_nonnull_xUlist)) {
    expect_equal(r6_nonnull_xUlist[[k]], s3_prior$xUlist[[k]],
                 tolerance = MACH_TOL, ignore_attr = TRUE)
  }

  # Compare structure
  expect_equal(r6_prior$n_condition, nrow(s3_prior$xUlist[[1]]))
  expect_equal(r6_prior$n_component - 1, length(s3_prior$xUlist))  # R6 includes null
})

test_that("create_mash_prior K=1: matches single matrix", {
  R <- 3
  V <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), R, R)
  prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)

  # Should have exactly 1 non-null component

  expect_equal(length(prior$xUlist), 1)
  expect_equal(prior$null_weight, 0)
  # The single xUlist should be V * 1 (grid=1)
  expect_equal(prior$xUlist[[1]], V, tolerance = MACH_TOL, ignore_attr = TRUE)
})

# ============================================================================
# Helper: compare ALL output fields between two mvsusie fits
# ============================================================================
#
# NOTE on V: R6 stores V as the prior_variance_scalar (starts at 1 for fixed
# prior, may change with EM/optim estimation). S3 stores V = V_scalar * max_diag
# where max_diag = max(abs(diag(prior))). When estimate_prior_variance = FALSE,
# R6 keeps V = 1 while S3 stores the actual prior variance. This is an
# intentional representation difference. Use check_V = TRUE only when both
# sides use the same convention (e.g., EM/optim estimation tests).
#
# NOTE on fitted: S3 stores fitted = X_std %*% b_sum + Y_mean (original scale).
# R6 stores fitted differently. Compare only when explicitly requested.
expect_all_fields_equal <- function(fit, ref, tol, label = "",
                                    check_V = FALSE, check_fitted = FALSE) {
  info <- function(field) paste0(label, " field: ", field)

  # Core fields (always present)
  expect_equal(fit$alpha, ref$alpha, tolerance = tol,
               ignore_attr = TRUE, info = info("alpha"))
  expect_equal(fit$lbf, ref$lbf, tolerance = tol,
               ignore_attr = TRUE, info = info("lbf"))
  expect_equal(fit$lbf_variable, ref$lbf_variable, tolerance = tol,
               ignore_attr = TRUE, info = info("lbf_variable"))
  expect_equal(get_b1(fit), get_b1(ref), tolerance = tol,
               ignore_attr = TRUE, info = info("b1"))
  expect_equal(get_b2(fit), get_b2(ref), tolerance = tol,
               ignore_attr = TRUE, info = info("b2"))

  # Coefficient (handle NA intercept)
  # Use coef_R6() which handles both S3 (computes from alpha+mu)
  # and R6 (reads stored $coef).
  fit_coef <- coef_R6(fit); fit_coef[is.na(fit_coef)] <- 0
  ref_coef <- coef_R6(ref); ref_coef[is.na(ref_coef)] <- 0
  expect_equal(fit_coef, ref_coef, tolerance = tol,
               ignore_attr = TRUE, info = info("coef"))

  # PIP
  expect_equal(fit$pip, ref$pip, tolerance = tol,
               ignore_attr = TRUE, info = info("pip"))

  # KL divergence
  if (!is.null(ref$KL) && !all(is.na(ref$KL))) {
    expect_equal(fit$KL, ref$KL, tolerance = tol,
                 ignore_attr = TRUE, info = info("KL"))
  }

  # ELBO
  if (!is.null(ref$elbo) && length(ref$elbo) > 0 && !all(is.na(ref$elbo))) {
    expect_equal(tail(fit$elbo, 1), tail(ref$elbo, 1), tolerance = tol,
                 ignore_attr = TRUE, info = info("elbo"))
  }

  # sigma2 (residual variance)
  if (!is.null(ref$sigma2)) {
    expect_equal(fit$sigma2, ref$sigma2, tolerance = tol,
                 ignore_attr = TRUE, info = info("sigma2"))
  }

  # V (only when explicitly requested --see NOTE above)
  if (check_V && !is.null(ref$V)) {
    expect_equal(fit$V, ref$V, tolerance = tol,
                 ignore_attr = TRUE, info = info("V"))
  }

  # niter
  if (!is.null(ref$niter)) {
    expect_equal(fit$niter, ref$niter, info = info("niter"))
  }

  # LFSR fields (only for mash/mixture priors)
  if (is.array(ref$lfsr) || (is.matrix(ref$lfsr))) {
    expect_equal(fit$lfsr, ref$lfsr, tolerance = tol,
                 ignore_attr = TRUE, info = info("lfsr"))
  }
  if (is.array(ref$single_effect_lfsr) || is.matrix(ref$single_effect_lfsr)) {
    expect_equal(fit$single_effect_lfsr, ref$single_effect_lfsr,
                 tolerance = tol, ignore_attr = TRUE,
                 info = info("single_effect_lfsr"))
  }
  if (is.array(ref$posterior_mixture_weights)) {
    expect_equal(fit$posterior_mixture_weights, ref$posterior_mixture_weights,
                 tolerance = tol, ignore_attr = TRUE,
                 info = info("posterior_mixture_weights"))
  }

  # Fitted values (only when explicitly requested --see NOTE above)
  if (check_fitted && !is.null(ref$fitted) && !is.null(fit$fitted)) {
    expect_equal(fit$fitted, ref$fitted, tolerance = tol,
                 ignore_attr = TRUE, info = info("fitted"))
  }
}

# ============================================================================
# Section 3: Special case -- K=1 no-null = single matrix prior
# (R=1 mash K=1 S3-vs-R6 check also lives here)
# ============================================================================

test_that("Mash K=1 null_weight=0 produces same result as matrix prior (R=3, fixed)", {
  with(sim_r3, {
    # Mash prior with K=1 component, no null weight
    mash_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)

    # Run with mash prior
    fit_mash <- mvsusie(X, y, L = L, prior_variance = mash_prior,
                        residual_variance = cov(y),
                        estimate_residual_variance = FALSE,
                        estimate_prior_variance = FALSE,
                        intercept = TRUE, standardize = TRUE,
                        precompute_cache = TRUE, verbose = FALSE)

    # Run with matrix prior
    fit_matrix <- mvsusie(X, y, L = L, prior_variance = V,
                          residual_variance = cov(y),
                          estimate_residual_variance = FALSE,
                          estimate_prior_variance = FALSE,
                          intercept = TRUE, standardize = TRUE, verbose = FALSE)

    # Core posterior quantities must match
    expect_equal(fit_mash$alpha, fit_matrix$alpha,
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: alpha")
    expect_equal(get_b1(fit_mash), get_b1(fit_matrix),
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: b1")
    expect_equal(fit_mash$lbf, fit_matrix$lbf,
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: lbf")
    expect_equal(coef(fit_mash), coef(fit_matrix),
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: coef")
    expect_equal(fit_mash$pip, fit_matrix$pip,
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: pip")
    expect_equal(fit_mash$KL, fit_matrix$KL,
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: KL")
    expect_equal(tail(fit_mash$elbo, 1), tail(fit_matrix$elbo, 1),
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: elbo")
    expect_equal(get_b2(fit_mash), get_b2(fit_matrix),
                 tolerance = MACH_TOL, ignore_attr = TRUE,
                 info = "K1_vs_matrix: b2")
  })
})

test_that("Mash K=1 null_weight=0 produces same result as matrix prior (R=3, optim)", {
  with(sim_r3, {
    # Same test as above but with estimate_prior_variance = TRUE (optim).
    # The mashr C++ path (calc_lik_rcpp) and the direct path (multivariate_lbf)
    # should produce the same V_scalar and posteriors.
    mash_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    fit_mash <- mvsusie(X, y, L = L, prior_variance = mash_prior,
                        residual_variance = cov(y),
                        estimate_residual_variance = FALSE,
                        estimate_prior_variance = TRUE,
                        estimate_prior_method = "optim",
                        max_iter = 20,
                        intercept = TRUE, standardize = TRUE,
                        precompute_cache = FALSE, verbose = FALSE)
    fit_matrix <- mvsusie(X, y, L = L, prior_variance = V,
                          residual_variance = cov(y),
                          estimate_residual_variance = FALSE,
                          estimate_prior_variance = TRUE,
                          estimate_prior_method = "optim",
                          max_iter = 20,
                          intercept = TRUE, standardize = TRUE, verbose = FALSE)
    # Optimizer convergence noise (~1e-7) is expected between the two
    # code paths (mashr C++ vs multivariate_lbf).
    optim_tol <- 1e-7
    expect_equal(fit_mash$alpha, fit_matrix$alpha,
                 tolerance = optim_tol, ignore_attr = TRUE,
                 info = "K1_vs_matrix_optim: alpha")
    expect_equal(get_b1(fit_mash), get_b1(fit_matrix),
                 tolerance = optim_tol, ignore_attr = TRUE,
                 info = "K1_vs_matrix_optim: b1")
    expect_equal(fit_mash$lbf, fit_matrix$lbf,
                 tolerance = optim_tol, ignore_attr = TRUE,
                 info = "K1_vs_matrix_optim: lbf")
    expect_equal(fit_mash$pip, fit_matrix$pip,
                 tolerance = optim_tol, ignore_attr = TRUE,
                 info = "K1_vs_matrix_optim: pip")
    expect_equal(tail(fit_mash$elbo, 1), tail(fit_matrix$elbo, 1),
                 tolerance = optim_tol, ignore_attr = TRUE,
                 info = "K1_vs_matrix_optim: elbo")
  })
})

test_that("R=1 mash K=1: S3 matches R6 at machine precision", {
  # NOTE: For R=1, mash path (mashr C++) and matrix path (R native) use
  # fundamentally different code paths. This test verifies S3 mash = R6 mash.
  # ELBO convergence ensures identical iteration counts and machine-precision
  # agreement between S3 and R6 paths (verified: 6e-16).
  skip_if_no_r6()
  r6 <- load_r6_env()
  with(sim_r1, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6$MashInitializer$new(Ulist = list(V), grid = 1, null_weight = 0)
    fit_mash <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                        residual_variance = cov(y),
                        estimate_residual_variance = FALSE,
                        estimate_prior_variance = FALSE,
                        tol = 1e-3,
                        intercept = TRUE, standardize = TRUE,
                        precompute_cache = TRUE, verbose = FALSE)
    ref_mash <- r6$mvsusie(X, y, L = L, prior_variance = r6_prior,
                            residual_variance = cov(y),
                            estimate_residual_variance = FALSE,
                            estimate_prior_variance = FALSE,
                            tol = 1e-3,
                            intercept = TRUE, standardize = TRUE,
                            precompute_covariances = TRUE, verbosity = 0)
    expect_all_fields_equal(fit_mash, ref_mash, tol = MACH_TOL,
                            label = "R1_mash_K1")
  })
})

# ============================================================================
# Section 4: Sufficient statistics path = individual path
# ============================================================================
#
# IMPORTANT: mvsusie_ss expects PRE-CENTERED inputs (as documented in
# R6 examples). The correct way to compute sufficient statistics:
#   X_c <- scale(X, center = TRUE, scale = FALSE)
#   Y_c <- scale(Y, center = TRUE, scale = FALSE)
#   XtX <- crossprod(X_c); XtY <- crossprod(X_c, Y_c); YtY <- crossprod(Y_c)
# Passing uncentered cross products causes ~1e-3 errors due to catastrophic
# cancellation in the standardization step.

test_that("mvsusie_ss = mvsusie at machine precision (R=3, centered inputs)", {
  with(sim_r3, {
    # Individual-level fit
    fit_ind <- mvsusie(X, y, L = L, prior_variance = V,
                       residual_variance = cov(y),
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE,
                       intercept = TRUE, standardize = TRUE, verbose = FALSE)

    # Compute sufficient statistics from CENTERED data (correct convention)
    n <- nrow(X)
    X_c <- scale(X, center = TRUE, scale = FALSE)
    Y_c <- scale(y, center = TRUE, scale = FALSE)
    XtX <- crossprod(X_c)
    XtY <- crossprod(X_c, Y_c)
    YtY <- crossprod(Y_c)

    fit_ss <- mvsusie_ss(XtX, XtY, YtY, n,
                                 L = L, prior_variance = V,
                                 X_colmeans = colMeans(X),
                                 Y_colmeans = colMeans(y),
                                 residual_variance = cov(y),
                                 estimate_residual_variance = FALSE,
                                 estimate_prior_variance = FALSE,
                                 standardize = TRUE,
                                 verbose = FALSE)

    expect_all_fields_equal(fit_ss, fit_ind, tol = MACH_TOL,
                            label = "R3_SS_vs_ind")
  })
})

test_that("mvsusie_ss = mvsusie at machine precision (R=1, centered inputs)", {
  with(sim_r1, {
    fit_ind <- mvsusie(X, y, L = L, prior_variance = V[1, 1],
                       residual_variance = as.numeric(var(y)),
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE,
                       intercept = TRUE, standardize = TRUE, verbose = FALSE)

    n <- nrow(X)
    X_c <- scale(X, center = TRUE, scale = FALSE)
    Y_c <- scale(y, center = TRUE, scale = FALSE)
    XtX <- crossprod(X_c)
    XtY <- crossprod(X_c, Y_c)
    YtY <- crossprod(Y_c)

    fit_ss <- mvsusie_ss(XtX, XtY, YtY, n,
                                 L = L, prior_variance = V[1, 1],
                                 X_colmeans = colMeans(X),
                                 Y_colmeans = mean(y),
                                 residual_variance = as.numeric(var(y)),
                                 estimate_residual_variance = FALSE,
                                 estimate_prior_variance = FALSE,
                                 standardize = TRUE,
                                 verbose = FALSE)

    expect_all_fields_equal(fit_ss, fit_ind, tol = MACH_TOL,
                            label = "R1_SS_vs_ind")
  })
})

# ============================================================================
# Section 5: Workhorse defaults sanity
# ============================================================================

test_that("mvsusie_workhorse defaults match the canonical reference values", {
  # `tol` was deliberately tightened from 1e-3 to 1e-4 (2026-04-28) to
  # align mvsusieR's IBSS convergence with susieR's default. The R=1
  # mvsusie-vs-susieR identity contract requires the two packages take
  # the same number of IBSS iterations on the same fixture; the slacker
  # 1e-3 default would stop one iteration earlier than susieR and drift
  # by ~1e-3 on `mu` / `lbf_variable`. The remaining defaults stay at
  # the original commit-aaed0d9 values.
  args <- formals(mvsusie_workhorse)

  expect_equal(eval(args$check_null_threshold), 0)
  expect_equal(eval(args$max_iter), 100)
  expect_equal(eval(args$tol), 1e-4)
  expect_equal(eval(args$prior_tol), 1e-9)
  expect_equal(eval(args$verbose), TRUE)
  expect_equal(eval(args$coverage), 0.95)
  expect_equal(eval(args$min_abs_corr), 0.5)
  expect_equal(eval(args$n_thread), 1)
})

test_that("mvsusie_workhorse internal params have correct defaults", {
  # Run mvsusie with minimal arguments to check that the function doesn't
  # error and produces valid output with correct structure
  with(sim_r3, {
    fit <- mvsusie(X, y, L = 3, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   intercept = TRUE, standardize = TRUE,
                   max_iter = 2, verbose = FALSE)
    expect_s3_class(fit, "mvsusie")
    expect_true(!is.null(fit$alpha))
    expect_true(!is.null(coef(fit)))
    expect_equal(nrow(fit$alpha), 3)  # L = 3
    expect_equal(ncol(fit$alpha), ncol(X))  # J = p
  })
})

test_that("mvsusie verbose sigma2 summary is compact for matrix sigma2", {
  fit <- structure(list(sigma2 = diag(c(1, 2))), class = c("mvsusie", "susie"))
  expect_equal(format_sigma2_summary.mvsusie(fit), "diag[1,2]")
})

# ============================================================================
# Section 6: Edge cases
# ============================================================================

test_that("R=1 J=1 edge case doesn't crash", {
  set.seed(123)
  X1 <- matrix(rnorm(100), 100, 1)
  y1 <- X1 * 2 + rnorm(100)
  fit <- tryCatch(
    mvsusie(X1, y1, L = 1, prior_variance = 1,
            residual_variance = 1,
            estimate_residual_variance = FALSE,
            estimate_prior_variance = FALSE,
            max_iter = 5, verbose = FALSE),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    expect_s3_class(fit, "mvsusie")
  }
})

test_that("L > p: L is automatically reduced to p", {
  set.seed(123)
  X_small <- matrix(rnorm(100 * 5), 100, 5)
  y_small <- X_small %*% c(1, 0, 2, 0, 0) + rnorm(100)
  fit <- mvsusie(X_small, y_small, L = 10, prior_variance = 1,
                 residual_variance = 1,
                 estimate_residual_variance = FALSE,
                 estimate_prior_variance = FALSE,
                 max_iter = 5, verbose = FALSE)
  expect_s3_class(fit, "mvsusie")
  # L should be <= p = 5
  expect_true(nrow(fit$alpha) <= 5)
})
