# ============================================================================
# Reference tests: S3 path vs R6 (master branch)
#
# These tests compare the S3 refactored code against the original R6
# implementation by loading the master branch source via pkgload.
# The R6 package is loaded first to capture reference functions,
# then the S3 (development) package is loaded for testing.
#
# NOTE on convergence:
#   S3 and R6 have identical per-iteration math (verified at machine
#   precision in test_math_components.R).  However, the R6 IBSS loop
#   includes a check_null_threshold warmup (skipped for iterations
#   1-10, active after) that can snap small V scalars to 0 and
#   accelerate convergence.  The S3 path (susieR workhorse) does not
#   implement this warmup, so niter can differ.  This is a known
#   procedural difference, not a math bug, and niter is NOT compared
#   in these tests.
#
# NOTE on fitted values:
#   S3's fitted = Xr + Y_mean (susieR convention); R6's fitted = Xr
#   (centered scale).  predict() matches at machine precision.  The
#   tests compare fitted after removing the constant Y_mean offset.
#
# NOTE on V output:
#   The V output field is stored differently: S3 stores the V_scalar
#   multiplier while R6 stores the effective prior variance.  V is
#   only compared for R=1 scalar prior where the conventions coincide.
# ============================================================================

# --------------------------------------------------------------------------
# Setup: load R6 reference from master, capture functions, then load S3
# --------------------------------------------------------------------------

repo_dir <- tryCatch({
  gcd <- system2("git", c("rev-parse", "--git-common-dir"),
                 stdout = TRUE, stderr = FALSE)
  normalizePath(file.path(gcd, ".."))
}, error = function(e) NULL)

# Cached R6 reference functions
.r6_funcs <- NULL

load_r6_reference <- function() {
  if (!is.null(.r6_funcs)) return(.r6_funcs)
  if (is.null(repo_dir)) stop("Cannot determine git repository root")

  ref_source <- file.path(tempdir(), "mvsusieR_r6_ref")
  if (!dir.exists(ref_source)) {
    system2("git", c("-C", repo_dir, "worktree", "remove", "--force", ref_source),
            stdout = FALSE, stderr = FALSE)
    ret <- system2("git", c("-C", repo_dir, "worktree", "add", "--detach",
                            ref_source, "master"),
                   stdout = FALSE, stderr = FALSE)
    if (ret != 0) stop("Failed to create worktree from master")
  }

  # Load R6 package into a separate environment
  # We use a fresh environment to avoid namespace conflicts
  ref <- pkgload::load_all(ref_source, export_all = TRUE, quiet = TRUE)

  # Capture R6 functions into a standalone environment that won't
  # be affected when we reload the S3 package
  r6_env <- new.env(parent = ref$env)
  funcs <- list(
    mvsusie = ref$env$mvsusie,
    MashInitializer = ref$env$MashInitializer,
    create_mixture_prior = ref$env$create_mixture_prior,
    mvsusie_sim1 = ref$env$mvsusie_sim1,
    env = ref$env
  )

  # Reload the S3 (development) package to restore S3 functions
  dev_source <- tryCatch(
    rprojroot::find_root(rprojroot::is_r_package),
    error = function(e) normalizePath(file.path(getwd(), "../.."))
  )
  pkgload::load_all(dev_source, export_all = TRUE, quiet = TRUE)

  .r6_funcs <<- funcs
  funcs
}

# Load R6 eagerly --these tests require the master branch R6 reference.
# Fail hard if R6 can't be loaded (master should always be available).
ensure_r6_loaded <- function() {
  load_r6_reference()
}

# R6 wrapper functions
# Translates S3 parameter names to R6 equivalents:
#   precompute_cache -> precompute_covariances (R6 uses old name)
#   verbose -> verbosity (R6 uses integer verbosity, S3 uses boolean verbose)
r6_mvsusie <- function(...) {
  r6 <- load_r6_reference()
  args <- list(...)
  # Translate S3 parameter names to R6 equivalents
  if ("precompute_cache" %in% names(args)) {
    args$precompute_covariances <- args$precompute_cache
    args$precompute_cache <- NULL
  }
  # R6 master uses verbosity (integer); S3 uses verbose (boolean)
  if ("verbose" %in% names(args)) {
    args$verbosity <- if (isTRUE(args$verbose)) 2 else 0
    args$verbose <- NULL
  }
  # R6 master doesn't have estimate_prior_mixture_weights; strip it
  args$estimate_prior_mixture_weights <- NULL
  do.call(r6$mvsusie, args)
}

r6_MashInitializer <- function(...) {
  r6 <- load_r6_reference()
  r6$MashInitializer$new(...)
}

r6_create_mixture_prior <- function(...) {
  r6 <- load_r6_reference()
  r6$create_mixture_prior(...)
}

r6_fitted_g_prior <- function(fitted_g, null_weight = NULL, weights_tol = 1e-10,
                              mixture_length = 40, include_indices = NULL) {
  # Replicate R6's create_mash_prior(fitted_g=...) logic:
  # Extract weights, then call MashInitializer$new(Ulist, grid, ...)
  if (fitted_g$usepointmass) {
    prior_weights <- fitted_g$pi[-1]
    if (is.null(null_weight)) null_weight <- fitted_g$pi[1]
  } else {
    prior_weights <- fitted_g$pi
    if (is.null(null_weight)) null_weight <- 0
  }
  r6_MashInitializer(Ulist = fitted_g$Ulist, grid = fitted_g$grid,
                     prior_weights = prior_weights, null_weight = null_weight,
                     weights_tol = weights_tol, top_mixtures = mixture_length,
                     include_conditions = include_indices)
}

# Generate simulation data deterministically
set.seed(1)
sim1 <- mvsusie_sim1(n = 100, p = 100, r = 1, s = 3)
sim1$L <- 10
set.seed(1)
sim3 <- mvsusie_sim1(n = 100, p = 100, r = 3, s = 3)
sim3$L <- 10
set.seed(1)
sim_miss <- mvsusie_sim1(n = 100, p = 100, r = 1, s = 3, y_missing = 0.2)
sim_miss$L <- 10

# Tolerance levels
tol_tight   <- 1e-8   # Fixed prior, same code path
tol_em      <- 5e-2   # EM/optim: procedural convergence differs
tol_em_miss <- 2.5e-1 # EM with missing data
tol_mash    <- 5e-3   # Mash: R6 MashRegression vs S3 multivariate path

# Helper: compare core model outputs
#
# Known convention differences (NOT compared):
#   - niter: different procedural convergence (check_null warmup)
#   - fitted: S3 includes Y_mean offset, R6 doesn't. We compare
#     fitted after removing the constant offset (mean difference).
#   - V: different parameterisation for R>1 matrix priors
#   - coef intercept: follows from fitted convention
expect_ref_equal <- function(fit, ref, tol = tol_tight,
                             check_elbo = FALSE, check_V = FALSE,
                             check_lfsr = TRUE) {
  # Core fields
  expect_equal(fit$alpha, ref$alpha, tolerance = tol, check.attributes = FALSE)
  expect_equal(fit$lbf,   ref$lbf,   tolerance = tol, check.attributes = FALSE)
  expect_equal(fit$b1,    ref$b1,    tolerance = tol, check.attributes = FALSE)

  # Coef: compare slopes only (intercept convention differs).
  # When intercept=FALSE, S3 returns 0 and R6 returns NA --treat NAs as 0.
  fit_coef <- fit$coef; fit_coef[is.na(fit_coef)] <- 0
  ref_coef <- ref$coef; ref_coef[is.na(ref_coef)] <- 0
  if (is.matrix(fit_coef) && nrow(fit_coef) > 1) {
    # Skip intercept row (first row)
    expect_equal(fit_coef[-1, , drop = FALSE], ref_coef[-1, , drop = FALSE],
                 tolerance = tol, check.attributes = FALSE)
  } else {
    expect_equal(fit_coef, ref_coef, tolerance = tol, check.attributes = FALSE)
  }

  # b2 (alpha-weighted second moment diagonal)
  if (!is.null(ref$b2))
    expect_equal(fit$b2, ref$b2, tolerance = tol, check.attributes = FALSE)

  # PIP (posterior inclusion probability)
  if (!is.null(ref$pip))
    expect_equal(fit$pip, ref$pip, tolerance = tol, check.attributes = FALSE)

  # KL divergence
  if (!is.null(ref$KL) && !all(is.na(ref$KL)))
    expect_equal(fit$KL, ref$KL, tolerance = tol, check.attributes = FALSE)

  # sigma2 (residual variance)
  if (!is.null(ref$sigma2))
    expect_equal(fit$sigma2, ref$sigma2, tolerance = tol, check.attributes = FALSE)

  # niter: NOT compared. S3 and R6 have identical per-iteration math
  # (verified in test_math_components.R) but differ in procedural
  # convergence (check_null warmup in R6 IBSS loop).

  # Fitted: compare on centered scale (remove constant Y_mean offset).
  # S3 fitted = Xr + Y_mean (susieR convention); R6 fitted = Xr.
  # predict() matches at machine precision (verified in test_math_components).
  if (!is.null(ref$fitted) && !is.null(fit$fitted)) {
    f <- as.matrix(fit$fitted)
    r <- as.matrix(ref$fitted)
    fit_centered <- sweep(f, 2, colMeans(f))
    ref_centered <- sweep(r, 2, colMeans(r))
    expect_equal(fit_centered, ref_centered, tolerance = tol,
                 check.attributes = FALSE)
  }

  # LFSR fields (mash/mixture priors only).
  # Skipped when check_lfsr = FALSE --downstream of procedural trimming
  # differences that cause R6 to zero out borderline effects.
  if (check_lfsr) {
    if (is.array(ref$lfsr) || is.matrix(ref$lfsr))
      expect_equal(fit$lfsr, ref$lfsr, tolerance = tol, check.attributes = FALSE)
    if (is.array(ref$single_effect_lfsr) || is.matrix(ref$single_effect_lfsr))
      expect_equal(fit$single_effect_lfsr, ref$single_effect_lfsr,
                   tolerance = tol, check.attributes = FALSE)
    if (is.array(ref$mixture_weights))
      expect_equal(fit$mixture_weights, ref$mixture_weights,
                   tolerance = tol, check.attributes = FALSE)
  }

  # V (prior variance) and ELBO.
  # S3 V output convention: for K=1 matrix prior with estimate_prior_variance,
  # out$V = V_scalar * V_structure[[1]][1,1] (absolute scale).
  # R6 V output: V_scalar (the multiplier).
  # check_V = TRUE: direct comparison (works for R=1 scalar priors).
  # check_V = "effective": compare effective prior V*V_structure matrices;
  #   reconstructs S3 V_scalar from output, compares against R6 V_scalar.
  if (identical(check_V, TRUE) && !is.null(ref$V))
    expect_equal(fit$V, ref$V, tolerance = tol, check.attributes = FALSE)
  if (identical(check_V, "effective") && !is.null(ref$V) &&
      !is.null(fit$V_structure)) {
    V_struct <- fit$V_structure
    if (length(V_struct) == 1 && is.matrix(V_struct[[1]])) {
      # K=1 matrix prior: S3 V = V_scalar * V_structure[1,1]
      s3_V_scalar <- fit$V / V_struct[[1]][1, 1]
      r6_V_scalar <- ref$V
      expect_equal(s3_V_scalar, r6_V_scalar, tolerance = tol,
                   check.attributes = FALSE,
                   info = "V_scalar (from effective V comparison)")
      # Also compare effective prior matrices for each effect
      for (l in seq_along(s3_V_scalar)) {
        s3_eff <- s3_V_scalar[l] * V_struct[[1]]
        r6_eff <- r6_V_scalar[l] * V_struct[[1]]
        expect_equal(s3_eff, r6_eff, tolerance = tol,
                     check.attributes = FALSE,
                     info = paste0("Effective prior V*V_structure, effect ", l))
      }
    }
  }
  if (check_elbo && !is.null(ref$elbo))
    expect_equal(tail(fit$elbo, 1), tail(ref$elbo, 1),
                 tolerance = tol, check.attributes = FALSE)
}

# ============================================================================
# R=1 missing data
# ============================================================================

test_that("R=1, missing data, fixed variance matches R6", {
  ensure_r6_loaded()
  with(sim_miss, {
    fit <- mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   intercept = TRUE, standardize = TRUE, verbose = FALSE)
    ref <- r6_mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      intercept = TRUE, standardize = TRUE, verbose = FALSE)
    # Core math fields match at tight tolerance; fitted differs for
    # missing observations (convention: R6 zeros missing rows).
    # Compare math fields individually, skip fitted for missing data.
    expect_equal(fit$alpha, ref$alpha, tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$lbf,   ref$lbf,   tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$b1,    ref$b1,    tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$pip,   ref$pip,   tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$sigma2, ref$sigma2, tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(tail(fit$elbo, 1), tail(ref$elbo, 1),
                 tolerance = tol_tight, check.attributes = FALSE)
    # Fitted: compare only observed rows
    obs <- !is.na(y_missing)
    expect_equal(fit$fitted[obs] - mean(fit$fitted[obs]),
                 ref$fitted[obs] - mean(ref$fitted[obs]),
                 tolerance = tol_tight, check.attributes = FALSE)
  })
})

test_that("R=1, missing data, EM (10 iter) matches R6 at tight tol", {
  ensure_r6_loaded()
  with(sim_miss, {
    fit <- mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   max_iter = 10,
                   intercept = TRUE, standardize = TRUE, verbose = FALSE)
    ref <- r6_mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      max_iter = 10,
                      intercept = TRUE, standardize = TRUE, verbose = FALSE)
    # Compare math fields; fitted only for observed rows (convention differs
    # for missing observations: R6 zeros them, S3 doesn't).
    expect_equal(fit$alpha, ref$alpha, tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$lbf,   ref$lbf,   tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$b1,    ref$b1,    tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(fit$pip,   ref$pip,   tolerance = tol_tight,
                 check.attributes = FALSE)
    expect_equal(tail(fit$elbo, 1), tail(ref$elbo, 1),
                 tolerance = tol_tight, check.attributes = FALSE)
    obs <- !is.na(y_missing)
    expect_equal(fit$fitted[obs] - mean(fit$fitted[obs]),
                 ref$fitted[obs] - mean(ref$fitted[obs]),
                 tolerance = tol_tight, check.attributes = FALSE)
  })
})

# ============================================================================
# fitted_g pathway (create_mixture_prior with fitted_g from mashr::mash)
#
# The R6 had create_mash_prior(fitted_g=...) which accepted fitted_g from
# mashr::mash() output. The S3 merges this into create_mixture_prior(fitted_g=...).
# These tests verify faithful translation of the fitted_g pathway by:
#   1. Creating a fitted_g via mashr::mash()
#   2. Building priors from fitted_g via both R6 and S3 pathways
#   3. Running mvsusie with both and comparing outputs
# ============================================================================

# Helper: create a mock fitted_g directly (no mashr::mash() needed)
# This mimics what mashr::mash() returns in fitted_g.
create_mock_fitted_g <- function(R, K = 3, n_grid = 3) {
  # Create K base covariance matrices
  Ulist <- list()
  for (k in seq_len(K)) {
    # Simple: identity, singleton_1, rank-1 shared
    if (k == 1) {
      Ulist[[k]] <- diag(R)
    } else if (k == 2) {
      mat <- matrix(0, R, R)
      mat[1, 1] <- 1
      Ulist[[k]] <- mat
    } else {
      Ulist[[k]] <- matrix(1, R, R)
    }
  }
  grid <- seq(0.5, 2, length.out = n_grid)
  n_components <- K * n_grid
  # Simulate mash-like weights: non-uniform with point mass
  set.seed(42)
  raw_wts <- runif(n_components)
  null_wt <- 0.15
  pi <- c(null_wt, raw_wts / sum(raw_wts) * (1 - null_wt))
  list(
    Ulist = Ulist,
    grid = grid,
    pi = pi,
    usepointmass = TRUE
  )
}

test_that("R=3, fitted_g prior, fixed variance matches R6", {
  ensure_r6_loaded()
  with(sim3, {
    fg <- create_mock_fitted_g(R = 3)
    s3_prior <- create_mixture_prior(fitted_g = fg)
    r6_prior <- r6_fitted_g_prior(fg)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   estimate_prior_mixture_weights = FALSE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_cache = TRUE, verbose = FALSE)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = FALSE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_cache = TRUE, verbose = FALSE)
    expect_ref_equal(fit, ref, tol = tol_tight, check_elbo = TRUE)
  })
})

test_that("R=3, fitted_g prior, EM (10 iter) matches R6", {
  ensure_r6_loaded()
  with(sim3, {
    fg <- create_mock_fitted_g(R = 3)
    s3_prior <- create_mixture_prior(fitted_g = fg)
    r6_prior <- r6_fitted_g_prior(fg)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   estimate_prior_mixture_weights = FALSE,
                   max_iter = 10,
                   intercept = TRUE, standardize = TRUE,
                   precompute_cache = FALSE, verbose = FALSE)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      estimate_prior_mixture_weights = FALSE,
                      max_iter = 10,
                      intercept = TRUE, standardize = TRUE,
                      precompute_cache = FALSE, verbose = FALSE)
    expect_ref_equal(fit, ref, tol = tol_tight, check_elbo = TRUE)
  })
})

test_that("R=3, fitted_g prior, null_weight override matches R6", {
  ensure_r6_loaded()
  with(sim3, {
    fg <- create_mock_fitted_g(R = 3)
    # Override null_weight to 0 (ignore mash-estimated null weight)
    s3_prior <- create_mixture_prior(fitted_g = fg, null_weight = 0)
    r6_prior <- r6_fitted_g_prior(fg, null_weight = 0)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   estimate_prior_mixture_weights = FALSE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_cache = TRUE, verbose = FALSE)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      estimate_prior_mixture_weights = FALSE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_cache = TRUE, verbose = FALSE)
    expect_ref_equal(fit, ref, tol = tol_tight, check_elbo = TRUE)
  })
})

test_that("fitted_g prior structure matches R6 MashInitializer", {
  # Test that create_mixture_prior(fitted_g=...) produces a prior with
  # identical structure to the R6's MashInitializer$new(Ulist, grid, ...)
  #
  # R6 convention: prior_variance$xUlist includes null at position 1,
  #   prior_variance$pi includes null weight at position 1.
  # S3 convention: xUlist has only non-null components,
  #   null_weight stored separately.
  ensure_r6_loaded()
  fg <- create_mock_fitted_g(R = 3)
  s3_prior <- create_mixture_prior(fitted_g = fg)
  r6_prior <- r6_fitted_g_prior(fg)
  r6_pv <- r6_prior$prior_variance

  # R6 xUlist[1] is null, R6 xUlist[2:end] are non-null
  r6_nonnull <- r6_pv$xUlist[-1]
  expect_equal(s3_prior$n_component, length(r6_nonnull))

  # Compare null weights: R6 pi[1] = null weight
  expect_equal(s3_prior$null_weight, unname(r6_pv$pi[1]), tolerance = 1e-10)

  # Compare non-null mixture weights (normalized)
  r6_pi_nonnull <- r6_pv$pi[-1]
  r6_pi_nonnull <- r6_pi_nonnull / sum(r6_pi_nonnull)
  expect_equal(as.vector(s3_prior$pi), as.vector(r6_pi_nonnull),
               tolerance = 1e-10, check.attributes = FALSE)

  # Compare actual prior matrices (non-null only)
  for (k in seq_along(s3_prior$xUlist)) {
    expect_equal(s3_prior$xUlist[[k]], r6_nonnull[[k]],
                 tolerance = 1e-14, check.attributes = FALSE)
  }
})

test_that("fitted_g with usepointmass=FALSE works", {
  # Edge case: mash output without point mass
  fg <- create_mock_fitted_g(R = 3)
  fg$usepointmass <- FALSE
  fg$pi <- fg$pi[-1]  # Remove null weight from pi
  fg$pi <- fg$pi / sum(fg$pi)
  s3_prior <- create_mixture_prior(fitted_g = fg)
  expect_equal(s3_prior$null_weight, 0)
  expect_equal(sum(s3_prior$pi), 1, tolerance = 1e-10)
})

test_that("create_mixture_prior rejects invalid fitted_g", {
  # Missing required fields
  expect_error(create_mixture_prior(fitted_g = list(pi = 1)),
               "Cannot find Ulist")
  expect_error(create_mixture_prior(fitted_g = list(Ulist = list())),
               "Cannot find pi")
  # Cannot combine fitted_g with other pathway
  expect_error(create_mixture_prior(fitted_g = list(pi = 1, Ulist = list(),
                                                    grid = 1,
                                                    usepointmass = TRUE),
                                    R = 3),
               "exactly one")
})

# ============================================================================
# Non-common-cov C++ eigendecomposition path
#
# NOTE: The C++ non-common-cov eigen path (precompute_eigen_cache_non_common_rcpp,
# loglik_non_common_rcpp, posterior_non_common_rcpp) is verified against the
# pure R implementation at machine precision in test_missing_y_methods.R
# (30+ C++ vs R comparison tests across multiple R/K/J configurations).
#
# S3-vs-R6 comparison for non-common-cov is not done here because the
# S3 eigendecomposition path and the R6 mashr Cholesky path use different
# numerical algorithms (eigendecomposition vs mashr's dmvnorm/posterior_cov)
# that can diverge in LFSR and mixture weights for borderline components,
# even when core quantities (alpha, pip, elbo) agree.
# ============================================================================

# ============================================================================
# API backward compatibility audit tests
#
# These tests verify that the S3 refactored API maintains backward
# compatibility with the R6 master branch:
#   - Output object fields that users may access
#   - Parameter names and defaults
#   - Function exports
# ============================================================================

test_that("S3 output has expected core fields", {
  ensure_r6_loaded()
  R <- sim3$r
  s3_prior <- create_mixture_prior(R = R)
  dev <- mvsusie(sim3$X, sim3$y, L = sim3$L, prior_variance = s3_prior,
                 estimate_residual_variance = FALSE,
                 max_iter = 5, verbose = FALSE)

  # Core output fields that must always be present
  required <- c("alpha", "b1", "b2", "coef", "fitted", "intercept",
                "pip", "sets", "sigma2", "V", "elbo", "niter",
                "lbf", "lbf_variable", "KL", "convergence",
                "outcome_names", "variable_names")
  missing <- setdiff(required, names(dev))
  expect_equal(length(missing), 0,
               info = paste("Missing fields:", paste(missing, collapse = ", ")))
})

test_that("track_fit produces trace in output", {
  ensure_r6_loaded()
  R <- sim3$r
  s3_prior <- create_mixture_prior(R = R)

  # track_fit = FALSE (default): no trace
  fit_no <- mvsusie(sim3$X, sim3$y, L = sim3$L, prior_variance = s3_prior,
                    estimate_residual_variance = FALSE,
                    max_iter = 5, verbose = FALSE, track_fit = FALSE)
  expect_null(fit_no$trace)

  # track_fit = TRUE: trace present with per-iteration snapshots
  fit_yes <- mvsusie(sim3$X, sim3$y, L = sim3$L, prior_variance = s3_prior,
                     estimate_residual_variance = FALSE,
                     max_iter = 5, verbose = FALSE, track_fit = TRUE)
  expect_false(is.null(fit_yes$trace))
  expect_true(length(fit_yes$trace) > 0)
  # Each trace entry should have key model fields
  entry <- fit_yes$trace[[1]]
  expect_true(all(c("alpha", "V", "sigma2") %in% names(entry)))
})

test_that("S3 output has outcome_names (replaces R6 condition_names)", {
  ensure_r6_loaded()
  R <- sim3$r
  s3_prior <- create_mixture_prior(R = R)
  dev <- mvsusie(sim3$X, sim3$y, L = sim3$L, prior_variance = s3_prior,
                 estimate_residual_variance = FALSE,
                 max_iter = 5, verbose = FALSE)
  expect_false(is.null(dev$outcome_names))
  expect_equal(length(dev$outcome_names), R)
})

test_that("estimate_residual_variance defaults: TRUE for individual, FALSE for ss", {
  # Verify by checking the function formals directly
  mvsusie_formals <- formals(mvsusie)
  expect_true(mvsusie_formals$estimate_residual_variance)

  mvsusie_ss_formals <- formals(mvsusie_ss)
  expect_false(mvsusie_ss_formals$estimate_residual_variance)
})

test_that("create_mixture_prior backward compatibility: R and mixture_prior paths", {
  ensure_r6_loaded()
  # R path: R6 and S3 should produce compatible priors
  r6_prior <- r6_create_mixture_prior(R = 3)
  s3_prior <- create_mixture_prior(R = 3)

  # Both should have the same number of components
  expect_equal(s3_prior$n_component,
               length(r6_prior$prior_variance$xUlist) - 1)  # R6 includes null

  # S3 xUlist should match R6 xUlist (excluding null at position 1)
  r6_xUlist <- r6_prior$prior_variance$xUlist[-1]
  for (k in seq_along(r6_xUlist)) {
    expect_equal(s3_prior$xUlist[[k]], r6_xUlist[[k]],
                 tolerance = 1e-12,
                 info = paste("Component", k))
  }
})
