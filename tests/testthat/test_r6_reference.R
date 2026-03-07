# ============================================================================
# Reference tests: S3 path vs R6 (master branch)
#
# These tests compare the S3 refactored code against the original R6
# implementation by loading the master branch source via pkgload.
# The R6 package is loaded first to capture reference functions,
# then the S3 (development) package is loaded for testing.
#
# NOTE on V (prior variance) comparison:
#   R6 stores V as a scaling factor; S3 stores V as max(abs(diag(prior))).
#   Both produce the same effective prior. For R=1 scalar prior, both
#   store the absolute variance and can be compared directly.
#
# NOTE on EM tolerance:
#   EM/optim tests use looser tolerance because R6 and S3 use different
#   EM update formulas, converging to similar but not identical optima.
# ============================================================================

context("S3 vs R6 reference tests")

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

skip_if_no_r6 <- function() {
  tryCatch(
    load_r6_reference(),
    error = function(e)
      skip(paste("R6 reference not available:", e$message))
  )
}

# R6 wrapper functions
r6_mvsusie <- function(...) {
  r6 <- load_r6_reference()
  r6$mvsusie(...)
}

r6_MashInitializer <- function(...) {
  r6 <- load_r6_reference()
  r6$MashInitializer$new(...)
}

r6_create_mixture_prior <- function(...) {
  r6 <- load_r6_reference()
  r6$create_mixture_prior(...)
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
tol_em      <- 5e-2   # EM/optim: different update formulas
tol_em_miss <- 2.5e-1 # EM with missing data
tol_mash    <- 5e-3   # Mash: R6 MashRegression vs S3 multivariate path

# Helper: compare core model outputs
expect_ref_equal <- function(fit, ref, tol = tol_tight,
                             check_elbo = FALSE, check_V = FALSE) {
  expect_equal(fit$alpha, ref$alpha, tolerance = tol, check.attributes = FALSE)
  expect_equal(fit$lbf,   ref$lbf,   tolerance = tol, check.attributes = FALSE)
  expect_equal(fit$b1,    ref$b1,    tolerance = tol, check.attributes = FALSE)
  # When intercept=FALSE, S3 returns 0 and R6 returns NA for the intercept
  # coefficient. Treat NAs as 0 for comparison purposes.
  fit_coef <- fit$coef; fit_coef[is.na(fit_coef)] <- 0
  ref_coef <- ref$coef; ref_coef[is.na(ref_coef)] <- 0
  expect_equal(fit_coef, ref_coef, tolerance = tol, check.attributes = FALSE)
  if (check_V && !is.null(ref$V))
    expect_equal(fit$V, ref$V, tolerance = tol, check.attributes = FALSE)
  if (check_elbo && !is.null(ref$elbo))
    expect_equal(tail(fit$elbo, 1), tail(ref$elbo, 1),
                 tolerance = tol, check.attributes = FALSE)
}

# ============================================================================
# R=1 matrix prior (scalar prior_variance)
# ============================================================================

test_that("R=1, matrix prior, fixed variance matches R6", {
  skip_if_no_r6()
  with(sim1, {
    fit <- mvsusie(X, y, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight,
                     check_elbo = TRUE, check_V = TRUE)
  })
})

test_that("R=1, matrix prior, EM prior estimation matches R6", {
  skip_if_no_r6()
  with(sim1, {
    fit <- mvsusie(X, y, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em, check_V = TRUE)
  })
})

# ============================================================================
# R=3 matrix prior
# ============================================================================

test_that("R=3, matrix prior, fixed variance matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = L, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight,
                     check_elbo = TRUE, check_V = FALSE)
  })
})

test_that("R=3, matrix prior, EM prior estimation matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = L, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em, check_V = FALSE)
  })
})

# ============================================================================
# R=3 mash / mixture prior
# ============================================================================

test_that("R=3, mash prior K=1, fixed variance matches R6", {
  skip_if_no_r6()
  with(sim3, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6_MashInitializer(Ulist = list(V), grid = 1, null_weight = 0)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_covariances = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_covariances = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight,
                     check_elbo = TRUE, check_V = FALSE)
  })
})

test_that("R=3, mash prior K=1, EM prior estimation matches R6", {
  skip_if_no_r6()
  with(sim3, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6_MashInitializer(Ulist = list(V), grid = 1, null_weight = 0)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_covariances = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_covariances = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em, check_V = FALSE)
  })
})

test_that("R=3, create_mixture_prior K=1, EM matches R6", {
  skip_if_no_r6()
  with(sim3, {
    s3_prior <- create_mixture_prior(mixture_prior = list(matrices = list(V)))
    r6_prior <- r6_create_mixture_prior(mixture_prior = list(matrices = list(V)))
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_covariances = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_covariances = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em,
                     check_elbo = TRUE, check_V = FALSE)
  })
})

# ============================================================================
# R=1 mash prior
# ============================================================================

test_that("R=1, mash prior K=1, fixed variance matches R6", {
  skip_if_no_r6()
  with(sim1, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6_MashInitializer(Ulist = list(V), grid = 1, null_weight = 0)
    fit <- mvsusie(X, y, L = L, prior_variance = s3_prior,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE,
                   precompute_covariances = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = r6_prior,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE,
                      precompute_covariances = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_mash, check_V = FALSE)
  })
})

# ============================================================================
# R=1 missing data
# ============================================================================

test_that("R=1, missing data, fixed variance matches R6", {
  skip_if_no_r6()
  with(sim_miss, {
    fit <- mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight,
                     check_elbo = TRUE, check_V = TRUE)
  })
})

test_that("R=1, missing data, EM prior estimation matches R6", {
  skip_if_no_r6()
  with(sim_miss, {
    fit <- mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                   residual_variance = as.numeric(var(y)),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "EM",
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y_missing, L = L, prior_variance = V[1, 1],
                      residual_variance = as.numeric(var(y)),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "EM",
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em_miss, check_V = TRUE)
  })
})

# ============================================================================
# Residual variance estimation
# ============================================================================

test_that("R=3, matrix prior, estimate residual variance matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = L, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight, check_V = FALSE)
    expect_equal(fit$sigma2, ref$sigma2,
                 tolerance = tol_tight, check.attributes = FALSE)
  })
})

# ============================================================================
# Center/scale options
# ============================================================================

test_that("R=3, no centering, no standardization matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = L, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = TRUE,
                   intercept = FALSE, standardize = FALSE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = TRUE,
                      intercept = FALSE, standardize = FALSE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight, check_V = FALSE)
  })
})

# ============================================================================
# L=1 single effect
# ============================================================================

test_that("R=3, L=1, matrix prior matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = 1, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = FALSE,
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = 1, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = FALSE,
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_tight, check_V = FALSE)
  })
})

# ============================================================================
# Optim prior estimation
# ============================================================================

test_that("R=3, matrix prior, optim estimation matches R6", {
  skip_if_no_r6()
  with(sim3, {
    fit <- mvsusie(X, y, L = L, prior_variance = V,
                   residual_variance = cov(y),
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = "optim",
                   compute_objective = FALSE,
                   intercept = TRUE, standardize = TRUE, verbosity = 0)
    ref <- r6_mvsusie(X, y, L = L, prior_variance = V,
                      residual_variance = cov(y),
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = "optim",
                      compute_objective = FALSE,
                      intercept = TRUE, standardize = TRUE, verbosity = 0)
    expect_ref_equal(fit, ref, tol = tol_em, check_V = FALSE)
  })
})
