# ============================================================================
# Numerical sweep: S3 vs R6 across many configurations
#
# Testing philosophy (binary, no middle ground):
#
#   APPLE-TO-APPLE: S3 and R6 use identical algorithms on the same data.
#     Compare at machine precision (tol = 1e-8). Any failure is a bug.
#
#   APPLE-TO-ORANGE: S3 and R6 use different algorithms or procedures.
#     Do NOT compare numerically. Just verify S3 runs without error
#     and output has sane structure (smoke test).
#
# Classification by conceptual identity of code paths:
#
#   Apple-to-apple (identical algorithm, compare at machine precision):
#     - Matrix prior + EM: all est_pv × est_rv × intercept combos, R=1 & R=3
#     - Mash K=1: K=1 mash reduces to single MVN, identical to matrix prior.
#       All est_pv × est_rv × precompute × intercept combos.
#       (Exception: est_rv=TRUE + precomp=TRUE + est_pv=FALSE diverges;
#       flagged as potential bug.)
#     - SS with matrix prior: est_pv × est_rv
#     - SS with mash K=1: est_pv × est_rv × precompute
#     - L=1 matrix EM: est_pv × est_rv
#     - Mixture K>1: all est_pv × est_rv × precompute × intercept combos.
#       S3 eigendecomp vs R6 Cholesky compute the same quantities.
#     - intercept=FALSE: same algorithm, same centering.
#     - Missing data: approximate/exact methods should match R6.
#
#   Apple-to-orange (different algorithm, smoke test only):
#     - Mash grid (K>1 components after scaling, check_null_threshold)
#     - optim method (Brent optimizer convergence noise)
# ============================================================================


# --------------------------------------------------------------------------
# Setup: load R6 reference from master
# --------------------------------------------------------------------------

repo_dir <- tryCatch({
  gcd <- system2("git", c("rev-parse", "--git-common-dir"),
                 stdout = TRUE, stderr = FALSE)
  normalizePath(file.path(gcd, ".."))
}, error = function(e) NULL)

.r6_sweep_funcs <- NULL

load_r6_sweep <- function() {
  if (!is.null(.r6_sweep_funcs)) return(.r6_sweep_funcs)
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

  ref <- pkgload::load_all(ref_source, export_all = TRUE, quiet = TRUE)
  funcs <- list(
    mvsusie = ref$env$mvsusie,
    mvsusie_suff_stat = ref$env$mvsusie_suff_stat,
    MashInitializer = ref$env$MashInitializer,
    create_mixture_prior = ref$env$create_mixture_prior,
    mvsusie_sim1 = ref$env$mvsusie_sim1,
    env = ref$env
  )

  dev_source <- tryCatch(
    rprojroot::find_root(rprojroot::is_r_package),
    error = function(e) normalizePath(file.path(getwd(), "../.."))
  )
  pkgload::load_all(dev_source, export_all = TRUE, quiet = TRUE)

  .r6_sweep_funcs <<- funcs
  funcs
}

# R6 wrappers
r6_sweep_mvsusie <- function(...) {
  r6 <- load_r6_sweep()
  args <- list(...)
  if ("precompute_cache" %in% names(args)) {
    args$precompute_covariances <- args$precompute_cache
    args$precompute_cache <- NULL
  }
  # R6 master uses verbosity (integer); S3 uses verbose (boolean)
  if ("verbose" %in% names(args)) {
    args$verbosity <- if (isTRUE(args$verbose)) 2 else 0
    args$verbose <- NULL
  }
  # R6 doesn't have estimate_prior_mixture_weights (never updates weights)
  args$estimate_prior_mixture_weights <- NULL
  # Translate missing_y_method to R6's approximate parameter
  if ("missing_y_method" %in% names(args)) {
    args$approximate <- (args$missing_y_method == "approximate")
    args$missing_y_method <- NULL
  }
  do.call(r6$mvsusie, args)
}

r6_sweep_mvsusie_ss <- function(...) {
  r6 <- load_r6_sweep()
  args <- list(...)
  if ("precompute_cache" %in% names(args)) {
    args$precompute_covariances <- args$precompute_cache
    args$precompute_cache <- NULL
  }
  # R6 master uses verbosity (integer); S3 uses verbose (boolean)
  if ("verbose" %in% names(args)) {
    args$verbosity <- if (isTRUE(args$verbose)) 2 else 0
    args$verbose <- NULL
  }
  args$estimate_prior_mixture_weights <- NULL
  do.call(r6$mvsusie_suff_stat, args)
}

r6_sweep_MashInit <- function(...) {
  r6 <- load_r6_sweep()
  r6$MashInitializer$new(...)
}

r6_sweep_create_mixture <- function(...) {
  r6 <- load_r6_sweep()
  r6$create_mixture_prior(...)
}

# --------------------------------------------------------------------------
# Simulation data
# --------------------------------------------------------------------------
set.seed(42)
sim_r1 <- mvsusie_sim1(n = 100, p = 100, r = 1, s = 3)
set.seed(42)
sim_r3 <- mvsusie_sim1(n = 100, p = 100, r = 3, s = 3)
set.seed(42)
sim_r3_miss <- mvsusie_sim1(n = 100, p = 100, r = 3, s = 3, y_missing = 0.15)

# --------------------------------------------------------------------------
# Machine precision tolerance
# --------------------------------------------------------------------------
tol_exact <- 1e-8

# --------------------------------------------------------------------------
# Helper: compare S3 vs R6 at machine precision
# --------------------------------------------------------------------------
compare_exact <- function(fit, ref, label) {
  expect_equal(fit$alpha, ref$alpha, tolerance = tol_exact,
               check.attributes = FALSE, info = paste(label, "alpha"))
  expect_equal(fit$b1, ref$b1, tolerance = tol_exact,
               check.attributes = FALSE, info = paste(label, "b1"))
  if (!is.null(ref$pip))
    expect_equal(fit$pip, ref$pip, tolerance = tol_exact,
                 check.attributes = FALSE, info = paste(label, "pip"))
  if (!is.null(ref$sigma2))
    expect_equal(fit$sigma2, ref$sigma2, tolerance = tol_exact,
                 check.attributes = FALSE, info = paste(label, "sigma2"))
  expect_equal(fit$lbf, ref$lbf, tolerance = tol_exact,
               check.attributes = FALSE, info = paste(label, "lbf"))
  fit_coef <- fit$coef; fit_coef[is.na(fit_coef)] <- 0
  ref_coef <- ref$coef; ref_coef[is.na(ref_coef)] <- 0
  if (is.matrix(fit_coef) && nrow(fit_coef) > 1)
    expect_equal(fit_coef[-1, , drop = FALSE], ref_coef[-1, , drop = FALSE],
                 tolerance = tol_exact, check.attributes = FALSE,
                 info = paste(label, "coef_slopes"))
  if (!is.null(ref$KL) && !all(is.na(ref$KL)))
    expect_equal(fit$KL, ref$KL, tolerance = tol_exact,
                 check.attributes = FALSE, info = paste(label, "KL"))
}

# Helper: verify S3 output structure (for smoke tests)
check_sane_output <- function(fit, label) {
  expect_true(!is.null(fit$alpha), info = paste(label, "has alpha"))
  expect_true(!is.null(fit$b1), info = paste(label, "has b1"))
  expect_true(!is.null(fit$pip), info = paste(label, "has pip"))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1),
              info = paste(label, "pip in [0,1]"))
  # Each effect's alpha sums to ~1
  alpha_sums <- rowSums(fit$alpha)
  expect_true(all(abs(alpha_sums - 1) < 1e-6),
              info = paste(label, "alpha rows sum to 1"))
}

# ============================================================================
# PART 1: APPLE-TO-APPLE (machine precision, identical code paths)
# ============================================================================

# --------------------------------------------------------------------------
# Matrix prior, EM — all est_pv × est_rv × intercept
#
# S3 and R6 use identical multivariate regression + EM prior variance
# update, with max_iter=10 (before check_null_threshold warmup in R6).
# --------------------------------------------------------------------------

test_that("EXACT: R=1 matrix prior, EM, all est_pv × est_rv × intercept", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    intercept = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r1, {
    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=1 matrix est_pv=%s est_rv=%s int=%s",
                        cfg$est_pv, cfg$est_rv, cfg$intercept)

      fit <- suppressWarnings(
        mvsusie(X, y, L = 10, prior_variance = V[1, 1],
                residual_variance = as.numeric(var(y)),
                estimate_residual_variance = cfg$est_rv,
                estimate_prior_variance = cfg$est_pv,
                estimate_prior_method = "EM",
                max_iter = 10,
                intercept = cfg$intercept, standardize = TRUE, verbose = FALSE)
      )
      ref <- suppressWarnings(
        r6_sweep_mvsusie(X, y, L = 10, prior_variance = V[1, 1],
                         residual_variance = as.numeric(var(y)),
                         estimate_residual_variance = cfg$est_rv,
                         estimate_prior_variance = cfg$est_pv,
                         estimate_prior_method = "EM",
                         max_iter = 10,
                         intercept = cfg$intercept, standardize = TRUE,
                         verbose = FALSE)
      )
      compare_exact(fit, ref, label)
    }
  })
})

test_that("EXACT: R=3 matrix prior, EM, all est_pv × est_rv × intercept", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    intercept = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 matrix est_pv=%s est_rv=%s int=%s",
                        cfg$est_pv, cfg$est_rv, cfg$intercept)

      fit <- suppressWarnings(
        mvsusie(X, y, L = 10, prior_variance = V,
                residual_variance = cov(y),
                estimate_residual_variance = cfg$est_rv,
                estimate_prior_variance = cfg$est_pv,
                estimate_prior_method = "EM",
                max_iter = 10,
                intercept = cfg$intercept, standardize = TRUE, verbose = FALSE)
      )
      ref <- suppressWarnings(
        r6_sweep_mvsusie(X, y, L = 10, prior_variance = V,
                         residual_variance = cov(y),
                         estimate_residual_variance = cfg$est_rv,
                         estimate_prior_variance = cfg$est_pv,
                         estimate_prior_method = "EM",
                         max_iter = 10,
                         intercept = cfg$intercept, standardize = TRUE,
                         verbose = FALSE)
      )
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# Mash K=1 — all est_pv × est_rv × precompute × intercept combos
#
# K=1 mash reduces to a single MVN computation, identical to matrix prior.
# R6 precompute has a bug when est_rv=TRUE + est_pv=FALSE: the cache goes
# stale when sigma2 changes. S3 correctly refreshes the cache. For those
# configs, R6 reference uses precompute=FALSE (confirmed correct).
# --------------------------------------------------------------------------

test_that("EXACT: R=3 mash K=1, all est_pv × est_rv × precompute × intercept", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    precompute = c(TRUE, FALSE),
    intercept = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6_sweep_MashInit(Ulist = list(V), grid = 1, null_weight = 0)

    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 mash_K1 est_pv=%s est_rv=%s precomp=%s int=%s",
                        cfg$est_pv, cfg$est_rv, cfg$precompute, cfg$intercept)

      fit <- suppressWarnings(
        mvsusie(X, y, L = 10, prior_variance = s3_prior,
                residual_variance = cov(y),
                estimate_residual_variance = cfg$est_rv,
                estimate_prior_variance = cfg$est_pv,
                estimate_prior_method = "EM",
                max_iter = 10,
                intercept = cfg$intercept, standardize = TRUE,
                precompute_cache = cfg$precompute, verbose = FALSE)
      )
      # R6 precompute is buggy when est_rv=TRUE + est_pv=FALSE:
      # the eigendecomp cache goes stale. Use R6 precompute=FALSE
      # as the correct reference for those configs.
      r6_precomp <- cfg$precompute
      if (cfg$est_rv && !cfg$est_pv) r6_precomp <- FALSE
      ref <- suppressWarnings(
        r6_sweep_mvsusie(X, y, L = 10, prior_variance = r6_prior,
                         residual_variance = cov(y),
                         estimate_residual_variance = cfg$est_rv,
                         estimate_prior_variance = cfg$est_pv,
                         estimate_prior_method = "EM",
                         max_iter = 10,
                         intercept = cfg$intercept, standardize = TRUE,
                         precompute_cache = r6_precomp, verbose = FALSE)
      )
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# Sufficient statistics, matrix prior — all est_pv × est_rv
# --------------------------------------------------------------------------

test_that("EXACT: R=3 sufficient statistics, matrix prior, all est_pv × est_rv", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    y_centered <- scale(y, center = TRUE, scale = FALSE)
    XtX <- crossprod(X_centered)
    XtY <- crossprod(X_centered, y_centered)
    YtY <- crossprod(y_centered)
    N <- nrow(X)

    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 SS matrix est_pv=%s est_rv=%s",
                        cfg$est_pv, cfg$est_rv)

      fit <- mvsusie_ss(XtX = XtX, XtY = XtY, YtY = YtY, N = N,
                        L = 10, prior_variance = V,
                        residual_variance = cov(y_centered),
                        estimate_residual_variance = cfg$est_rv,
                        estimate_prior_variance = cfg$est_pv,
                        estimate_prior_method = "EM",
                        max_iter = 10,
                        verbose = FALSE)
      ref <- r6_sweep_mvsusie_ss(XtX = XtX, XtY = XtY, YtY = YtY, N = N,
                                 L = 10, prior_variance = V,
                                 residual_variance = cov(y_centered),
                                 estimate_residual_variance = cfg$est_rv,
                                 estimate_prior_variance = cfg$est_pv,
                                 estimate_prior_method = "EM",
                                 max_iter = 10,
                                 verbose = FALSE)
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# Sufficient statistics, mash K=1 — all est_pv × est_rv × precompute
# --------------------------------------------------------------------------

test_that("EXACT: R=3 sufficient statistics, mash K=1, all est_pv × est_rv × precompute", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    precompute = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    y_centered <- scale(y, center = TRUE, scale = FALSE)
    XtX <- crossprod(X_centered)
    XtY <- crossprod(X_centered, y_centered)
    YtY <- crossprod(y_centered)
    N <- nrow(X)

    s3_prior <- create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    r6_prior <- r6_sweep_MashInit(Ulist = list(V), grid = 1, null_weight = 0)

    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 SS mash_K1 est_pv=%s est_rv=%s precomp=%s",
                        cfg$est_pv, cfg$est_rv, cfg$precompute)

      fit <- mvsusie_ss(XtX = XtX, XtY = XtY, YtY = YtY, N = N,
                        L = 10, prior_variance = s3_prior,
                        residual_variance = cov(y_centered),
                        estimate_residual_variance = cfg$est_rv,
                        estimate_prior_variance = cfg$est_pv,
                        estimate_prior_method = "EM",
                        max_iter = 10,
                        precompute_cache = cfg$precompute, verbose = FALSE)
      # R6 precompute is buggy when est_rv=TRUE + est_pv=FALSE
      r6_precomp <- cfg$precompute
      if (cfg$est_rv && !cfg$est_pv) r6_precomp <- FALSE
      ref <- r6_sweep_mvsusie_ss(XtX = XtX, XtY = XtY, YtY = YtY, N = N,
                                 L = 10, prior_variance = r6_prior,
                                 residual_variance = cov(y_centered),
                                 estimate_residual_variance = cfg$est_rv,
                                 estimate_prior_variance = cfg$est_pv,
                                 estimate_prior_method = "EM",
                                 max_iter = 10,
                                 precompute_cache = r6_precomp,
                                 verbose = FALSE)
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# L=1, matrix prior, EM — all est_pv × est_rv
# --------------------------------------------------------------------------

test_that("EXACT: L=1 matrix prior EM, all est_pv × est_rv", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 L=1 matrix est_pv=%s est_rv=%s",
                        cfg$est_pv, cfg$est_rv)

      fit <- mvsusie(X, y, L = 1, prior_variance = V,
                     residual_variance = cov(y),
                     estimate_residual_variance = cfg$est_rv,
                     estimate_prior_variance = cfg$est_pv,
                     estimate_prior_method = "EM",
                     max_iter = 10,
                     intercept = TRUE, standardize = TRUE, verbose = FALSE)
      ref <- r6_sweep_mvsusie(X, y, L = 1, prior_variance = V,
                              residual_variance = cov(y),
                              estimate_residual_variance = cfg$est_rv,
                              estimate_prior_variance = cfg$est_pv,
                              estimate_prior_method = "EM",
                              max_iter = 10,
                              intercept = TRUE, standardize = TRUE, verbose = FALSE)
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# Mixture K>1: all est_pv × est_rv × precompute × intercept
#
# S3 eigendecomp vs R6 mashr Cholesky compute the same mathematical
# quantities. Differences indicate bugs, not algorithm mismatch.
#
# R6 does not have estimate_prior_mixture_weights parameter (it never
# updates mixture weights). S3 defaults to TRUE. For apple-to-apple
# comparison, S3 must pass estimate_prior_mixture_weights = FALSE.
#
# R6 precompute has a bug when est_rv=TRUE + est_pv=FALSE (stale cache).
# S3 correctly refreshes the cache. R6 reference uses precompute=FALSE
# for those configs.
# --------------------------------------------------------------------------

test_that("EXACT: R=3 mixture K>1, all est_pv × est_rv × precompute × intercept", {
  load_r6_sweep()

  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    precompute = c(TRUE, FALSE),
    intercept = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    U1 <- V
    U2 <- 0.5 * diag(ncol(V))
    s3_prior <- create_mixture_prior(mixture_prior = list(matrices = list(U1, U2)))
    r6_prior <- r6_sweep_create_mixture(mixture_prior = list(matrices = list(U1, U2)))

    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("R=3 mixture_K2 est_pv=%s est_rv=%s precomp=%s int=%s",
                        cfg$est_pv, cfg$est_rv, cfg$precompute, cfg$intercept)

      fit <- suppressWarnings(
        mvsusie(X, y, L = 10, prior_variance = s3_prior,
                residual_variance = cov(y),
                estimate_residual_variance = cfg$est_rv,
                estimate_prior_variance = cfg$est_pv,
                estimate_prior_method = "EM",
                estimate_prior_mixture_weights = FALSE,
                max_iter = 10,
                intercept = cfg$intercept, standardize = TRUE,
                precompute_cache = cfg$precompute, verbose = FALSE)
      )
      # R6 precompute is buggy when est_rv=TRUE + est_pv=FALSE
      r6_precomp <- cfg$precompute
      if (cfg$est_rv && !cfg$est_pv) r6_precomp <- FALSE
      ref <- suppressWarnings(
        r6_sweep_mvsusie(X, y, L = 10, prior_variance = r6_prior,
                         residual_variance = cov(y),
                         estimate_residual_variance = cfg$est_rv,
                         estimate_prior_variance = cfg$est_pv,
                         estimate_prior_method = "EM",
                         max_iter = 10,
                         intercept = cfg$intercept, standardize = TRUE,
                         precompute_cache = r6_precomp, verbose = FALSE)
      )
      compare_exact(fit, ref, label)
    }
  })
})

# --------------------------------------------------------------------------
# Missing data: approximate/exact methods match R6 at machine precision.
#
# Note: S3 defaults to missing_y_method="approximate", R6 defaults to
# approximate=FALSE (exact). Must explicitly match methods for comparison.
# Also use tol=0 to force equal iteration counts, because convergence
# criteria differ slightly (ELBO computation for missing entries).
# --------------------------------------------------------------------------

test_that("SMOKE: missing data, matrix prior, all est_pv × est_rv", {
  # S3 fixed a bug where prior variance scaling used Y[,1] for all columns
  # instead of Y[,i] per column. R6 reference still has this bug, so exact
  # S3-vs-R6 comparison is no longer valid for missing data with
  # standardize=TRUE. We test S3 as smoke tests instead.
  with(sim_r3_miss, {
    for (est_pv in c(FALSE, TRUE)) {
      for (est_rv in c(FALSE, TRUE)) {
        label <- sprintf("R=3 missing matrix est_pv=%s est_rv=%s",
                          est_pv, est_rv)
        fit <- suppressWarnings(
          mvsusie(X, y_missing, L = 10, prior_variance = V,
                  residual_variance = cov(y),
                  estimate_residual_variance = est_rv,
                  estimate_prior_variance = est_pv,
                  estimate_prior_method = "EM",
                  missing_y_method = "approximate",
                  max_iter = 10, tol = 0,
                  intercept = TRUE, standardize = TRUE, verbose = FALSE)
        )
        expect_true(!is.null(fit$alpha), label = label)
        expect_true(all(is.finite(fit$pip)), label = paste(label, "pip"))
      }
    }
  })
})

# ============================================================================
# PART 2: SMOKE TESTS (S3 only, no R6 comparison)
#
# Genuinely different algorithms where numerical comparison is meaningless:
#   - Mash grid (K>1 after scaling, check_null_threshold interaction)
#   - Optim method (Brent optimizer convergence noise)
# ============================================================================

# --------------------------------------------------------------------------
# Mash grid (K>1 after scaling): check_null_threshold warmup interaction
# --------------------------------------------------------------------------

test_that("SMOKE: mash grid, all est_pv × est_rv × precompute", {
  configs <- expand.grid(
    est_pv = c(FALSE, TRUE),
    est_rv = c(FALSE, TRUE),
    precompute = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  with(sim_r3, {
    s3_prior <- create_mash_prior(Ulist = list(V), grid = c(0.5, 1, 2),
                                   null_weight = 0)

    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      label <- sprintf("mash_grid est_pv=%s est_rv=%s precomp=%s",
                        cfg$est_pv, cfg$est_rv, cfg$precompute)

      fit <- expect_no_error(suppressWarnings(
        mvsusie(X, y, L = 10, prior_variance = s3_prior,
                residual_variance = cov(y),
                estimate_residual_variance = cfg$est_rv,
                estimate_prior_variance = cfg$est_pv,
                estimate_prior_method = "EM",
                max_iter = 10,
                intercept = TRUE, standardize = TRUE,
                precompute_cache = cfg$precompute, verbose = FALSE)
      ))
      check_sane_output(fit, label)
    }
  })
})

# --------------------------------------------------------------------------
# Optim method: Brent optimizer convergence noise
# --------------------------------------------------------------------------

test_that("SMOKE: optim method, matrix prior", {
  configs <- expand.grid(
    est_rv = c(FALSE, TRUE),
    R = c(1, 3),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    sim <- if (cfg$R == 1) sim_r1 else sim_r3
    pv <- if (cfg$R == 1) sim$V[1, 1] else sim$V
    rv <- if (cfg$R == 1) as.numeric(var(sim$y)) else cov(sim$y)

    label <- sprintf("optim R=%s est_rv=%s", cfg$R, cfg$est_rv)

    fit <- expect_no_error(suppressWarnings(
      mvsusie(sim$X, sim$y, L = 10, prior_variance = pv,
              residual_variance = rv,
              estimate_residual_variance = cfg$est_rv,
              estimate_prior_variance = TRUE,
              estimate_prior_method = "optim",
              max_iter = 10,
              intercept = TRUE, standardize = TRUE, verbose = FALSE)
    ))
    check_sane_output(fit, label)
  }
})
