# =============================================================================
# HELPER FUNCTIONS FOR REFERENCE PACKAGE COMPARISON
# =============================================================================
#
# These functions compare the S3 mvsusieR implementation against the reference
# R6 implementation (stephenslab/mvsusieR@1971d14) to ensure results are
# identical.
#
# Strategy: Use pkgload to load the reference package into a separate
# environment, then compare outputs field by field.
library(pkgload)
library(rprojroot)

# Reference package details
.mvsusie_ref_repo <- "stephenslab/mvsusieR"
.mvsusie_ref_commit <- "1971d14"

# Cached environments
.mvsusie_ref_env <- NULL
.mvsusie_dev_env <- NULL
.mvsusie_ref_source_path <- NULL

# Get reference package source (download once, cache path)
get_mvsusie_reference_source <- function() {
  if (!is.null(.mvsusie_ref_source_path) &&
      dir.exists(.mvsusie_ref_source_path)) {
    return(.mvsusie_ref_source_path)
  }

  ref_source <- file.path(tempdir(), "mvsusieR_reference_source")

  if (!dir.exists(ref_source)) {
    message("Downloading reference mvsusieR source from GitHub...")

    result <- system(sprintf(
      "git clone -q https://github.com/%s.git %s 2>&1",
      .mvsusie_ref_repo, ref_source
    ), intern = FALSE)

    if (result != 0) {
      stop("Failed to clone reference package")
    }

    result <- system(sprintf(
      "cd %s && git checkout -q %s 2>&1",
      ref_source, .mvsusie_ref_commit
    ), intern = FALSE)

    if (result != 0) {
      stop("Failed to checkout commit ", .mvsusie_ref_commit)
    }

    message("Reference source downloaded")
  }

  .mvsusie_ref_source_path <<- ref_source
  return(ref_source)
}

# Load reference package using pkgload
load_mvsusie_reference_env <- function() {
  if (!is.null(.mvsusie_ref_env)) {
    return(.mvsusie_ref_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required for reference tests.")
  }

  ref_source <- get_mvsusie_reference_source()

  message("Loading reference mvsusieR with pkgload...")
  env <- pkgload::load_all(ref_source, export_all = TRUE, quiet = TRUE)

  .mvsusie_ref_env <<- env
  return(env)
}

# Load development package using pkgload
load_mvsusie_development_env <- function() {
  if (!is.null(.mvsusie_dev_env)) {
    return(.mvsusie_dev_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required for reference tests.")
  }

  dev_source <- tryCatch({
    rprojroot::find_root(rprojroot::is_r_package)
  }, error = function(e) {
    normalizePath(file.path(getwd(), "../.."))
  })

  message("Loading development mvsusieR with pkgload...")
  env <- pkgload::load_all(dev_source, export_all = TRUE, quiet = TRUE)

  .mvsusie_dev_env <<- env
  return(env)
}

# Skip test if reference not available
skip_if_no_mvsusie_reference <- function() {
  tryCatch({
    load_mvsusie_reference_env()
    load_mvsusie_development_env()
  }, error = function(e) {
    skip(paste("Reference comparison not available:", e$message))
  })
}

# Run mvsusie via the R6 path in the reference package
run_reference_mvsusie <- function(ref_env, X, Y, L, prior_variance, ...) {
  ref_env$env$mvsusie(
    X, Y, L = L,
    prior_variance = prior_variance,
    ...
  )
}

# Run mvsusie via the S3 path in the development package
run_development_mvsusie <- function(dev_env, X, Y, L, prior_variance, ...) {
  dev_env$env$mvsusie(
    X, Y, L = L,
    prior_variance = prior_variance,
    ...
  )
}

# Deep comparison of mvsusie objects
#
# @param Y_mean Column means of Y, used to compare fitted values at
#   standardized scale.  The S3 path returns fitted on the original Y scale
#   (X_std %*% b + Y_mean, matching susieR), while the R6 reference returns
#   fitted on the standardized scale (X_std %*% b, no intercept).
#   When Y_mean is supplied the comparison is:
#     dev_fitted - Y_mean  ==  ref_fitted
#   When Y_mean is NULL the fitted comparison is skipped.
expect_equal_mvsusie_objects <- function(dev_obj, ref_obj,
                                          tolerance = 1e-8,
                                          Y_mean = NULL) {
  # ELBO
  expect_equal(dev_obj$elbo, ref_obj$elbo, tolerance = tolerance,
               info = "ELBO values differ")

  # Number of iterations
  expect_equal(dev_obj$niter, ref_obj$niter,
               info = "Number of iterations differs")

  # alpha (posterior inclusion probabilities)
  expect_equal(unname(dev_obj$alpha), unname(ref_obj$alpha),
               tolerance = tolerance,
               info = "alpha (posterior inclusion probabilities) differ")

  # b1 (posterior first moments)
  expect_equal(unname(dev_obj$b1), unname(ref_obj$b1),
               tolerance = tolerance,
               info = "b1 (posterior first moments) differ")

  # b2 (posterior second moments)
  expect_equal(unname(dev_obj$b2), unname(ref_obj$b2),
               tolerance = tolerance,
               info = "b2 (posterior second moments) differ")

  # KL divergence
  expect_equal(dev_obj$KL, ref_obj$KL, tolerance = tolerance,
               info = "KL divergence differs")

  # lbf (log Bayes factors)
  expect_equal(dev_obj$lbf, ref_obj$lbf, tolerance = tolerance,
               info = "lbf (log Bayes factors) differ")

  # pip
  if (!is.null(dev_obj$pip) && !is.null(ref_obj$pip)) {
    expect_equal(unname(dev_obj$pip), unname(ref_obj$pip),
                 tolerance = tolerance,
                 info = "pip differs")
  }

  # coef (compare non-intercept rows; intercept may differ in NA handling)
  if (!is.null(dev_obj$coef) && !is.null(ref_obj$coef)) {
    if (is.null(dim(dev_obj$coef))) {
      expect_equal(unname(dev_obj$coef[-1]), unname(ref_obj$coef[-1]),
                   tolerance = tolerance,
                   info = "coef (effects) differ")
    } else {
      expect_equal(unname(dev_obj$coef[-1, ]), unname(ref_obj$coef[-1, ]),
                   tolerance = tolerance,
                   info = "coef (effects) differ")
    }
  }

  # intercept
  if (!is.null(dev_obj$intercept) && !is.null(ref_obj$intercept)) {
    # Both may have NA for no-intercept models; compare non-NA values
    non_na <- !is.na(dev_obj$intercept) & !is.na(ref_obj$intercept)
    if (any(non_na)) {
      expect_equal(unname(dev_obj$intercept[non_na]),
                   unname(ref_obj$intercept[non_na]),
                   tolerance = tolerance,
                   info = "intercept differs")
    }
  }

  # Fitted values: compare at standardized scale.
  # S3 returns fitted on original Y scale (X_std %*% b + Y_mean, like susieR).
  # R6 returns fitted on standardized scale (X_std %*% b, no intercept).
  # We compare by subtracting Y_mean from the S3 fitted.
  if (!is.null(Y_mean) &&
      !is.null(dev_obj$fitted) && !is.null(ref_obj$fitted)) {
    n <- nrow(dev_obj$fitted)
    R <- ncol(dev_obj$fitted)
    dev_fitted_std <- dev_obj$fitted -
      matrix(Y_mean, n, R, byrow = TRUE)
    expect_equal(unname(dev_fitted_std), unname(ref_obj$fitted),
                 tolerance = tolerance,
                 info = "fitted values differ (compared at standardized scale)")
  }

  # sigma2 (residual variance)
  if (!is.null(dev_obj$sigma2) && !is.null(ref_obj$sigma2)) {
    expect_equal(unname(dev_obj$sigma2), unname(ref_obj$sigma2),
                 tolerance = tolerance,
                 info = "sigma2 (residual variance) differs")
  }

  # NOTE: V is NOT compared directly because the S3 path stores V as a

  # numeric vector of scalar multipliers (length L), while the R6 path stores
  # V as a list of R x R matrices.  The full prior is V[l] * V_structure in
  # S3 vs V[[l]] in R6.  The EM update for scalar V also follows a slightly
  # different trajectory than the matrix EM update, so estimated V values can
  # differ even when all other quantities match closely.

  # lbf_variable (per-variable log Bayes factors)
  if (!is.null(dev_obj$lbf_variable) && !is.null(ref_obj$lbf_variable)) {
    expect_equal(unname(dev_obj$lbf_variable), unname(ref_obj$lbf_variable),
                 tolerance = tolerance,
                 info = "lbf_variable differs")
  }

  # Credible sets
  if (!is.null(dev_obj$sets) && !is.null(ref_obj$sets)) {
    expect_equal(dev_obj$sets$cs, ref_obj$sets$cs,
                 info = "Credible sets differ")
    # Compare purity if both have it
    if (!is.null(dev_obj$sets$purity) && !is.null(ref_obj$sets$purity)) {
      expect_equal(unname(as.matrix(dev_obj$sets$purity)),
                   unname(as.matrix(ref_obj$sets$purity)),
                   tolerance = tolerance,
                   info = "Credible set purity differs")
    }
  }

  # X_column_scale_factors
  if (!is.null(dev_obj$X_column_scale_factors) &&
      !is.null(ref_obj$X_column_scale_factors)) {
    expect_equal(unname(dev_obj$X_column_scale_factors),
                 unname(ref_obj$X_column_scale_factors),
                 tolerance = tolerance,
                 info = "X_column_scale_factors differ")
  }

  # b1_rescaled: compare non-intercept rows only.
  # S3 stores intercept row as 0, R6 uses NA; dimensions may also differ.
  # For L x (J+1) x R arrays, skip row 1 (intercept) of the 2nd dimension.
  if (!is.null(dev_obj$b1_rescaled) && !is.null(ref_obj$b1_rescaled) &&
      identical(dim(dev_obj$b1_rescaled), dim(ref_obj$b1_rescaled))) {
    if (length(dim(dev_obj$b1_rescaled)) == 3) {
      # L x (J+1) x R array: compare effects rows [, -1, ]
      dev_eff <- dev_obj$b1_rescaled[, -1, , drop = FALSE]
      ref_eff <- ref_obj$b1_rescaled[, -1, , drop = FALSE]
    } else if (is.matrix(dev_obj$b1_rescaled)) {
      # L x (J+1) matrix (R=1): compare effects rows [, -1]
      dev_eff <- dev_obj$b1_rescaled[, -1, drop = FALSE]
      ref_eff <- ref_obj$b1_rescaled[, -1, drop = FALSE]
    } else {
      dev_eff <- dev_obj$b1_rescaled
      ref_eff <- ref_obj$b1_rescaled
    }
    expect_equal(unname(dev_eff), unname(ref_eff),
                 tolerance = tolerance,
                 info = "b1_rescaled (effects) differs")
  }

  invisible(TRUE)
}
