#' Create a mash prior object
#'
#' Constructs an S3 object of class \code{mash_prior} containing a mixture of
#' prior covariance matrices for multivariate SuSiE.
#'
#' @param Ulist List of unscaled prior covariance matrices (used with grid).
#' @param grid Numeric vector of scaling factors (used with Ulist).
#' @param xUlist List of pre-scaled prior covariance matrices (alternative to
#'   Ulist + grid).
#' @param prior_weights Numeric vector of mixture weights for non-null
#'   components. If NULL, uniform weights are used.
#' @param null_weight Numeric scalar: weight of the null (zero) component.
#' @param weights_tol Filter out components with weight below this threshold.
#' @param null_tol Tolerance for detecting zero matrices.
#' @param top_mixtures Keep only the top components by weight. Use -1 for all.
#' @param include_conditions Indices of conditions to include (subset prior).
#'
#' @return An S3 object of class \code{mash_prior}.
#'
#' @keywords internal
create_mash_prior <- function(Ulist = NULL, grid = NULL, xUlist = NULL,
                               prior_weights = NULL, null_weight = 0,
                               weights_tol = 1e-10, null_tol = 5e-7,
                               top_mixtures = 20,
                               include_conditions = NULL) {
  all_zeros <- vector()

  # Build xUlist from Ulist + grid if not provided directly
  if (is.null(xUlist)) {
    if (is.null(Ulist)) {
      stop("Either xUlist or Ulist have to be non-null")
    }
    for (l in seq_along(Ulist)) {
      if (all(abs(Ulist[[l]]) < null_tol)) {
        stop(paste("Prior covariance", l, "is zero matrix. This is not allowed."))
      }
    }
    if (any(grid <= 0)) {
      stop("grid values should be greater than zero")
    }
    xUlist <- expand_cov(Ulist, grid, usepointmass = TRUE)
  } else {
    # Ensure null component is first
    if (!all(xUlist[[1]] == 0)) {
      xUlist <- c(
        list(null_model = matrix(0, nrow(xUlist[[1]]), ncol(xUlist[[1]]),
          dimnames = list(rownames(xUlist[[1]]), colnames(xUlist[[1]]))
        )),
        xUlist
      )
    }
  }

  # Subset conditions if requested
  if (!is.null(include_conditions)) {
    for (l in seq_along(xUlist)) {
      xUlist[[l]] <- xUlist[[l]][include_conditions, include_conditions]
      if (l > 1) {
        all_zeros[l - 1] <- all(abs(xUlist[[l]]) < null_tol)
      }
    }
  }

  # Set up prior weights for non-null components
  plen <- length(xUlist) - 1
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / plen, plen)
  }
  if (length(prior_weights) != plen) {
    stop(paste(
      "Invalid prior_weights setting: expect length", plen,
      "but input is of length", length(prior_weights)
    ))
  }

  # Filter by weights lower bound (keep null component)
  if (weights_tol > 0) {
    which.comp <- which(prior_weights > weights_tol)
    prior_weights <- prior_weights[which.comp]
    xUlist <- xUlist[c(1, which.comp + 1)]
  }

  # Remove all-zero priors after condition subsetting
  if (length(which(all_zeros)) > 0) {
    which.comp <- which(sapply(
      2:length(xUlist),
      function(l) !all(xUlist[[l]] == 0)
    ))
    prior_weights <- prior_weights[which.comp]
    xUlist <- xUlist[c(1, which.comp + 1)]
  }

  # Keep only top components by weight
  if (top_mixtures > 0 && top_mixtures < length(prior_weights)) {
    which.comp <- head(
      sort(prior_weights, index.return = TRUE, decreasing = TRUE)$ix,
      top_mixtures
    )
    prior_weights <- prior_weights[which.comp]
    xUlist <- xUlist[c(1, which.comp + 1)]
  }

  # Validate all matrices
  u_rows <- vapply(xUlist, nrow, integer(1))
  for (i in seq_along(xUlist)) {
    check_covmat_basics(xUlist[[i]])
    if (!issemidef(xUlist[[i]])) {
      stop(paste("The prior matrices", i, "should be positive semi-definite"))
    }
  }
  if (length(unique(u_rows)) > 1) {
    stop("Ulist contains matrices of different dimensions")
  }

  # Normalize non-null weights
  prior_weights <- prior_weights / sum(prior_weights)

  # Separate null from non-null: xUlist[[1]] is always the null (zero matrix)
  # Store only non-null components in the S3 object
  xUlist_nonnull <- xUlist[-1]  # K non-null components
  K <- length(xUlist_nonnull)

  structure(list(
    null_weight  = null_weight,
    pi           = setNames(prior_weights, names(xUlist_nonnull)),
    xUlist       = xUlist_nonnull,
    xUlist_3d    = matlist2array(xUlist_nonnull),
    xUlist_inv   = NULL,
    xUlist_rank  = NULL,
    n_condition  = nrow(xUlist_nonnull[[1]]),
    n_component  = K
  ), class = "mash_prior")
}

#' Scale prior variance matrices by residual standard deviations
#'
#' @param prior A \code{mash_prior} object.
#' @param sigma Numeric vector of length R (residual standard deviations).
#'
#' @return Modified \code{mash_prior} object with scaled matrices.
#' @keywords internal
scale_prior_variance.mash_prior <- function(prior, sigma) {
  prior$xUlist <- lapply(prior$xUlist, function(U) scale_covariance(U, sigma))
  prior$xUlist_3d <- matlist2array(prior$xUlist)
  # Invalidate cached inverses since matrices changed
  prior$xUlist_inv <- NULL
  prior$xUlist_rank <- NULL
  return(prior)
}

#' Compute pseudo-inverses and ranks for prior matrices
#'
#' Required for EM update of prior variance scalar.
#'
#' @param prior A \code{mash_prior} object.
#'
#' @return Modified \code{mash_prior} object with \code{xUlist_inv} and
#'   \code{xUlist_rank} populated.
#' @keywords internal
compute_prior_inv.mash_prior <- function(prior) {
  my_lapply <- if (requireNamespace("future.apply", quietly = TRUE))
    future.apply::future_lapply else lapply
  inv_list <- my_lapply(prior$xUlist, pseudo_inverse)
  prior$xUlist_inv <- matlist2array(lapply(inv_list, `[[`, "inv"))
  prior$xUlist_rank <- vapply(inv_list, `[[`, numeric(1), "rank")
  return(prior)
}


#' Construct full prior arrays for mashr C++ (prepend null component)
#'
#' @param V_structure_3d R x R x K array of non-null prior components.
#' @param pi_V K-vector of non-null mixture weights.
#' @param null_weight Scalar null component weight.
#' @param V_scalar Scalar to scale the non-null components.
#' @param R Number of conditions.
#'
#' @return List with \code{xUlist_full_3d} (R x R x (K+1)) and
#'   \code{pi_full} ((K+1)-vector).
#' @keywords internal
build_mashr_prior <- function(V_structure, pi_V, null_weight, V_scalar, R) {
  null_mat <- array(0, c(R, R, 1))
  # V_structure can be an R×R×K array or a list of K R×R matrices.
  # Convert list to array on-the-fly (avoids storing the 3d array permanently).
  if (is.list(V_structure)) V_structure <- matlist2array(V_structure)
  V_scaled <- V_structure * V_scalar
  list(
    xUlist_full_3d = abind::abind(null_mat, V_scaled, along = 3),
    pi_full = c(null_weight, pi_V * (1 - null_weight))
  )
}
