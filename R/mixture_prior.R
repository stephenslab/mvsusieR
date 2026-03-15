#' Create mash prior object
#'
#' Constructs a mixture prior for use with \code{mvsusie()}.
#' Accepts one of three input types:
#' \itemize{
#'   \item \code{fitted_g}: Output from \code{mashr::mash()}, which
#'     provides data-driven mixture weights and covariance matrices.
#'     This is the recommended approach for large R.
#'   \item \code{mixture_prior}: A list with \code{matrices} (list of
#'     covariance matrices) and optional \code{weights}.
#'   \item \code{R}: Number of outcomes, to auto-generate canonical
#'     covariance matrices (singletons + shared effects).
#' }
#'
#' @param fitted_g The \code{fitted_g} element from \code{mashr::mash()}
#'   output. Must contain \code{Ulist}, \code{grid}, \code{pi}, and
#'   \code{usepointmass}. The estimated mixture weights are used directly.
#'
#' @param mixture_prior A list of \code{(weights = vector(), matrices =
#'   list())} where matrices is a list of prior covariance matrices.
#'
#' @param R Number of outcomes. Generates canonical covariance matrices
#'   via \code{create_cov_canonical(R)}.
#'
#' @param null_weight Weight for the null component in single effect
#'   models. For \code{fitted_g}, defaults to the mash-estimated null
#'   weight. Use \code{null_weight = 0} to override.
#'
#' @param weights_tol Filter out mixture components with weights
#'   smaller than \code{weights_tol}.
#'
#' @param max_mixture_len Only keep the top priors by weight so that
#'   the list of mixture prior is of length \code{max_mixture_len}. Use
#'   \code{max_mixture_len = -1} to include all input weights after
#'   filtering by \code{weights_tol}.
#'
#' @param grid Numeric vector of scaling factors for the canonical
#'   covariance matrices (used only with \code{R}). When provided, each
#'   canonical matrix is scaled by each grid value, producing
#'   \code{length(Ulist) * length(grid)} mixture components. When
#'   \code{NULL} (default), unscaled canonical matrices are used directly.
#'
#' @param include_indices Post-process input prior to only include
#'   outcomes at these indices.
#'
#' @param \dots Other parameters passed to
#'   \code{mvsusieR:::create_cov_canonical} (e.g., \code{singletons},
#'   \code{hetgrid}).
#'
#' @return A \code{mash_prior} object for use with \code{mvsusie()}.
#'
#' @examples
#' # Canonical prior for R=3 outcomes:
#' prior <- create_mixture_prior(R = 3)
#'
#' # Data-driven prior from mashr (recommended for large R):
#' # m <- mashr::mash(data, Ulist = c(U.ed, U.c))
#' # prior <- create_mixture_prior(fitted_g = m$fitted_g)
#'
#' @importFrom stats cov2cor
#' @importFrom stats setNames
#' @importFrom stats pnorm
#'
#' @export
#'
create_mixture_prior <- function(mixture_prior, R, null_weight = NULL,
                                 weights_tol = 1e-10,
                                 max_mixture_len = -1, grid = NULL,
                                 fitted_g = NULL,
                                 include_indices = NULL, ...) {
  n_provided <- sum(!missing(mixture_prior), !missing(R), !is.null(fitted_g))
  if (n_provided != 1) {
    stop("Require exactly one of fitted_g, mixture_prior, or R")
  }

  # fitted_g from mashr::mash()
  if (!is.null(fitted_g)) {
    for (item in c("pi", "Ulist", "grid", "usepointmass")) {
      if (!(item %in% names(fitted_g))) {
        stop(paste("Cannot find", item, "in fitted_g input"))
      }
    }
    if (fitted_g$usepointmass) {
      prior_weights <- fitted_g$pi[-1]
      # Inherit null weight from mash unless user overrides
      if (is.null(null_weight)) {
        null_weight <- fitted_g$pi[1]
      }
    } else {
      prior_weights <- fitted_g$pi
      if (is.null(null_weight)) null_weight <- 0
    }
    return(create_mash_prior(
      Ulist = fitted_g$Ulist, grid = fitted_g$grid,
      prior_weights = prior_weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_outcomes = include_indices
    ))
  }

  # mixture_prior (list of matrices + weights)
  if (is.null(null_weight)) null_weight <- 0
  if (!missing(mixture_prior)) {
    if (!is.null(grid)) {
      stop("grid cannot be used with mixture_prior (matrices are already scaled)")
    }
    if (!("matrices" %in% names(mixture_prior))) {
      stop("mixture_prior must contain 'matrices'.")
    }
    if (is.null(mixture_prior$weights)) {
      mixture_prior$weights <- rep(
        1 / length(mixture_prior$matrices),
        length(mixture_prior$matrices)
      )
    }
    return(create_mash_prior(
      xUlist = mixture_prior$matrices,
      prior_weights = mixture_prior$weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_outcomes = include_indices
    ))
  }

  # R (auto-generate canonical covariances)
  if (!missing(R)) {
    Ulist <- create_cov_canonical(R, ...)
    if (!is.null(include_indices)) {
      Ulist <- lapply(Ulist, function(mat) {
        rownames(mat) <- include_indices
        colnames(mat) <- include_indices
        return(mat)
      })
    }
    if (!is.null(grid)) {
      return(create_mash_prior(
        Ulist = Ulist, grid = grid,
        null_weight = null_weight,
        weights_tol = weights_tol,
        top_mixtures = max_mixture_len,
        include_outcomes = include_indices
      ))
    }
    weights <- rep(1 / length(Ulist), length(Ulist))
    weights <- setNames(weights, names(Ulist))
    if (max_mixture_len < length(Ulist) && max_mixture_len > 0) {
      stop(paste0(
        "Automatically generated uniform mixture prior is of ",
        "length ", length(Ulist), " and is greater than currently ",
        "specified max_mixture_len ", max_mixture_len,
        ". Please set max_mixture_len = -1 to allow using all ",
        "of them (although computational speed will suffer)."
      ))
    }
    return(create_mash_prior(
      xUlist = Ulist,
      prior_weights = weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_outcomes = include_indices
    ))
  }
}

#' Compute canonical covariance matrices
#'
#' Generates canonical prior covariance matrices (singletons and
#' heterogeneity matrices) for use in mixture priors.
#'
#' @param R Integer specifying the number of outcomes.
#'
#' @param singletons \code{TRUE} or \code{FALSE} indicating whether
#'   the singleton matrices are computed.
#'
#' @param hetgrid A vector of numbers between -1 and 1, each
#'   representing the off-diagonal elements of matrices with 1s on the
#'   diagonal. If 0 is included, the identity matrix will be returned
#'   which corresponds to assuming effects are independent across
#'   outcomes. If \code{hetgrid = NULL}, these matrices are not
#'   returned.
#'
#' @return A list of canonical covariance matrices.
#'
#' @examples
#' mvsusieR:::create_cov_canonical(3)
#' mvsusieR:::create_cov_canonical(3, singletons = FALSE)
#' mvsusieR:::create_cov_canonical(3, hetgrid = NULL)
#'
#' @keywords internal
#'
create_cov_canonical <- function(R, singletons = TRUE,
                                 hetgrid = seq(0, 1, 0.25)) {
  mats <- list()
  nms <- vector("double")
  s_idx <- 0

  # Singleton matrices.
  if (singletons) {
    for (i in 1:R) {
      mats[[i]] <- matrix(0, R, R)
      mats[[i]][i, i] <- 1
      nms[i] <- paste0("singleton_", i)
    }
    s_idx <- R
  }

  # Heterogeneity matrices.
  if (!is.null(hetgrid)) {
    for (j in 1:length(hetgrid)) {
      mats[[s_idx + j]] <- matrix(1, R, R)
      mats[[s_idx + j]][lower.tri(mats[[s_idx + j]], diag = FALSE)] <- hetgrid[j]
      mats[[s_idx + j]][upper.tri(mats[[s_idx + j]], diag = FALSE)] <- hetgrid[j]
      nms[s_idx + j] <- paste0("shared_", j)
    }
  }
  names(mats) <- nms
  return(mats)
}

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
#' @param include_outcomes Indices of outcomes to include (subset prior).
#'
#' @return An S3 object of class \code{mash_prior}.
#'
#' @keywords internal
create_mash_prior <- function(Ulist = NULL, grid = NULL, xUlist = NULL,
                               prior_weights = NULL, null_weight = 0,
                               weights_tol = 1e-10, null_tol = 5e-7,
                               top_mixtures = -1,
                               include_outcomes = NULL) {

  # --- Build xUlist from Ulist + grid, or use xUlist directly ---
  if (!is.null(Ulist)) {
    if (is.null(grid)) stop("grid is required when Ulist is provided")
    if (any(grid <= 0)) stop("grid values must be positive")
    for (l in seq_along(Ulist))
      check_covmat_basics(Ulist[[l]])
    # expand_cov prepends a null component; strip it
    xUlist <- expand_cov(Ulist, grid, usepointmass = TRUE)[-1]
  } else if (!is.null(xUlist)) {
    for (i in seq_along(xUlist))
      check_covmat_basics(xUlist[[i]])
  } else {
    stop("Either Ulist or xUlist must be provided")
  }

  # --- Default weights (before any filtering) ---
  if (is.null(prior_weights)) prior_weights <- rep(1 / length(xUlist), length(xUlist))
  if (length(prior_weights) != length(xUlist))
    stop(paste("prior_weights length", length(prior_weights),
               "!= number of matrices", length(xUlist)))

  # --- Subset outcomes first (may turn some matrices near-zero) ---
  if (!is.null(include_outcomes)) {
    xUlist <- lapply(xUlist, function(m) m[include_outcomes, include_outcomes])
  }

  # --- Dimension check ---
  dims <- vapply(xUlist, nrow, integer(1))
  if (length(unique(dims)) > 1)
    stop("Prior matrices have different dimensions")

  # --- Remove near-zero matrices (any position) and drop their weights ---
  is_nonzero <- vapply(xUlist,
                       function(m) !all(abs(m) < null_tol), logical(1))
  xUlist <- xUlist[is_nonzero]
  prior_weights <- prior_weights[is_nonzero]
  if (length(xUlist) == 0)
    stop("All prior covariance matrices are zero or near-zero")

  # --- Filter by weight threshold ---
  if (weights_tol > 0) {
    keep <- prior_weights > weights_tol
    xUlist <- xUlist[keep]
    prior_weights <- prior_weights[keep]
  }

  # --- Keep top components by weight ---
  if (top_mixtures > 0 && top_mixtures < length(prior_weights)) {
    top_idx <- order(prior_weights, decreasing = TRUE)[seq_len(top_mixtures)]
    xUlist <- xUlist[top_idx]
    prior_weights <- prior_weights[top_idx]
  }

  # --- Normalize and build S3 object ---
  prior_weights <- prior_weights / sum(prior_weights)
  K <- length(xUlist)

  structure(list(
    null_weight  = null_weight,
    pi           = setNames(prior_weights, names(xUlist)),
    xUlist       = xUlist,
    xUlist_3d    = matlist2array(xUlist),
    xUlist_inv   = NULL,
    xUlist_rank  = NULL,
    n_outcome    = nrow(xUlist[[1]]),
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
#' @param V_structure List of K R x R matrices or R x R x K array.
#' @param pi_V K-vector of non-null mixture weights.
#' @param null_weight Scalar null component weight.
#' @param V_scalar Scalar to scale the non-null components.
#' @param R Number of outcomes.
#'
#' @return List with \code{xUlist_full_3d} (R x R x (K+1)) and
#'   \code{pi_full} ((K+1)-vector).
#' @keywords internal
build_mashr_prior <- function(V_structure, pi_V, null_weight, V_scalar, R) {
  null_mat <- array(0, c(R, R, 1))
  # V_structure can be an R x R x K array or a list of K R x R matrices.
  # Convert list to array on-the-fly (avoids storing the 3d array permanently).
  if (is.list(V_structure)) V_structure <- matlist2array(V_structure)
  V_scaled <- V_structure * V_scalar
  list(
    xUlist_full_3d = abind::abind(null_mat, V_scaled, along = 3),
    pi_full = c(null_weight, pi_V * (1 - null_weight))
  )
}
