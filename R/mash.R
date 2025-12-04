#' @title Create mash prior object.
#'
#' @param mixture_prior a list of (weights = vector(), matrices =
#'   list()) where matrices is a list of prior matrices and have same
#'   length as weights.
#'
#' @param R number of traits
#'
#' @param null_weight whether or not to add a weight for null in
#'   single effect models. By default it takes the null weight from
#'   fitted_g if available. Use \code{null_weight = 0} to override this.
#'
#' @param weights_tol Filter out mixture components with weights
#' smaller than \code{weights_tol}.
#'
#' @param max_mixture_len Only keep the top priors by weight so that
#'   the list of mixture prior is of length \code{max_mixture_len}. Use
#'   \code{max_mixture_len = -1} to include all input weights after
#'   filtering by \code{weights_tol}.
#'
#' @param include_indices Post-process input prior to only include
#'   conditions from this indices.
#'
#' @param \dots other parameters, for mvsusieR:::create_cov_canonical
#'
#' @return mash prior object for use with mvsusie() function
#'
#' @details Add details here.
#'
#' @examples
#' # Add examples here.
#'
#' @importFrom stats cov2cor
#' @importFrom stats setNames
#' @importFrom mashr mash
#'
#' @export
#'
create_mixture_prior <- function(mixture_prior, R, null_weight = NULL,
                                 weights_tol = 1e-10,
                                 max_mixture_len = -1, include_indices = NULL, ...) {
  if (sum(c(missing(mixture_prior), missing(R))) != 1) {
    stop("Require exactly one of mixture_prior and R")
  }
  if (is.null(null_weight)) {
    null_weight <- 0
  }
  if (!missing(mixture_prior)) {
    if (!("matrices" %in% names(mixture_prior))) {
      stop("mixture_prior must contain 'matrices'.")
    }
    if (is.null(mixture_prior$weights)) {
      mixture_prior$weights <- rep(
        1 / length(mixture_prior$matrices),
        length(mixture_prior$matrices)
      )
    }
    return(MashInitializer$new(NULL, NULL,
      xUlist = mixture_prior$matrices,
      prior_weights = mixture_prior$weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_conditions = include_indices
    ))
  }
  if (!missing(R)) {
    Ulist <- create_cov_canonical(R, ...)
    Ulist <- lapply(Ulist, function(mat) {
                           rownames(mat) <- include_indices
                           colnames(mat) <- include_indices
                           return(mat)
                         })
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
    return(MashInitializer$new(NULL, NULL,
      xUlist = Ulist,
      prior_weights = weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_conditions = include_indices
    ))
  }
}

#' @title Compute List of Canonical Covariance Matrices
#'
#' @description This function computes canonical covariance matrices
#'   to be provided to mash.
#'
#' @param R Integer specifying the number of conditions.
#'
#' @param singletons \code{TRUE} or \code{FALSE} indicating whether
#'   the singleton matrices are computed.
#'
#' @param hetgrid A vector of numbers between -1 and 1, each
#'   representing the off-diagonal elements of matrices with 1s on the
#'   diagonal. If 0 is included, the identity matrix will be returned
#'   which corresponds to assuming effects are independent across
#'   conditions. IF \code{hetgrid = NULL}, these matrices are not
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
