# Utility functions for mvsusieR.
#
# Includes matrix operations, numerical helpers, benchmarking tools,
# and lfsr computation functions.

# chol decomposition without warning message.
muffled_chol <- function(x, ...) {
  withCallingHandlers(chol(x, ...),
    warning = function(w) {
      if (grepl("the matrix is either rank-deficient or indefinite", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Invert a symmetric, positive definite square matrix via its Cholesky
# decomposition.
invert_via_chol <- function(x) {
  if (all(x == 0)) {
    return(list(inv = x, rank = 0))
  } else {
    return(list(inv = chol2inv(muffled_chol(x)), rank = nrow(x)))
  }
}

# Invert SPD via triangular back-fitting.
invert_chol_tri <- function(x) {
  list(inv = t(backsolve(muffled_chol(x), diag(nrow(x)))), rank = nrow(x))
}

# Pseudoinverse of matrix.
pseudo_inverse <- function(x, tol = sqrt(.Machine$double.eps)) {
  xsvd <- svd(x)
  Positive <- xsvd$d > max(tol * xsvd$d[1L], 0)
  if (all(Positive)) {
    xinv <- xsvd$v %*% (1 / xsvd$d * t(xsvd$u))
  } else {
    xinv <- xsvd$v[, Positive, drop = FALSE] %*%
      ((1 / xsvd$d[Positive]) * t(xsvd$u[, Positive, drop = FALSE]))
  }
  return(list(inv = xinv, rank = sum(Positive)))
}

# Check if x is diagonal matrix.
isDiagonal <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (is.matrix(x)) {
    diag(x) <- rep(0, nrow(x))
    return(all(abs(x) < tol))
  } else {
    return(TRUE)
  }
}

# Find trace of diag matrix.
tr <- function(m) {
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) {
    stop("Input to tr() function must be a square matrix")
  }
  return(sum(diag(m), na.rm = TRUE))
}

# Convert a list of matrices to array without losing dimension.
matlist2array <- function(l) {
  if (!inherits(l, "list")) {
    return(l)
  }
  l <- simplify2array(l)
  if (is.null(dim(l))) {
    l <- array(l, c(1, 1, length(l)))
  }
  return(l)
}

# Cannot use "unique" directly here --- for perfectly identical rows
# (by computation) due to possible numerical issues, "unique" and
# "duplicated" function report that they are not identical.
almost.unique <- function(x, tolerance = sqrt(.Machine$double.eps), ...) {
  if (is.matrix(x)) {
    y <- round(x / tolerance, 0)
  } else {
    y <- lapply(1:length(x), function(i) round(x[[i]] / tolerance, 0))
  }
  d <- duplicated(y, ...)
  if (is.matrix(x)) {
    return(x[!d, , drop = FALSE])
  } else {
    return(x[!d])
  }
}

# Duplicated function with a tolerance.
almost.duplicated <- function(x, tolerance = sqrt(.Machine$double.eps), ...) {
  y <- round(x / tolerance, 0)
  return(duplicated(y, ...))
}

# Check if all elements are the same in matrix of J by R, J >> R.
is_mat_common <- function(mat) {
  nrow(almost.unique(mat)) == 1
}

# Check if all elements are the same in list.
is_list_common <- function(lst) {
  length(almost.unique(lst)) == 1
}

# Check if matrix has constant columns.
is_zero_variance <- function(x) {
  if (length(unique(x)) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Faster way to implement diag(sigma) %*% mat %*% diag(sigma)
scale_covariance <- function(mat, sigma) {
  t(mat * sigma) * sigma
}

# Check if input is numeric matrix.
is_numeric_matrix <- function(X, name) {
  if (!((is.double(X) || is.integer(X)) & is.matrix(X))) {
    stop(paste("Input", name, "must be a numeric matrix."))
  }
  if (anyNA(X)) {
    stop(paste("Input", name, "must not contain any missing values."))
  }
  if (any(dim(X) == 0)) {
    stop(paste("Input", name, "dimension is invalid."))
  }
  return(NULL)
}

#' @title Computes the z-scores (t-statistics) for association
#'   between Y and each column of X.
#'
#' @param X N by J matrix of covariates.
#'
#' @param Y Vector of length N, or N by R matrix of response
#'   variables.
#'
#' @param center If \code{center = TRUE}, center X and Y.
#' 
#' @param scale If \code{scale = TRUE}, scale X and Y.
#'
#' @return A matrix of z-scores.
#' 
#' @importFrom susieR univariate_regression
#'
#' @export
#' 
calc_z = function (X, Y, center = FALSE, scale = FALSE) {
  univariate_z = function(X,Y,center,scale) {
    out = univariate_regression(X,Y,center = center,scale = scale)
    return(out$betahat/out$sebetahat)
  }
  if (is.null(dim(Y)))
    return(univariate_z(X,Y,center,scale))
  else
    return(do.call(cbind,lapply(1:ncol(Y),
                                function(i) univariate_z(X,Y[,i],
                                                         center = center,
                                                         scale = scale))))
}

# Weighted log-sum-exp (softmax).
#
# Given log-values and weights, compute the normalized weights
# (softmax) and the log of their weighted sum.
#
# @param value A numeric vector of log-values (or values if log = FALSE).
# @param weight A numeric vector of weights (same length as value).
# @param log If TRUE (default), value is already on log scale.
#
# @return A list with \code{weights} (normalized) and \code{log_sum}.
#
# @keywords internal
compute_softmax <- function(value, weight, log = TRUE) {
  if (length(value) != length(weight))
    stop("Values and their weights should have equal length")
  if (!log) value <- log(value)
  mvalue <- max(value)
  w <- exp(value - mvalue)
  w_weighted <- w * weight
  weighted_sum_w <- sum(w_weighted)
  list(weights = as.vector(w_weighted / weighted_sum_w),
       log_sum = log(weighted_sum_w) + mvalue)
}

# Remove duplicated columns in the matrix while keeping track of what
# columns are removed for duplicates with what other columns. This
# function is only used for evaluation purposes.
rm_collinear <- function(mat, ...) {
  # "duplicated" will only work for matrix, not data frame.
  mat <- as.matrix(mat)
  dimmat <- dim(mat)
  bool_coll <- almost.duplicated(mat, MARGIN = 2, ...)
  if (any(bool_coll)) {
    # Get columns to be removed.
    rmvd_coll <- which(bool_coll)

    # Now find columns they are collinear with. The idea is, when
    # using fromLast = TRUE, the previously NOT duplicated columns (FALSE)
    # will now become duplicated (TRUE) then we can find these columns
    # and use them as the columns that has some duplicated associated
    # with them.
    bool_with_coll <-
      almost.duplicated(mat, MARGIN = 2, fromLast = TRUE, ...) & !bool_coll
    mat_with_coll <- mat[, bool_with_coll, drop = FALSE]
    mat_coll <- mat[, bool_coll, drop = FALSE]

    # These are columns with which the removed columns are collinear with
    # "match"; will only work for data frames.
    assoc_coll <- which(bool_with_coll)[match(
      data.frame(mat_coll),
      data.frame(mat_with_coll)
    )]
    rmvd_coll <- cbind(assoc_coll, rmvd_coll)
    colnames(rmvd_coll) <- c("associated", "removed")
    mat <- mat[, !bool_coll, drop = FALSE]

    # Now generate index to recover the original.
  } else {
    rmvd_coll <- NULL
  }
  attr(mat, "original_dim") <- dimmat
  attr(mat, "collinear_cols") <- rmvd_coll
  if (is.null(rmvd_coll)) {
    attr(mat, "collinear_counts") <- NULL
  } else {
    attr(mat, "collinear_counts") <- table(rmvd_coll[, "associated"])
    options(stringsAsFactors = FALSE)
    attr(mat, "collinear_counts") <-
      cbind(
        as.integer(names(attr(mat, "collinear_counts"))),
        attr(mat, "collinear_counts") + 1
      )
    colnames(attr(mat, "collinear_counts")) <- c("associated", "counts")
    rownames(attr(mat, "collinear_counts")) <- NULL
  }

  return(mat)
}

# Reconstruct complete matrix (with duplicates) using stripped matrix
# and information regarding duplicate pattern in original matrix.
#
# example usage:
#
#   m = rm_collinear(X1)
#   X2 = reconstruct_coll(m,attr(m,"collinear_cols"),
#                         attr(m,"collinear_counts"),
#                          attr(m,"original_dim"))
#   sum(X1 - X2) == 0
#
reconstruct_coll <- function(mat, coll_cols, coll_counts, original_dim,
                             adjust_counts = FALSE, transpose = FALSE) {
  get_count <- function(counts, idx, adjust_counts) {
    if (!adjust_counts || !(idx %in% counts[, "associated"])) {
      return(1)
    } else {
      print(idx)
      print(counts[, "counts"][which(counts[, "associated"] == idx)])
      return(counts[, "counts"][which(counts[, "associated"] == idx)])
    }
  }

  vectorise <- FALSE
  if (is.vector(mat)) {
    vectorise <- TRUE
    mat <- matrix(mat, ncol = length(mat), nrow = 1)
  }

  if (transpose && !vectorise) {
    mat <- t(mat)
  }

  # Create empty matrix to the original scale.
  res <- matrix(as.numeric(NA), original_dim[1], original_dim[2])

  # First column should always be good, and also duplicated columns
  # always can be found in columns already established.
  res[, 1] <- mat[, 1] / get_count(coll_counts, 1, adjust_counts)
  i <- 2
  for (j in 2:ncol(res)) {
    if (j %in% coll_cols[, "removed"]) {
      # A duplicate column, just add it from before.
      j0 <- coll_cols[, "associated"][which(coll_cols[, "removed"] == j)]
      res[, j] <- res[, j0]
    } else {
      # A new column; have to take it from the next inline in input mat.
      res[, j] <- mat[, i] / get_count(coll_counts, j, adjust_counts)
      i <- i + 1
    }
  }

  if (transpose && !vectorise) {
    res <- t(res)
  }
  if (vectorise) {
    res <- as.vector(res)
  }

  return(res)
}

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
    return(create_mash_prior(
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
    return(create_mash_prior(
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