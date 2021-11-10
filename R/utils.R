# chol decomposition without warning message.
muffled_chol = function (x, ...)
  withCallingHandlers(chol(x, ...),
    warning = function (w) {
      if (grepl("the matrix is either rank-deficient or indefinite",w$message))
        invokeRestart("muffleWarning")
  })

# Invert a symmetric, positive definite square matrix via its Cholesky
# decomposition.
invert_via_chol = function (x) {
  if (all(x == 0))
    return(list(inv = x,rank = 0))
  else
    return(list(inv = chol2inv(muffled_chol(x)),rank = nrow(x)))
}

# Invert SPD via triangular back-fitting.
invert_chol_tri = function (x)
  list(inv = t(backsolve(muffled_chol(x),diag(nrow(x)))),rank = nrow(x))

# Pseudoinverse of matrix.
pseudo_inverse = function (x, tol = sqrt(.Machine$double.eps)) {
  xsvd     <- svd(x)
  Positive <- xsvd$d > max(tol * xsvd$d[1L],0)
  if (all(Positive)) {
    xinv <- xsvd$v %*% (1/xsvd$d * t(xsvd$u))
  } else {
    xinv <- xsvd$v[,Positive,drop = FALSE] %*%
              ((1/xsvd$d[Positive]) * t(xsvd$u[,Positive,drop = FALSE]))
  }
  return(list(inv = xinv,rank = sum(Positive)))
}

# Check if x is diagonal matrix.
isDiagonal = function (x, tol = sqrt(.Machine$double.eps)) {
  if (is.matrix(x)) {
    diag(x) <- rep(0,nrow(x))
    return(all(abs(x) < tol))
  } else
    return(TRUE)
}

# Find trace of diag matrix.
tr = function (m) {
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2]))
    stop("Input to tr() function must be a square matrix")
  return(sum(diag(m),na.rm = TRUE))
}

# Convert a list of matrices to array without losing dimension.
matlist2array = function (l) {
  if (class(l) != "list")
    return(l)
  l = simplify2array(l)
  if (is.null(dim(l)))
    l = array(l,c(1,1,length(l)))
  return(l)
}

# Compute value_j * weight_j / sum(value_j * weight_j).
compute_softmax = function (value, weight, log = TRUE) {
  if (length(value) != length(weight))
    stop("Values and their weights should have equal length")
  if (!log)
    value = log(value)
  mvalue = max(value)
  w = exp(value - mvalue)
  w_weighted = w * weight
  weighted_sum_w = sum(w_weighted)
  return(list(weights = as.vector(w_weighted / weighted_sum_w),
              log_sum = log(weighted_sum_w) + mvalue))
}

# Cannot use "unique" directly here --- for perfectly identical rows
# (by computation) due to possible numerical issues, "unique" and
# "duplicated" function report that they are not identical.
almost.unique <- function (x,  tolerance = sqrt(.Machine$double.eps), ...) {
  if (is.matrix(x)) {
    y <- round(x/tolerance,0)
  } else {
    y <- lapply(1:length(x),function(i) round(x[[i]]/tolerance,0))
  }
  d <- duplicated(y,...)
  if (is.matrix(x))
    return(x[!d,,drop = FALSE])
  else
    return(x[!d])
}

# Duplicated function with a tolerance.
almost.duplicated <- function (x, tolerance = sqrt(.Machine$double.eps), ...) {
  y <- round(x/tolerance, 0)
  return(duplicated(y, ...))
}

# A null progressbar, because currently "progressbar_enabled" feature
# does not work for "progress_bar".
#
#' @importFrom R6 R6Class
null_progress_bar = R6Class("null_progress_bar",
                            public = list(tick = function(...) {}))

# Check if all elements are the same in matrix of J by R, J >> R.
is_mat_common = function (mat)
  nrow(almost.unique(mat)) == 1

# Check if all elements are the same in list.
is_list_common = function (lst)
  length(almost.unique(lst)) == 1

# Check if matrix has constant columns.
is_zero_variance <- function (x) {
  if (length(unique(x)) == 1)
    return(TRUE)
  else
    return(FALSE)
}

# Scale prior matrix.
scale_covariance <- function (mat, sigma) {
  # faster way to implement diag(sigma) %*% mat %*% diag(sigma)
  t(mat * sigma) * sigma
}

# Check if input is numeric matrix.
is_numeric_matrix <- function (X, name) {
  if (!((is.double(X) || is.integer(X)) & is.matrix(X)))
    stop(paste("Input",name,"must be a numeric matrix."))
  if (any(is.na(X)))
    stop(paste("Input",name,"must not contain any missing values."))
  if (any(dim(X) == 0))
    stop(paste("Input",name,"dimension is invalid."))
  return(NULL)
}

#' @title remove duplicated columns in matrix while keeping track of what columns are removed for duplicate with what other column
#' @keywords internal
rm_collinear = function(mat, ...) {
    # `duplicated` will only work for matrix not data frame
    mat = as.matrix(mat)
    dimmat = dim(mat)
    bool_coll = almost.duplicated(mat, MARGIN = 2, ...)
    if (any(bool_coll)) {
        # these are columns to be removed
        rmvd_coll = which(bool_coll)
        # now find columns they are collinear with
        # the idea is, when using fromLast, the previously NOT duplicated columns (FALSE) will now become duplicated (TRUE)
        # then we can find these columns and use them as the columns that has some duplicated associated with them
        bool_with_coll = almost.duplicated(mat, MARGIN = 2, fromLast = TRUE, ...) & !bool_coll
        mat_with_coll = mat[, bool_with_coll, drop = FALSE]
        mat_coll = mat[, bool_coll, drop = FALSE]
        # these are columns with which the removed columns are collinear with
        # `match` will only work for data frames
        assoc_coll = which(bool_with_coll)[match(data.frame(mat_coll), data.frame(mat_with_coll))]
        rmvd_coll = cbind(assoc_coll, rmvd_coll)
        colnames(rmvd_coll) = c('associated', 'removed')
        mat = mat[, !bool_coll, drop = FALSE]
        # now generate index to recover the original
        # ...
    } else {
        rmvd_coll = NULL
    }
    attr(mat, "original_dim") = dimmat
    attr(mat, "collinear_cols") = rmvd_coll
    if (is.null(rmvd_coll)) {
        attr(mat, "collinear_counts") = NULL
    } else {
        attr(mat, "collinear_counts") = table(rmvd_coll[,'associated'])
        options(stringsAsFactors = FALSE)
        attr(mat, "collinear_counts") = cbind(as.integer(names(attr(mat, "collinear_counts"))), attr(mat, "collinear_counts") + 1)
        colnames(attr(mat, "collinear_counts")) = c('associated', 'counts')
        rownames(attr(mat, "collinear_counts")) = NULL
    }
    return(mat)
}

#' @title Reconstruct complete matrix (with duplicates) using stripped matrix and information regarding duplicate pattern in original matrix
#' @keywords internal
reconstruct_coll = function(mat, coll_cols, coll_counts, original_dim, adjust_counts=FALSE, transpose=FALSE) {
    # usage:
    # m = rm_collinear(X1)
    # X2 = reconstruct_coll(m, attr(m, 'collinear_cols'), attr(m, 'collinear_counts'), attr(m, 'original_dim'))
    # sum(X1 - X2) == 0
    get_count = function(counts, idx, adjust_counts) {
        if (!adjust_counts || ! (idx %in% counts[,"associated"])) {
            return (1)
        } else {
            print(idx)
            print(counts[,"counts"][which(counts[,"associated"] == idx)])
            return (counts[,"counts"][which(counts[,"associated"] == idx)])
        }
    }
    vectorise = FALSE
    if (is.vector(mat)) {
        vectorise = TRUE
        mat = matrix(mat, ncol = length(mat), nrow = 1)
    }
    if (transpose && !vectorise) mat = t(mat)
    # create empty matrix to the original scale
    res = matrix(NA, original_dim[1], original_dim[2])
    # first column should always be good,
    # and also duplicated columns always can be found in columns already established
    res[,1] = mat[,1] / get_count(coll_counts, 1, adjust_counts)
    i = 2
    for (j in 2:ncol(res)) {
        if (j %in% coll_cols[,'removed']) {
            # a duplicate column, just add it from before
            j0 = coll_cols[,'associated'][which(coll_cols[,'removed'] == j)]
            res[,j] = res[,j0]
        } else {
            # a new column; have to take it from the next inline in input mat
            res[,j] = mat[,i] / get_count(coll_counts, j, adjust_counts)
            i = i + 1
        }
    }
    if (transpose && !vectorise) res = t(res)
    if (vectorise) res = as.vector(res)
    return(res)
}
