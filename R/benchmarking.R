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
