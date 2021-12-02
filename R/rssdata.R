# Summary statistics object (storing Z and R).
#
#' @importFrom R6 R6Class
RSSData <- R6Class("RSSData",
  inherit = DenseData,
  portable = FALSE,
  public = list(
    initialize = function(Z, R=NULL, eigenR=NULL, tol) {
      if(any(is.infinite(Z))){
        stop('Z scores contain infinite value.')
      }
      if (is.null(dim(Z))) {
        Z = matrix(Z, length(Z), 1)
      }
      if(is.null(R) & is.null(eigenR)){
        stop('At least one of R and eigen-decomposition of R should be provided.')
      }
      if(!is.null(eigenR)){
        if(any(names(eigenR) != c("values", "vectors"))){
          stop('eigen-decomposition of R must contains 2 elements with the following names: values, vectors')
        }
        if(nrow(eigenR$vectors) != nrow(Z)){
          stop('The dimension of eigenvectors of R does not agree with expected.')
        }
        .eigenR <<- eigenR
      }
      if(!is.null(R)){
        is_numeric_matrix(R, "R")
        # Check input R.
        if (!susieR:::is_symmetric_matrix(R)) {
          stop('R is not a symmetric matrix.')
        }
        if (nrow(R) != nrow(Z)) {
          stop(paste0('The dimension of correlation matrix (',nrow(R),' by ',ncol(R),
                      ') does not agree with expected (',nrow(Z),' by ',nrow(Z),')'))
        }
        .XtX <<- R
      }

      # replace NA in z with 0
      if (anyNA(Z)) {
        warning('NA values in Z-scores are replaced with 0.')
        Z[is.na(Z)] = 0
      }
      .J <<- nrow(Z)
      .R <<- ncol(Z)
      private$check_semi_pd(tol)
      .X <<- t(.eigenvectors[, .eigenvalues !=0]) * .eigenvalues[.eigenvalues != 0] ^ (0.5)
      .Y <<- (t(.eigenvectors[, .eigenvalues != 0]) * .eigenvalues[.eigenvalues != 0] ^ (-0.5)) %*% Z
      .Y_has_missing <<- FALSE
      .X_has_missing <<- FALSE
      .XtY <<- .UUt %*% Z # = Z when Z is in eigen space of R
      .residual <<- .Y
    }
  ),
  active = list(
    # n_sample doesn't mean sample size here, it means the number of non zero eigenvalues
    n_sample = function() sum(.eigenvalues > 0)
  ),
  private = list(
    .UUt = NULL,
    .eigenR = NULL,
    .eigenvectors = NULL,
    .eigenvalues = NULL,
    check_semi_pd = function(tol) {
      if(is.null(private$.eigenR)){
        private$.eigenR <<- eigen(.XtX, symmetric = TRUE)
      }
      private$.eigenR$values[abs(.eigenR$values) < tol] = 0

      if (any(private$.eigenR$values < 0)) {
        private$.eigenR$values[private$.eigenR$values < 0] = 0
        warning('Negative eigenvalues are set to 0.')
      }

      .XtX <<- private$.eigenR$vectors %*% (t(private$.eigenR$vectors) * private$.eigenR$values)
      .csd <<- rep(1, length = .J)
      .d <<- diag(.XtX)
      private$.eigenvectors <<- private$.eigenR$vectors
      private$.eigenvalues <<- private$.eigenR$values
      private$.UUt <<- tcrossprod(private$.eigenvectors[, which(private$.eigenvalues > 0)])
    }
  )
)
