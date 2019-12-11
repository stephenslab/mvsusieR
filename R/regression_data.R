#' @title The regular regression data object
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseData <- R6Class("DenseData",
  public = list(
    initialize = function(X,Y,center=TRUE,scale=TRUE,missing_code=NA) {
      # Check input X.
      if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix"))
        stop("Input X must be a double-precision matrix, or a sparse matrix.")
      if (any(is.na(X))) {
        stop("Input X must not contain missing values.")
      }
      if (any(dim(X) == 0)) stop('X input dimension is invalid.')
      private$.X = X
      private$.X_has_missing = any(is.na(private$.X))
      # FIXME: might want to allow for missing in X later?
      # see gaow/mmbr/#5
      if (private$.X_has_missing)
        stop("Missing data in input matrix X is not allowed at this point.")
      if (is.null(dim(Y))) private$.Y = matrix(Y,length(Y),1)
      else private$.Y = Y
      private$R = ncol(private$.Y)
      private$N = nrow(private$.Y)
      private$J = ncol(X)
      if (is.na(missing_code) || is.null(missing_code)) Y_missing = is.na(Y)
      else Y_missing = (Y == missing_code)
      private$Y_non_missing = !Y_missing
      private$.Y_has_missing = any(Y_missing)
      # a backdoor to set Y missing status always to true, to compare computational routines
      if (is.null(missing_code)) private$.Y_has_missing = TRUE
      if(private$.Y_has_missing){
        private$X_for_Y_missing = array(X, dim = c(private$N, private$J, private$R))
        for(r in 1:private$R){
          private$X_for_Y_missing[Y_missing[,r],,r] = NA
        }
      }
      private$standardize(center,scale)
      private$residual = private$.Y
    },
    compute_Xb = function(b) {
      # tcrossprod(A,B) performs A%*%t(B) but faster
      if(private$.Y_has_missing){
        sapply(1:private$R, function(r) private$X_for_Y_missing[,,r] %*% b[,r])
      }else{
        tcrossprod(private$.X,t(b))
      }
    },
    compute_MXt = function(M) {
      # tcrossprod(A,B) performs A%*%t(B) but faster
      if(private$.Y_has_missing){
        t(sapply(1:private$R, function(r) private$X_for_Y_missing[,,r] %*% M[r,]))
      }else{
        tcrossprod(M, private$.X)
      }
    },
    remove_from_residual = function(value) {
      private$residual = private$residual - value
    },
    add_to_residual = function(value) {
      private$residual = private$residual + value
    },
    compute_residual = function(fitted) {
      private$residual = private$.Y - fitted
    },
    rescale_coef = function(b) {
      coefs = b/private$csd
      if (is.null(dim(coefs))) {
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - sum(private$cm * coefs)
        else intercept = 0
        c(intercept, coefs)
      } else {
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - colSums(private$cm * coefs)
        else intercept = 0
        mat = as.matrix(rbind(intercept, coefs))
        rownames(mat) = NULL
        return(mat)
      }
    },
    # Compute multivariate summary statistics in the presence of missing data
    get_sumstats = function(residual_variances, residual_correlation=NULL) {
      if (is.null(residual_correlation)) {
        if (!is.matrix(residual_variances))
          stop("residual variance has to be a matrix if residual correlation is not specified")
        residual_correlation = cov2cor(residual_variances)
        residual_variances = diag(residual_variances)
      }
      # private$d is either vector or matrix
      bhat = self$XtY/private$d
      bhat[which(is.nan(bhat))] = 0
      if (private$.Y_has_missing) {
        Sigma = sqrt(residual_variances) * t(sqrt(residual_variances) * residual_correlation)
        SVS = list()
        XtX_inv = list()
        for(j in 1:private$J) {
          SVS[[j]] = matrix(NA, private$R, private$R)
          XtX_inv[[j]] = matrix(NA, private$R, private$R)
          for(r1 in 1:private$R){
            for(r2 in r1:private$R){
              common = as.logical(private$Y_non_missing[,r1] * private$Y_non_missing[,r2])
              # if `common` is all FALSE the sum below will return zero
              # FIXME: here X_for_Y_missing, after scaling per column, is not 100% correct on the numerator.
              # not sure what to do for now.
              XtX_inv[[j]][r1,r2] = ifelse(private$d[j,r1]*private$d[j,r2] != 0, sum(private$X_for_Y_missing[common,j,r1] * private$X_for_Y_missing[common,j,r2])/(private$d[j,r1]*private$d[j,r2]), ifelse(r1==r2, 1E6, 0))
              SVS[[j]][r1,r2] = Sigma[r1,r2] * XtX_inv[[j]][r1,r2]
              if (r1 != r2) SVS[[j]][r2,r1] = SVS[[j]][r1,r2]
            }
          }
        }
        is_common_sbhat = is_list_common(SVS)
        if (is_common_sbhat) {
            SVS = SVS[[1]]
        }
        sbhat = NA
      } else {
        sbhat = sqrt(do.call(rbind, lapply(1:length(private$d), function(j) residual_variances / private$d[j])))
        sbhat[which(is.nan(sbhat) | is.infinite(sbhat))] = 1E3
        is_common_sbhat = is_mat_common(sbhat)
        if (is_common_sbhat) SVS = sbhat[1,] * t(residual_correlation * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
        else SVS = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(residual_correlation * sbhat[j,]))
      }
      return(list(svs=SVS, sbhat=sbhat, is_common_sbhat = is_common_sbhat, bhat=bhat))
    }
  ),
  active = list(
    X = function() private$.X,
    Y = function() private$.Y,
    X2_sum = function() private$d,
    XtY = function() {
      if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$X_for_Y_missing[private$Y_non_missing[,r],,r], private$.Y[private$Y_non_missing[,r],r]))
      else crossprod(private$.X, private$.Y)
    },
    XtX = function() {
      if(private$.Y_has_missing){
        sapply(1:private$R, function(r) crossprod(private$X_for_Y_missing[private$Y_non_missing[,r],,r]), simplify = "array")
      }else{
        crossprod(private$.X)
      }
    },
    XtR = function() {
      if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$X_for_Y_missing[private$Y_non_missing[,r],,r], private$residual[private$Y_non_missing[,r],r]))
      else crossprod(private$.X, private$residual)
    },
    n_sample = function() private$N,
    n_condition = function() private$R,
    n_effect = function() private$J,
    X_has_missing = function() private$.X_has_missing,
    Y_has_missing = function() private$.Y_has_missing
  ),
  private = list(
    .X = NULL,
    X_for_Y_missing = NULL,
    .Y = NULL,
    d = NULL,
    N = NULL,
    J = NULL,
    R = NULL,
    residual = NULL,
    csd = NULL,
    cm = NULL,
    Y_mean = NULL,
    Y_non_missing = NULL,
    .Y_has_missing = NULL,
    .X_has_missing = NULL,
    standardize = function(center, scale) {
      # Credit: This is heavily based on code from
      # https://www.r-bloggers.com/a-faster-scale-function/
      # The only change from that code is its treatment of columns with 0 variance.
      # This "safe" version scales those columns by 1 instead of 0.
      if(center){
        # center X
        private$cm = colMeans(private$.X, na.rm=T)
        # center Y
        if (private$R == 1) private$Y_mean = mean(private$.Y, na.rm = TRUE)
        else private$Y_mean = colMeans(private$.Y, na.rm = TRUE)
        private$.Y = t(t(private$.Y) - private$Y_mean)
      }else{
        private$cm = rep(0, length = private$J)
      }
      # scale X
      private$csd = rep(1, length = private$J)
      if (scale) {
        private$csd = colSds(private$.X, center = private$cm)
        private$csd[private$csd==0] = 1
      }
      private$.X = t( (t(private$.X) - private$cm) / private$csd )
      # scale X when Y has missing, and compute colSums(X^2)
      if(private$.Y_has_missing){
        if(center){
          private$cm = colMeans(private$X_for_Y_missing, na.rm=T) # J by R
        }else{
          private$cm = matrix(0, private$J, private$R) # J by R
        }
        private$csd = matrix(1, private$J, private$R)
        for(r in 1:private$R){
          if (scale) {
            private$csd[,r] = colSds(private$X_for_Y_missing[,,r], center = private$cm[,r], na.rm = TRUE)
            private$csd[private$csd[,r]==0, r] = 1
          }
          private$X_for_Y_missing[,,r] = t( (t(private$X_for_Y_missing[,,r]) - private$cm[,r]) / private$csd[,r] )
        }
        # For missing Y, d is a J by R matrix
        private$d = sapply(1:private$R, function(r) colSums(private$X_for_Y_missing[,,r]^2, na.rm = T))
      } else {
        # For non-missing Y, d is a J vector
        private$d = colSums(private$.X ^ 2)
      }
    }
  )
)