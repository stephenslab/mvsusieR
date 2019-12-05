#' @title The regular regression data object
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseData <- R6Class("DenseData",
  public = list(
    initialize = function(X,Y,center=TRUE,scale=TRUE) {
      private$.X = X
      if (is.null(dim(Y))) private$.Y = matrix(Y,length(Y),1)
      else private$.Y = Y
      if (any(dim(X) == 0)) stop('X input dimension is invalid.')
      private$R = ncol(private$.Y)
      private$N = nrow(private$.Y)
      private$J = ncol(private$.X)
      Y_missing = is.na(Y)
      private$Y_non_missing = !Y_missing
      private$.Y_has_missing = any(Y_missing)
      private$standardize(center,scale)
      private$residual = private$.Y
    },
    compute_Xb = function(b) {
      # tcrossprod(A,B) performs A%*%t(B) but faster
      tcrossprod(private$.X,t(b))
    },
    compute_MXt = function(M) {
      # tcrossprod(A,B) performs A%*%t(B) but faster
      tcrossprod(M,private$.X)
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
    }
  ),
  active = list(
    X = function() private$.X,
    Y = function() private$.Y,
    X2_sum = function() private$d,
    XtY = function() {
      if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$.X[private$Y_non_missing[,r],], private$.Y[private$Y_non_missing[,r],r]))
      else crossprod(private$.X, private$.Y)
    },
    XtX = function() crossprod(private$.X, private$.X),
    XtR = function() {
      if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$.X[private$Y_non_missing[,r],], private$residual[private$Y_non_missing[,r],r]))
      else crossprod(private$.X, private$residual)
    },
    n_sample = function() private$N,
    n_condition = function() private$R,
    n_effect = function() private$J,
    X_has_missing = function() any(is.na(private$.X)),
    Y_has_missing = function() private$.Y_has_missing
  ),
  private = list(
    .X = NULL,
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
    standardize = function(center, scale) {
      # Credit: This is heavily based on code from
      # https://www.r-bloggers.com/a-faster-scale-function/
      # The only change from that code is its treatment of columns with 0 variance.
      # This "safe" version scales those columns by 1 instead of 0.
      private$cm = colMeans(private$.X, na.rm = TRUE)
      if (scale) {
        private$csd = colSds(private$.X, center = private$cm)
        private$csd[private$csd==0] = 1
      } else {
        # just divide by 1 if not
        private$csd = rep(1, length = private$J)
      }
      if (!center) {
        # just subtract 0
        private$cm = rep(0, length = private$J)
      } else {
        if (private$R == 1) private$Y_mean = mean(private$.Y, na.rm = TRUE)
        else private$Y_mean = colMeans(private$.Y, na.rm = TRUE)
        private$.Y = private$.Y - private$Y_mean
      }
      private$.X = t( (t(private$.X) - private$cm) / private$csd )
      # For non-missing Y, d is a J vector
      # For missing Y, d is a J by R matrix
      if (!private$.Y_has_missing) private$d = colSums(private$.X ^ 2)
      else private$d = sapply(1:private$R, function(r) colSums(private$.X[private$Y_non_missing[,r],]^2))
    }
  )
)

#' @title Sufficient statistics object
# XtX, XtY, YtY and N
# X and Y have to be centered
# before computing these sufficient statistics
#' @importFrom R6 R6Class
#' @keywords internal
SSData <- R6Class("SSData",
  public = list(
    initialize = function(XtX,XtY,YtY,N,scale=TRUE) {
      private$.XtX = XtX
      private$.XtY = XtY
      private$.YtY = YtY
      private$N = N
      private$J = nrow(XtX)
      if (is.null(dim(YtY))) private$R = 1
      else private$R = nrow(YtY)
      private$standardize(scale)
      private$residual = private$.XtY
    },
    compute_Xb = function(b) {
      tcrossprod(private$.XtX,t(b))
    },
    remove_from_residual = function(value) {
      private$residual = private$residual - value
    },
    add_to_residual = function(value) {
      private$residual = private$residual + value
    },
    compute_residual = function(fitted) {
      private$residual = private$.XtY - fitted
    },
    rescale_coef = function(b) {
      c(0, b/private$csd)
    }
  ),
  active = list(
    XtX = function() private$.XtX,
    XtY = function() private$.XtY,
    YtY = function() private$.YtY,
    X2_sum = function() private$d,
    XtR = function() private$residual,
    n_sample = function() private$N,
    n_condition = function() private$R,
    n_effect = function() private$J
  ),
  private = list(
    .XtX = NULL,
    .XtY = NULL,
    .YtY = NULL,
    d = NULL,
    N = NULL,
    J = NULL,
    R = NULL,
    residual = NULL,
    csd = NULL,
    standardize = function(scale) {
      if (scale) {
          d = diag(private$.XtX)
          private$csd = sqrt(d/(private$N-1))
          private$csd[private$csd == 0] = 1
          private$.XtX = (1/private$csd) * t((1/private$csd) * private$.XtX)
          private$.XtY = (1/private$csd) * private$.XtY
      } else {
          private$csd = rep(1, length = private$J)
      }
      private$d = diag(private$.XtX)
    }
  )
)

#' @title Compute multivariate summary statistics in the presence of missing data
#' @keywords internal
get_sumstats_missing_data = function(X, Y, residual_variances, residual_correlation, alpha) {
  J = ncol(X)
  R = ncol(Y)
  M = !is.na(Y)
  # this is same as DenseData$new()$Xty
  # recomputed here for clarity
  Xty = sapply(1:R, function(r) crossprod(X[M[,r],], Y[M[,r],r]) ) # J by R
  # this is same as DenseData$new()$X2_sum
  # recomputed here for clarity
  X2 = sapply(1:R, function(r) colSums(X[M[,r],]^2 )) # J by R
  bhat = Xty/X2
  S = lapply(1:J, function(j) diag(1/X2[j,]) * residual_variances)
  sbhat0 = sqrt(do.call(rbind, lapply(1:length(S), function(j) diag(S[[j]]))))
  bhat[which(is.nan(bhat))] = 0
  sbhat0[which(is.nan(sbhat0) | is.infinite(sbhat0))] = 1E6
  Sigma = sqrt(residual_variances ^ (1 - alpha)) * t(sqrt(residual_variances ^ (1 - alpha)) * residual_correlation)
  # this is for MASH EE and EZ model
  if (alpha != 0) {
    for (j in 1:J) diag(S[[j]]) = diag(S[[j]]) ^ (1 - alpha)
  }
  # FIXME: may want to do this in parallel
  for(j in 1:J){
    S[[j]][which(is.nan(S[[j]]) | is.infinite(S[[j]]))] = 1E6
    for(r in 1:(R-1)){
      for(d in (r+1):R){
        common = as.logical(M[,r] * M[,d])
        S[[j]][r,d] = Sigma[r,d] * ifelse(X2[j,r]*X2[j,d] != 0, sum((X[common,j])^2)/(X2[j,r]*X2[j,d]), 0) ^ (1 - alpha)
        S[[j]][d,r] = S[[j]][r,d]
      }
    }
    if (alpha == 1) {
      S = S[[1]]
      break
    }
  }
  return(list(svs=S, sbhat0=sbhat0, bhat=bhat))
}