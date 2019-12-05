#' @title The regular regression data object
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseData <- R6Class("DenseData",
  public = list(
    initialize = function(X,Y,center=TRUE,scale=TRUE,missing_code=NULL) {
      if (any(dim(X) == 0)) stop('X input dimension is invalid.')
      private$.X = X
      if (is.null(dim(Y))) private$.Y = matrix(Y,length(Y),1)
      else private$.Y = Y
      private$R = ncol(private$.Y)
      private$N = nrow(private$.Y)
      private$J = ncol(X)
      if (is.null(missing_code)) Y_missing = is.na(Y)
      else Y_missing = (Y == missing_code)
      private$Y_non_missing = !Y_missing
      private$.Y_has_missing = any(Y_missing)
      if (!is.null(missing_code) && missing_code=='test') private$.Y_has_missing = TRUE
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
    get_sumstats = function(residual_variances, residual_correlation, alpha = 0) {
      # compute and assign Xty
      Xty = self$XtY
      bhat = Xty/private$d
      bhat[which(is.nan(bhat))] = 0
      if (private$.Y_has_missing) {
        S = lapply(1:private$J, function(j) diag(1/private$d[j,]) * residual_variances)
        sbhat0 = sqrt(do.call(rbind, lapply(1:length(S), function(j) diag(S[[j]]))))
        sbhat0[which(is.nan(sbhat0) | is.infinite(sbhat0))] = 1E6
        sbhat = sbhat0 ^ (1 - alpha)
        is_common_sbhat = is_mat_common(sbhat)
        Sigma = sqrt(residual_variances ^ (1 - alpha)) * t(sqrt(residual_variances ^ (1 - alpha)) * residual_correlation)
        # FIXME: may want to do this in parallel
        for(j in 1:private$J){
          if (alpha != 0) diag(S[[j]]) = diag(S[[j]]) ^ (1 - alpha)
          S[[j]][which(is.nan(S[[j]]) | is.infinite(S[[j]]))] = 1E6
          for(r in 1:(private$R-1)){
            for(p in (r+1):private$R){
              common = as.logical(private$Y_non_missing[,r] * private$Y_non_missing[,p])
              S[[j]][r,p] = Sigma[r,p] * ifelse(private$d[j,r]*private$d[j,p] != 0, sum((private$.X[common,j])^2)/(private$d[j,r]*private$d[j,p]), 0) ^ (1 - alpha)
              S[[j]][p,r] = S[[j]][r,p]
            }
          }
          if (is_common_sbhat) {
            S = S[[1]]
            break
          }
        }
      } else {
        sbhat0 = sqrt(do.call(rbind, lapply(1:length(private$d), function(j) residual_variances / private$d[j])))
        sbhat0[which(is.nan(sbhat0) | is.infinite(sbhat0))] = 1E6
        sbhat = sbhat0 ^ (1 - alpha)
        is_common_sbhat = is_mat_common(sbhat)
        if (is_common_sbhat) S = sbhat[1,] * t(residual_correlation * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
        else S = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(residual_correlation * sbhat[j,]))
      }
      return(list(svs=S, sbhat0=sbhat0, sbhat=sbhat, is_common_sbhat = is_common_sbhat, bhat=bhat))
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
    XtX = function() {
      if(private$.Y_has_missing){
        sapply(1:private$R, function(r) crossprod(private$X_for_Y_missing[private$Y_non_missing[,r],,r]), simplify = "array")
      }else{
        crossprod(private$.X)
      }
    },
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