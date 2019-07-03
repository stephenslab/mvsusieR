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
    X_has_missing = function() {
      any(is.na(private$.X))
    },
    Y_has_missing = function() {
      private$.Y_has_missing
    },
    rescale_coef = function(b) {
      coefs = b/private$csd
      if (is.null(dim(coefs))) {
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - sum(private$cm * coefs)
        else intercept = 0
        c(intercept, coefs)
      } else {
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - colSums(private$cm * coefs)
        else intercept = c(0,0)
        mat = as.matrix(rbind(intercept, coefs))
        rownames(mat) = NULL
        return(mat)
      }
    }
  ),
  private = list(
    .X = NULL,
    .Y = NULL,
    .d = NULL,
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
      if (!private$.Y_has_missing) private$.d = colSums(private$.X ^ 2)
      else private$.d = sapply(1:private$R, function(r) colSums(private$.X[private$Y_non_missing[,r],]^2))
    },
    denied = function(v) stop(paste0('$', v, ' is read-only'), call. = FALSE)
  ),
  active = list(
    X = function(value) {
      if (missing(value)) private$.X
      else private$denied('X')
    },
    Y = function(value) {
      if (missing(value)) private$.Y
      else private$denied('Y')
    },
    d = function(value) {
      if (missing(value)) private$.d
      else private$denied('d')
    },
    XtY = function(value) {
      if (missing(value)) {
        if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$.X[private$Y_non_missing[,r],], private$.Y[private$Y_non_missing[,r],r]))
        else crossprod(private$.X, private$.Y)
      } else {
        private$denied('XtY')
      }
    },
    XtR = function(value) {
      if (missing(value)) {
        if (private$.Y_has_missing) sapply(1:private$R, function(r) crossprod(private$.X[private$Y_non_missing[,r],], private$residual[private$Y_non_missing[,r],r]))
        else crossprod(private$.X, private$residual)
      } else {
        private$denied('XtR')
      }
    },
    n_sample = function(value) private$N,
    n_condition = function(value) private$R,
    n_effect = function(value) private$J
  )
)

#' @title Sufficient statistics object
# XtX, XtY, YtY and N
# X and Y have to be centered
# before computing these sufficient statistics
#' @importFrom R6 R6Class
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
  private = list(
    .XtX = NULL,
    .XtY = NULL,
    .YtY = NULL,
    .d = NULL,
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
      private$.d = diag(private$.XtX)
    }
  ),
  active = list(
    XtX = function(value) {
      if (missing(value)) private$.XtX
      else private$.XtX = value
    },
    XtY = function(value) {
      if (missing(value)) private$.XtY
      else private$.XtY = value
    },
    YtY = function(value) {
      if (missing(value)) private$.YtY
      else private$.YtY = value
    },
    d = function(value) {
      if (missing(value)) private$.d
      else private$denied('d')
    },
    XtR = function(value) {
      if (missing(value)) private$residual
      else private$denied('XtR')
    },
    n_sample = function(value) private$N,
    n_condition = function(value) private$R,
    n_effect = function(value) private$J
  )
)