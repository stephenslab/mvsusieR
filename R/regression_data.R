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
      private$standardize(center,scale)
      private$.fitted = matrix(0, private$N, private$R) 
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
    remove_from_fitted = function(value) {
      # this is meant for SuSiE model
      # where a good fit is achieved bit by bit like this
      # (along with add_to_fitted method)
      private$.fitted = private$.fitted - value
    },
    add_to_fitted = function(value) {
      private$.fitted = private$.fitted + value
    },
    compute_residual = function() {
      private$residual = private$.Y - private$.fitted 
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
    .fitted = NULL,
    N = NULL,
    J = NULL,
    R = NULL,
    X2t = NULL,
    residual = NULL,
    csd = NULL,
    cm = NULL,
    Y_mean = NULL,
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
      private$X2t = t(private$.X * private$.X)
      private$.d = colSums(t(private$X2t))
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
    fitted = function(value) {
      if (missing(value)) private$.fitted
      else private$denied('fitted')
    },
    XtY = function(value) {
      if (missing(value)) crossprod(private$.X, private$.Y)
      else private$denied('XtY')
    },
    XtR = function(value) {
      if (missing(value)) crossprod(private$.X, private$residual)
      else private$denied('XtR')
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
      private$.fitted = matrix(0, private$J, private$R)
      private$residual = private$.XtY
    },
    compute_Xb = function(b) {
      tcrossprod(private$.XtX,t(b))
    },
    remove_from_fitted = function(value) {
      private$.fitted = private$.fitted - value
    },
    add_to_fitted = function(value) {
      private$.fitted = private$.fitted + value
    },
    compute_residual = function() {
      private$residual = private$.XtY - private$.fitted 
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
    .fitted = NULL,
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
    fitted = function(value) {
      if (missing(value)) private$.fitted
      else private$denied('fitted')
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
