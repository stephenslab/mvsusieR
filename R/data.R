# The regular regression data
DenseData <- R6Class("DenseData",
  public = list(
    X = NULL,
    Y = NULL,
    csd = NULL,
    cm = NULL,
    initialize = function(X,Y,center=TRUE,scale=TRUE) {
      self$X = X
      self$Y = Y
      if (is.null(dim(Y))) private$K = 1
      else private$K = ncol(Y)
      private$N = nrow(Y)
      private$J = ncol(X)
      private$standardize(center,scale)
    },
    get_n_sample = function() private$N,
    get_n_condition = function() private$K,
    get_n_effect = function() private$J
  ),
  private = list(
    N = NULL,
    J = NULL,
    K = NULL,
    d = NULL,
    X2t = NULL,
    standardize = function(center, scale) {
      # Credit: This is heavily based on code from
      # https://www.r-bloggers.com/a-faster-scale-function/
      # The only change from that code is its treatment of columns with 0 variance.
      # This "safe" version scales those columns by 1 instead of 0.
      self$cm = colMeans(self$X, na.rm = TRUE)
      if (scale) {
        self$csd = matrixStats::colSds(self$X, center = self$cm)
        self$csd[self$csd==0] = 1
      } else {
        # just divide by 1 if not
        self$csd = rep(1, length = private$J)
      }
      if (!center) {
        # just subtract 0
        self$cm = rep(0, length = private$J)
      }
      self$X = t( (t(self$X) - self$cm) / self$csd )
      private$X2t = t(self$X * self$X)
      private$d = colSums(t(private$X2t))
    }
  )
)

# Sufficient statistics
# XtX, XtY, YtY and N
# X and Y have to be centered
# before computing these sufficient statistics
SSData <- R6Class("SSData",
  public = list(
    XtX = NULL,
    XtY = NULL,
    YtY = NULL,
    csd = NULL,
    initialize = function(XtX,XtY,YtY,N,scale=TRUE) {
      self$XtX = XtX
      self$XtY = XtY
      self$YtY = YtY
      private$N = N
      private$J = nrow(XtX)
      if (is.null(dim(YtY))) private$K = 1
      else private$K = nrow(YtY)
      private$standardize(scale)
    },
    get_n_sample = function() private$N,
    get_n_condition = function() private$K,
    get_n_effect = function() private$J
  ),
  private = list(
    N = NULL,
    J = NULL,
    K = NULL,
    d = NULL,
    standardize = function(scale) {
      if (scale) {
          d = diag(self$XtX)
          self$csd = sqrt(d/(private$N-1))
          self$csd[self$csd == 0] = 1
          self$XtX = (1/self$csd) * t((1/self$csd) * self$XtX)
          self$XtY = (1/self$csd) * self$XtY
      } else {
          self$csd = rep(1, length = private$J)
      }
      private$d = diag(self$XtX)
    }
  )
)