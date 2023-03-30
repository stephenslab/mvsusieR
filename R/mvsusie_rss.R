#' @rdname mvsusie
#'
#' @param Z J x R matrix of z-scores.
#' 
#' @param R J x J LD matrix.
#' 
#' @param N sample size
#' 
#' @param Bhat Alternative summary data giving the estimated effects
#'   (J X R matrix). This, together with \code{Shat}, may be
#'   provided instead of \code{Z}.
#'
#' @param Shat Alternative summary data giving the standard errors of
#'   the estimated effects (J X R matrix). This, together with
#'   \code{Bhat}, may be provided instead of \code{Z}.
#'
#' @param varY The sample covariance of Y, defined as \eqn{Y'Y/(N-1)}.
#'   When the sample covariance is not provided, the coefficients
#'   (returned from \code{coef}) are computed on the
#'   \dQuote{standardized} X, y scale.
#'
#' @param prior_variance Can be either (1) a vector of length L, or a
#'   scalar, for scaled prior variance when Y is univariate (which
#'   should then be equivalent to \code{\link[susieR]{susie}}); or (2) a
#'   matrix for a simple multivariate regression; or (3) a mixture prior 
#'   from \code{\link{create_mixture_prior}}.
#'   
#' @param residual_variance The residual variance
#'
#' @param \dots Additional arguments passed to
#'   \code{\link{mvsusie_suff_stat}}.
#' 
#' @export
#' 
mvsusie_rss = function (Z, R, N, Bhat, Shat, varY, 
                        prior_variance=0.2, 
                        residual_variance=NULL,
                        ...) {
  is_numeric_prior =
    !(is.matrix(prior_variance) || inherits(prior_variance,"MashInitializer"))
  
  if (sum(c(missing(Z),missing(Bhat) || missing(Shat))) != 1)
    stop("Please provide either Z or (Bhat, Shat), but not both")
  
  # Check input.
  if(anyNA(R)){
    stop("Input R matrix contains NAs.")
  }
  if (missing(Z)){
    J <- ifelse(is.matrix(Bhat), nrow(Bhat), length(Bhat))
  }else{
    J <- ifelse(is.matrix(Z), nrow(Z), length(Z))
  }
  if (nrow(R) != J)
    stop(paste0("The dimension of R (",nrow(R)," x ",ncol(R),") does not ",
                "agree with expected (",J," x ",J,")"))
  
  # Check input N.
  if (!missing(N))
    if (N <= 1)
      stop("N must be greater than 1")
  
  if (missing(Z)) {
    if (length(Shat) == 1){
      if (is.matrix(Bhat)){
        Shat = matrix(Shat, nrow(Bhat), ncol(Bhat))
      }else{
        Shat = rep(Shat, length(Bhat))
      }
    }
    if(is.matrix(Bhat)){
      if (nrow(Bhat) != nrow(Shat))
        stop("The number of rows of Bhat and Shat do not agree")
      if (ncol(Bhat) != ncol(Shat))
        stop("The number of columns of Bhat and Shat do not agree")
    }else{
      if(length != length(Shat)){
        stop("The length of Bhat and Shat do not agree")
      }
    }
    if (anyNA(Bhat) || anyNA(Shat))
      stop("Bhat, Shat cannot have missing values")
    if (any(Shat <= 0))
      stop("Shat cannot have zero or negative elements")
    Z = Bhat/Shat
  }
  if (anyNA(Z)) {
    warning("NA values in z-scores are replaced with 0")
    Z[is.na(Z)] = 0
  }
  
  if (!missing(N)) {
    adj = (N-1)/(Z^2 + N - 2)
    Z   = sqrt(adj) * Z
  }
  
  if (!is.null(dim(Z)) && ncol(Z) > 1 && is_numeric_prior)
    stop("Please specify prior variance for the multivariate z-scores")
  
  is_numeric_matrix(R,"R")
  
  if (missing(N)) {
    warning("Providing the sample size (N), or even a rough estimate of N, ",
            "is highly recommended. Without N, the implicit assumption is ",
            "N is large (Inf) and the effect sizes are small (close to zero).")
    if(!missing(varY)){
      if(!is.null(dim(Z))){
        varY = cov2cor(varY)
      }
    }else{
      if(is.null(residual_variance)){
        if(is.null(dim(Z))){
          varY = 1
        }else{
          varY = diag(ncol(Z))
        }
      }else{
        if(is.null(dim(residual_variance))){
          varY = residual_variance
        }else{
          varY = cov2cor(residual_variance)
        }
      }
    }
    s = mvsusie_suff_stat(XtX = R,XtY = Z,YtY = varY,N = 2,
                          prior_variance = prior_variance, 
                          standardize = FALSE, 
                          residual_variance = residual_variance, ...)
  } else {
    # The sample size (N) is provided, so use PVE-adjusted z-scores.
    if (!missing(Shat) & !missing(varY)) {
      # var_y, shat (and bhat) are provided, so the effects are on the
      # *original scale*.
      if(is.null(dim(Z))){
        XtXdiag = varY * adj/(Shat^2)
        XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
        XtX = (XtX + t(XtX))/2
        XtY = Z * sqrt(adj) * varY / Shat
      }else{
        XtXdiag = colMeans(diag(varY) * adj/t(Shat^2))
        XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
        XtX = (XtX + t(XtX))/2
        XtY = Z * sqrt(adj) * diag(varY) / Shat
      }
    } else {
      # The effects are on the *standardized* X, y scale.
      XtX = (N-1)*R
      XtY = sqrt(N-1)*Z
      if(!missing(varY)){
        if(!is.null(dim(Z))){
          varY = cov2cor(varY)
        }
      }else{
        if(is.null(residual_variance)){
          if(is.null(dim(Z))){
            varY = 1
          }else{
            varY = diag(ncol(Z))
          }
        }else{
          if(is.null(dim(residual_variance))){
            varY = residual_variance
          }else{
            varY = cov2cor(residual_variance)
          }
        }
      }
    }
    s = mvsusie_suff_stat(XtX = XtX, XtY = XtY, YtY = (N-1)*varY, N = N, 
                          prior_variance = prior_variance, 
                          residual_variance = residual_variance, ...)
  }
  
  return(s)
}
