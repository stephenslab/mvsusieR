# The regular regression data object.
#
# See discussion in
# https://github.com/r-lib/R6/issues/213
# https://github.com/r-lib/R6/issues/201
# for why we use <<- instead of <- for assignnment.
#
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' 
DenseData <- R6Class("DenseData",
  portable = FALSE,
  public = list(
    initialize = function (X, Y) {
      is_numeric_matrix(X,"X")
      if (length(which(apply(X,2,is_zero_variance))))
        stop("Input X must not have constant columns (some columns have ",
             "standard deviation zero)")
      .X             <<- X
      .X_has_missing <<- anyNA(.X)
      
      # FIXME: might want to allow for missing in X later?
      # see stephenslab/mvsusieR/#5
      if (.X_has_missing)
        stop("Missing data in input matrix X is not allowed at this point")
      if (is.null(dim(Y)))
        .Y <<- matrix(Y,length(Y),1)
      else
        .Y <<- Y
      .Y_missing     <<- is.na(.Y)
      .Y_has_missing <<- any(.Y_missing)
      .residual      <<- .Y
      .R             <<- ncol(.Y)
      .N             <<- nrow(.Y)
      .J             <<- ncol(.X)
      
      # Quantities involved in center and scaling.
      .cm  <<- rep(0, length.out = .J)
      .csd <<- rep(1, length.out = .J)
      .d   <<- colSums(.X ^ 2)
      .d[.d == 0] <<- 1e-6
      return(invisible(self))
    },

    # This returns the R6 object invisibly.
    set_residual_variance = function (residual_variance = NULL,
                                      numeric = FALSE,
                                      precompute_covariances = TRUE,
                                      quantities = c("residual_variance",
                                                     "effect_variance")) {
      if ("residual_variance" %in% quantities) {
        if (is.null(residual_variance)) {
          if (.R > 1) {
            if (!.Y_has_missing)
              residual_variance = cov(.Y)
            else
              stop("Unspecified residual_variance is not allowed in the ",
                   "presence of missing data in Y")
          }
          else
            residual_variance = var(.Y, na.rm = TRUE)
        }
        if (numeric)
          residual_variance = as.numeric(residual_variance)
        if (is.matrix(residual_variance)) {
          if (nrow(residual_variance) != .R)
            stop(paste0("The residual variance is not a ",.R," by ",.R,
                        " matrix"))
          if (anyNA(diag(residual_variance)))
            stop("Diagonal of residual_variance cannot be NA")
          residual_variance[which(is.na(residual_variance))] = 0
          mashr:::check_positive_definite(residual_variance)
          .residual_correlation <<- cov2cor(as.matrix(residual_variance))
        } else {
          if (is.na(residual_variance) || is.infinite(residual_variance))
            stop("Invalid residual_variance")
          .residual_correlation <<- 1
        }
        .residual_variance <<- residual_variance
        tryCatch({
          .residual_variance_inv <<- invert_via_chol(residual_variance)$inv
        },error = function (e) {
          stop(paste0("Cannot compute inverse for residual_variance:\n",e))
        })
      }
      
      if ("effect_variance" %in% quantities) {
        if (precompute_covariances) {
          .svs <<- lapply(1:.J,
                          function (j) {
                            res = .residual_variance /.d[j]
                            res[which(is.nan(res) | is.infinite(res))] = 1e6
                            return(res)
                          })
          .svs_inv <<- lapply(1:.J,function (j) .residual_variance_inv * .d[j])
          .is_common_sbhat <<- is_list_common(.svs)
        } else {
          .sbhat <<- sqrt(do.call(rbind,lapply(1:.J,
                       function(j) diag(as.matrix(.residual_variance))/.d[j])))
          .sbhat[which(is.nan(.sbhat) | is.infinite(.sbhat))] <<- 1e3
          .is_common_sbhat <<- is_mat_common(.sbhat)
        }
      }

      return(invisible(self))
    },
      
    # Return XtY, J x R matrix
    get_coef = function (use_residual = FALSE) {
      if (use_residual)
        XtY = self$XtR
      else
        XtY = self$XtY
      bhat = XtY/.d  # .d is a length-J vector
      bhat[which(is.nan(bhat))] = 0
      return(bhat)
    },
      
    # This is heavily based on code from
    # www.r-bloggers.com/a-faster-scale-function. The only change from
    # that code is its treatment of columns with 0 variance. This
    # "safe" version scales those columns by 1 instead of 0.
    #  
    # This method returns the R6 object invisibly.
    standardize = function (center, scale) {
      if (center) {
          
        # Center X
        .cm <<- colMeans(.X,na.rm = TRUE)
        
        # Center Y.
        if (.R == 1)
          .Y_mean <<- mean(.Y,na.rm = TRUE)
        else
          .Y_mean <<- colMeans(.Y,na.rm = TRUE)
        .Y <<- t(t(.Y) - .Y_mean)
      }
      
      if (scale) {
        .csd <<- colSds(.X, center = .cm)
        .csd[.csd==0] <<- 1
      }
      
      .X <<- t((t(.X) - .cm)/.csd)
      .d <<- colSums(.X ^ 2)
      .d[.d == 0] <<- 1e-6
      
      return(invisible(self))
    },
      
    # J x R
    # Performs A %*% t(B) but faster.
    compute_Xb = function (b)
      tcrossprod(.X,t(b)),
        
    # Performs A %*% t(B) but faster.
    compute_MXt = function (M)
      tcrossprod(M,.X),
      
    remove_from_residual = function (value) {
      .residual <<- .residual - value
      return(.residual)
    },
      
    add_to_residual = function (value) {
      .residual <<- .residual + value
      return(.residual)
    },
      
    compute_residual = function (fitted) {
      .residual <<- .Y - fitted
      return(.residual)
    },
      
    rescale_coef = function(b) {
      coefs = b/.csd
      if (is.null(dim(coefs))) {
        if (!is.null(.Y_mean))
          intercept = .Y_mean - sum(.cm * coefs)
        else
          intercept = as.numeric(NA)
        return(c(intercept,coefs))
      } else {
        if (!is.null(.Y_mean))
          intercept = .Y_mean - colSums(.cm * coefs)
        else
          intercept = as.numeric(NA)
        mat = as.matrix(rbind(intercept, coefs))
        rownames(mat) = NULL
        return(mat)
      }
    }
  ),
                     
  active = list(
    X = function()
      .X,
      
    Y = function()
      .Y,
      
    X2_sum = function()
      .d,
      
    XtY = function() {
      if (is.null(.XtY))
        .XtY <<- crossprod(.X,.Y)
      return(.XtY)
    },
      
    XtX = function() {
      if (is.null(.XtX))
        .XtX <<- crossprod(.X)
      return(.XtX)
    },
      
    XtR = function()
      crossprod(.X,.residual),
    
    residual              = function() .residual,
    n_sample              = function() .N,
    n_condition           = function() .R,
    n_effect              = function() .J,
    X_has_missing         = function() .X_has_missing,
    Y_has_missing         = function() .Y_has_missing,
    residual_variance     = function() .residual_variance,
    residual_variance_inv = function() .residual_variance_inv,
    residual_correlation  = function() .residual_correlation,
    sbhat                 = function() .sbhat,
    svs                   = function() .svs,
    svs_inv               = function() .svs_inv,
    is_common_cov         = function() .is_common_sbhat
  ),
                     
  private = list(
    .X        = NULL,
    .Y        = NULL,
    .XtX      = NULL,
    .XtY      = NULL,
    .d        = NULL,
    .N        = NULL,
    .J        = NULL,
    .R        = NULL,
    .residual = NULL,
    .csd      = NULL,
    .cm       = NULL,
    .Y_mean   = NULL,
    .Y_has_missing         = FALSE,
    .Y_missing             = NULL,
    .X_has_missing         = NULL,
    .residual_variance     = NULL,
    .residual_variance_inv = NULL,
    .residual_correlation  = NULL,
    .sbhat   = matrix(0,0,0),
    .svs     = 0,
    .svs_inv = 0,
    .is_common_sbhat = FALSE
  )
)
