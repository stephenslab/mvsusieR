# Sufficient statistics object for XtX, XtY, YtY, N.
#
#' @importFrom R6 R6Class
SSData <- R6Class("SSData",
  inherit = DenseData,
  portable = FALSE,
  public = list(
    initialize = function(XtX, XtY, YtY, N, X_colmeans, Y_colmeans) {
      if (is.null(dim(XtY))) {
        .XtY <<- matrix(XtY, length(XtY), 1)
      } else {
        .XtY <<- XtY
      }
      if (ncol(XtX) != nrow(.XtY)) {
        stop(paste0(
          "The dimension of XtX (", nrow(XtX), " by ", ncol(XtX), ") ",
          ") does not agree with expected (", nrow(.XtY), " by ",
          nrow(.XtY), ")"
        ))
      }
      if (!susieR:::is_symmetric_matrix(XtX)) {
        stop("XtX is not a symmetric matrix")
      }
      if (any(is.infinite(.XtY))) {
        stop("XtY contains infinite values")
      }
      is_numeric_matrix(XtX, "XtX")

      .XtX <<- XtX
      .YtY <<- YtY
      .Y_has_missing <<- FALSE
      .Xtresidual <<- .XtY
      .R <<- ncol(.XtY)
      .N <<- N
      .J <<- ncol(.XtX)

      # Quantities used in centering and scaling.
      .cm <<- rep(0, length.out = .J)
      .csd <<- rep(1, length.out = .J)
      .d <<- diag(.XtX)
      .d[.d == 0] <<- 1e-6

      if (all(!is.null(X_colmeans), !is.null(Y_colmeans))) {
        if (length(X_colmeans) == 1 && X_colmeans == 0) {
          X_colmeans <- numeric(.J)
        }
        if (length(Y_colmeans) == 1 && Y_colmeans == 0) {
          y_colmeans <- numeric(.R)
        }
        if (length(X_colmeans) != .J) {
          stop(
            "The length of X_colmeans does not agree with number of ",
            "variables"
          )
        }
        if (length(Y_colmeans) != .R) {
          stop(
            "The length of Y_colmeans does not agree with number of ",
            "conditions"
          )
        }
      }
      .cm <<- X_colmeans
      .Y_mean <<- Y_colmeans

      return(invisible(self))
    },

    # This returns the R6 object invisibly.
    set_residual_variance = function(residual_variance = NULL,
                                     numeric = FALSE,
                                     precompute_covariances = TRUE,
                                     quantities = c(
                                       "residual_variance",
                                       "effect_variance"
                                     )) {
      if ("residual_variance" %in% quantities) {
        if (is.null(residual_variance)) {
          if (.R > 1) {
            residual_variance <- cov2cor(.YtY)
          } else {
            residual_variance <- .YtY / (.N - 1)
          }
        }
        if (numeric) {
          residual_variance <- as.numeric(residual_variance)
        }
        if (is.matrix(residual_variance)) {
          if (nrow(residual_variance) != .R) {
            stop(paste0(
              "The residual variance is not a ", .R, " by ", .R,
              " matrix"
            ))
          }
          if (anyNA(diag(residual_variance))) {
            stop("Diagonal of residual_variance cannot be NA")
          }
          residual_variance[which(is.na(residual_variance))] <- 0
          mashr:::check_positive_definite(residual_variance)
          .residual_correlation <<- cov2cor(as.matrix(residual_variance))
        } else {
          if (is.na(residual_variance) || is.infinite(residual_variance)) {
            stop("Invalid residual_variance")
          }
          .residual_correlation <<- 1
        }
        .residual_variance <<- residual_variance
        tryCatch(
          {
            .residual_variance_inv <<-
              invert_via_chol(residual_variance)$inv
          },
          error = function(e) {
            stop(paste0(
              "Cannot compute inverse for ",
              "residual_variance:\n", e
            ))
          }
        )
      }

      if ("effect_variance" %in% quantities) {
        if (precompute_covariances) {
          .svs <<- lapply(1:.J, function(j) {
            res <- .residual_variance / .d[j]
            res[which(is.nan(res) | is.infinite(res))] <- 1e6
            return(res)
          })
          .svs_inv <<- lapply(1:.J, function(j) .residual_variance_inv * .d[j])
          .is_common_sbhat <<- is_list_common(.svs)
        } else {
          .sbhat <<-
            sqrt(do.call(
              rbind,
              lapply(
                1:.J,
                function(j) diag(as.matrix(.residual_variance)) / .d[j]
              )
            ))
          .sbhat[which(is.nan(.sbhat) | is.infinite(.sbhat))] <<- 1e3
          .is_common_sbhat <<- is_mat_common(.sbhat)
        }
      }

      return(invisible(self))
    },
    standardize = function(scale) {
      if (scale) {
        dXtX <- diag(.XtX)
        .csd <<- sqrt(dXtX / (.N - 1))
        .csd[.csd == 0] <<- 1
        .XtX <<- (1 / .csd) * t((1 / .csd) * XtX)
        .XtY <<- (1 / .csd) * XtY
        .d <<- diag(.XtX)
        .d[.d == 0] <<- 1e-6
      }
      return(invisible(self))
    },

    # J x R
    # tcrossprod(A,B) performs A%*%t(B) but faster
    compute_Xb = function(b) {
      tcrossprod(.XtX, t(b))
    },

    # tcrossprod(A,B) performs A %*% t(B) but faster
    compute_MXt = function(M) {
      tcrossprod(M, .XtX)
    },
    remove_from_residual = function(value) {
      .Xtresidual <<- .Xtresidual - value
      return(.Xtresidual)
    },
    add_to_residual = function(value) {
      .Xtresidual <<- .Xtresidual + value
      return(.Xtresidual)
    },
    compute_residual = function(fitted) {
      .Xtresidual <<- .XtY - fitted
      return(.Xtresidual)
    },
    rescale_coef = function(b) {
      coefs <- b / .csd
      if (is.null(dim(coefs))) {
        if (any(is.null(.cm), is.null(.Y_mean))) {
          intercept <- as.numeric(NA)
        } else {
          intercept <- .Y_mean - sum(.cm * coefs)
        }
        return(c(intercept, coefs))
      } else {
        if (any(is.null(.cm), is.null(.Y_mean))) {
          intercept <- as.numeric(NA)
        } else {
          intercept <- .Y_mean - colSums(.cm * coefs)
        }
        mat <- as.matrix(rbind(intercept, coefs))
        rownames(mat) <- NULL
        return(mat)
      }
    }
  ),
  active = list(
    YtY      = function() .YtY,
    XtX      = function() .XtX,
    XtY      = function() .XtY,
    XtR      = function() .Xtresidual,
    residual = function() .Xtresidual
  )
)
