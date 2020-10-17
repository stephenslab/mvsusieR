#' @title The regular regression data object
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseData <- R6Class("DenseData",
  public = list(
    initialize = function(X,Y) {
      # Check input X.
      if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix"))
        stop("Input X must be a double-precision matrix, or a sparse matrix.")
      if (any(is.na(X))) {
        stop("Input X must not contain missing values.")
      }
      if (any(dim(X) == 0)) stop('Input X dimension is invalid.')
      if (length(which(apply(X, 2, is_zero_variance)))) stop('Input X must not have constant columns (some columns have standard deviation zero)')
      private$.X = X
      private$.X_has_missing = any(is.na(private$.X))
      # FIXME: might want to allow for missing in X later?
      # see stephenslab/mmbr/#5
      if (private$.X_has_missing)
        stop("Missing data in input matrix X is not allowed at this point.")
      if (is.null(dim(Y))) private$.Y = matrix(Y,length(Y),1)
      else private$.Y = Y
      private$.residual = private$.Y
      private$R = ncol(private$.Y)
      private$N = nrow(private$.Y)
      private$J = ncol(X)
      # quantities involved in center and scaling
      private$cm = rep(0, length = private$J)
      private$csd = rep(1, length = private$J)
      private$d = colSums(private$.X ^ 2)
      private$d[private$d == 0] = 1E-6
    },
    set_residual_variance = function(residual_variance=NULL, numeric = FALSE,
                                     precompute_covariances = TRUE,
                                     quantities = c('residual_variance','effect_variance')){
      if('residual_variance' %in% quantities){
        if (is.null(residual_variance)) {
          if (private$R > 1) {
            if (!private$.Y_has_missing) residual_variance = cov(private$.Y)
            else stop("Unspecified residual_variance is not allowed in the presence of missing data in Y")
          }
          else residual_variance = var(private$.Y, na.rm=T)
        }
        if(numeric){
          residual_variance = as.numeric(residual_variance)
        }
        if (is.matrix(residual_variance)) {
          if(nrow(residual_variance) != private$R){
            stop(paste0("The residual variance is not a ", private$R, ' by ', private$R, ' matrix.'))
          }
          if (any(is.na(diag(residual_variance))))
            stop("Diagonal of residual_variance cannot be NA")
          residual_variance[which(is.na(residual_variance))] = 0
          mashr:::check_positive_definite(residual_variance)
          private$.residual_correlation = cov2cor(as.matrix(residual_variance))
        }else {
          if (is.na(residual_variance) || is.infinite(residual_variance))
            stop("Invalid residual_variance")
          private$.residual_correlation = 1
        }
        private$.residual_variance = residual_variance
        tryCatch({
          private$.residual_variance_inv = invert_via_chol(residual_variance)
        }, error = function(e) {
          stop(paste0('Cannot compute inverse for residual_variance:\n', e))
        })
      }
      if('effect_variance' %in% quantities){
        if(precompute_covariances){
          private$.svs = lapply(1:private$J, function(j){
            res = private$.residual_variance /private$d[j]
            res[which(is.nan(res) | is.infinite(res))] = 1E6
            return(res)
          })
          private$.svs_inv = lapply(1:private$J, function(j) private$.residual_variance_inv * private$d[j])
          private$.is_common_sbhat = is_list_common(private$.svs)
        }else{
          private$.sbhat = sqrt(do.call(rbind, lapply(1:private$J, function(j) diag(as.matrix(private$.residual_variance)) / private$d[j])))
          private$.sbhat[which(is.nan(private$.sbhat) | is.infinite(private$.sbhat))] = 1E3
          private$.is_common_sbhat = is_mat_common(private$.sbhat)
        }
      }
    },
    get_coef = function(use_residual = FALSE){
      # XtY J by R matrix
      if (use_residual) XtY = self$XtR
      else XtY = self$XtY
      bhat = XtY / private$d # private$d length J vector
      bhat[which(is.nan(bhat))] = 0
      return(bhat)
    },
    standardize = function(center, scale) {
      # Credit: This is heavily based on code from
      # https://www.r-bloggers.com/a-faster-scale-function/
      # The only change from that code is its treatment of columns with 0 variance.
      # This "safe" version scales those columns by 1 instead of 0.
      if(center){
        # center X
        private$cm = colMeans(private$.X, na.rm=TRUE)
        # center Y
        if (private$R == 1) private$Y_mean = mean(private$.Y, na.rm=TRUE)
        else private$Y_mean = colMeans(private$.Y, na.rm=TRUE)
        private$.Y = t(t(private$.Y) - private$Y_mean)
      }
      if (scale) {
        private$csd = colSds(private$.X, center = private$cm)
        private$csd[private$csd==0] = 1
      }
      private$.X = t( (t(private$.X) - private$cm) / private$csd )
      private$d = colSums(private$.X ^ 2)
      private$d[private$d == 0] = 1E-6
    },
    compute_Xb = function(b) {
      # J by R
      # tcrossprod(A,B) performs A%*%t(B) but faster
      tcrossprod(private$.X,t(b))
    },
    compute_MXt = function(M) {
      # tcrossprod(A,B) performs A%*%t(B) but faster
      tcrossprod(M, private$.X)
    },
    remove_from_residual = function(value) {
      private$.residual = private$.residual - value
    },
    add_to_residual = function(value) {
      private$.residual = private$.residual + value
    },
    compute_residual = function(fitted) {
      private$.residual = private$.Y - fitted
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
      if (is.null(private$.XtY)) private$.XtY = crossprod(private$.X, private$.Y)
      return(private$.XtY)
    },
    XtX = function() {
      if (is.null(private$.XtX)) private$.XtX = crossprod(private$.X)
      return(private$.XtX)
    },
    XtR = function() {
      crossprod(private$.X, private$.residual)
    },
    residual = function() private$.residual,
    n_sample = function() private$N,
    n_condition = function() private$R,
    n_effect = function() private$J,
    X_has_missing = function() private$.X_has_missing,
    Y_has_missing = function() private$.Y_has_missing,
    residual_variance = function() private$.residual_variance,
    residual_variance_inv = function() private$.residual_variance_inv,
    residual_correlation = function() private$.residual_correlation,
    sbhat = function() private$.sbhat,
    svs = function() private$.svs,
    svs_inv = function() private$.svs_inv,
    is_common_cov = function() private$.is_common_sbhat
  ),
  private = list(
    .X = NULL,
    .Y = NULL,
    .XtX = NULL,
    .XtY = NULL,
    d = NULL,
    N = NULL,
    J = NULL,
    R = NULL,
    .residual = NULL,
    csd = NULL,
    cm = NULL,
    Y_mean = NULL,
    .Y_has_missing = FALSE,
    .X_has_missing = NULL,
    .residual_variance = NULL,
    .residual_variance_inv = NULL,
    .residual_correlation = NULL,
    .sbhat = matrix(0,0,0),
    .svs = 0,
    .svs_inv = 0,
    .is_common_sbhat = FALSE
  )
)

#' @title Regression data object with missing values in Y
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseDataYMissing <- R6Class("DenseDataYMissing",
  inherit = DenseData,
  public = list(
    initialize = function(X,Y,approximate=FALSE) {
      # initialize with super class but postpone center and scaling to later
      super$initialize(X,Y)
      Y_missing = is.na(private$.Y)
      private$.Y_has_missing = any(Y_missing)
      if(!private$.Y_has_missing) {
        warning("Y does not have any missing values in it. You should consider using DenseData class instead. Here we force set attribute Y_has_missing = TRUE")
        # To force use this class when there is no missing data in Y
        private$.Y_has_missing = TRUE
      }
      private$Y_non_missing = !Y_missing
      private$approximate = approximate
      # store missing pattern
      private$missing_pattern = unique(private$Y_non_missing)
      Y_missing_pattern_assign = numeric(private$N)
      for(k in 1:nrow(private$missing_pattern)){
        idx = which(apply(private$Y_non_missing, 1, function(x) identical(x, private$missing_pattern[k,])))
        Y_missing_pattern_assign[idx] = k
      }
      private$Y_missing_pattern_assign = Y_missing_pattern_assign
      private$.Y[Y_missing] = 0
      private$.residual = private$.Y
      if(private$approximate){
        X_for_Y_missing = array(private$.X, dim = c(private$N, private$J, private$R))
        for(r in 1:private$R) {
          X_for_Y_missing[Y_missing[,r],,r] = NA
        }
      }else{
        X_for_Y_missing = outer(private$.X,  diag(private$R))
      }
      private$X_for_Y_missing = X_for_Y_missing
    },
    set_residual_variance = function(residual_variance=NULL, numeric = FALSE,
                                     quantities = c('residual_variance','effect_variance')){
      if('residual_variance' %in% quantities){
        if (is.null(residual_variance)) {
          if (private$R > 1) {
            stop("Unspecified residual_variance is not allowed in the presence of missing data in Y")
          }
          else residual_variance = var(private$.Y[private$Y_non_missing], na.rm=T)
        }
        if (is.matrix(residual_variance)) {
          if (any(is.na(diag(residual_variance))))
            stop("Diagonal of residual_variance cannot be NA")
          residual_variance[which(is.na(residual_variance))] = 0
          mashr:::check_positive_definite(residual_variance)
        }else {
          if (is.na(residual_variance) || is.infinite(residual_variance))
            stop("Invalid residual_variance")
        }
        if(numeric){
          residual_variance = as.numeric(residual_variance)
        }
        residual_variance_inv = list()
        private$.residual_variance = residual_variance
        for(k in 1:nrow(private$missing_pattern)){
          if(private$R == 1){
            residual_variance_inv[[k]] = private$missing_pattern[k,]/residual_variance
          }else{
            eigenSigmak = eigen(t(residual_variance * private$missing_pattern[k,]) * private$missing_pattern[k,], symmetric = TRUE)
            eigenSigmak$values[eigenSigmak$values < 0] = 0
            dinv = 1/(eigenSigmak$values)
            dinv[is.infinite(dinv)] = 0
            residual_variance_inv[[k]] = eigenSigmak$vectors %*% (dinv * t(eigenSigmak$vectors))
          }
        }
        private$.residual_variance_inv = residual_variance_inv
        if(is.matrix(residual_variance)){
          private$.residual_correlation = cov2cor(residual_variance)
        }else{
          private$.residual_correlation = 1
        }
      }
      if('effect_variance' %in% quantities){
        svs_inv = list()
        svs = list()
        for(j in 1:private$J){
          # For variant j, sum_i private$X_for_Y_missing[i,j,,] Gamma_i Sigma_i^{-1} Gamma_i private$X_for_Y_missing[i,j,,], R by R matrix
          # when there is no missing, it is sum(x_j^2) * Sigma^{-1}
          if(private$approximate){
            svs_inv[[j]] = Reduce('+', lapply(1:private$N, function(i)
                                                         t(t((private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                            t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                private$missing_pattern[private$Y_missing_pattern_assign[i],])) *
                                                         private$X_for_Y_missing[i,j,])* private$X_for_Y_missing[i,j,])))
            svs[[j]] = tryCatch(invert_via_chol(svs_inv[[j]]), error = function(e){
              invert_via_chol(svs_inv[[j]] + 1e-8 * diag(private$R))} )
          }else{
            svs_inv[[j]] = Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*%
                                                         (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                            t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*%
                                                         private$X_for_Y_missing[i,j,,]))
            svs[[j]] = tryCatch(invert_via_chol(svs_inv[[j]]), error = function(e){
              invert_via_chol(svs_inv[[j]] + 1e-8 * diag(private$R))} )
          }
        }
        private$.svs_inv = svs_inv
        private$.svs = svs
        private$.is_common_sbhat = is_list_common(private$.svs)
      }
    },
    get_coef = function(use_residual = FALSE){
      if (use_residual) XtY = self$XtR
      else XtY = self$XtY
      bhat = t(sapply(1:private$J, function(j) private$.svs[[j]] %*% XtY[j,]))
      bhat[which(is.nan(bhat))] = 0
      if(private$R == 1){
        bhat = t(bhat)
      }
      return(bhat)
    },
    standardize = function(center, scale) {
      if(private$approximate){
        if(center){
          # center X
          private$cm = colMeans(private$X_for_Y_missing, na.rm=T) # J by R
          # center Y
          if (private$R == 1) private$Y_mean = mean(private$.Y[private$Y_non_missing])
          else private$Y_mean = sapply(1:private$R, function(r) mean(private$.Y[private$Y_non_missing[,r],r]))
          private$.Y = t(t(private$.Y) - private$Y_mean)
          private$.Y[!private$Y_non_missing] = 0
        }else{
          private$cm = matrix(0, private$J, private$R) # J by R
        }
        private$csd = matrix(1, private$J, private$R)
        # scale X when Y has missing, and compute colSums(X^2)
        X_for_Y_missing = private$X_for_Y_missing
        for(r in 1:private$R){
          if (scale) {
            private$csd[,r] = colSds(private$X_for_Y_missing[,,r], center = private$cm[,r], na.rm = TRUE)
            private$csd[private$csd[,r]==0, r] = 1
          }
          X_for_Y_missing[,,r] = t( (t(X_for_Y_missing[,,r]) - private$cm[,r]) / private$csd[,r] )
          X_for_Y_missing[,,r][is.na(X_for_Y_missing[,,r])] = 0
        }
        private$X_for_Y_missing = X_for_Y_missing
      }else{
        if(scale==TRUE){
          private$csd = colSds(private$.X, center = 0)
          private$csd[private$csd==0] = 1
          private$.X = t(t(private$.X) / private$csd )
        }
        if(center){
          # sum_i Gamma_i Sigma_i^{-1} Gamma_i R by R matrix
          A = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                       private$missing_pattern[private$Y_missing_pattern_assign[i],]) ))
          private$Ainv = invert_via_chol(A)
          
          # sum_i Gamma_i Sigma_i^{-1} Gamma_i y_i R by 1 matrix
          B = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                       private$missing_pattern[private$Y_missing_pattern_assign[i],]) %*% private$.Y[i,] ))
          
          # center Y
          private$Y_mean = as.numeric(private$Ainv %*% B)
          private$.Y = t(t(private$.Y) - private$Y_mean)
          private$.Y[!private$Y_non_missing] = 0
          
          # center X
          X_for_Y_missing = private$X_for_Y_missing
          for(j in 1:private$J){
            # For variant j, Ainv sum_i X_{i,j} Gamma_i Sigma_i^{-1} Gamma_i R by R matrix
            tmp = private$Ainv %*% Reduce('+', lapply(1:private$N, function(i) private$.X[i,j] * (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                                                                    t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                                                        private$missing_pattern[private$Y_missing_pattern_assign[i],])) ))
            X_for_Y_missing[,j,,] = sweep(X_for_Y_missing[,j,,,drop=F], 3:4, tmp)
          }
          private$X_for_Y_missing = X_for_Y_missing
        }
      }
    },
    compute_Xb = function(b) {
      if(is.vector(b)){
        b = matrix(b, length(b),1)
      }
      if(private$approximate){
        Xb = sapply(1:private$R, function(r) private$X_for_Y_missing[,,r] %*% b[,r])
      }else{
        Xb = t(sapply(1:private$N, function(i){
          Reduce('+', lapply(1:private$J, function(j) private$X_for_Y_missing[i,j,,] %*% b[j,]))
        }))
      }
      if(nrow(Xb) != private$N) Xb = t(Xb)
      return(Xb)
    },
    rescale_coef = function(b){
      if(private$R == 1){
        private$csd = as.vector(private$csd)
        private$cm = as.vector(private$cm)
      }
      coefs = b/private$csd
      if(private$approximate){
        if (is.null(dim(coefs))) {
          if (!is.null(private$Y_mean)){
            intercept = private$Y_mean - sum(private$cm * coefs)}
          else intercept = 0
          c(intercept, coefs)
        } else {
          if (!is.null(private$Y_mean)) intercept = private$Y_mean - colSums(private$cm * coefs)
          else intercept = 0
          mat = as.matrix(rbind(intercept, coefs))
          rownames(mat) = NULL
          return(mat)
        }
      }else{
        # Length R
        # sum_i Gamma_i Sigma_i^{-1} Gamma_i b^T X[i,]
        D = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                 t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                     private$missing_pattern[private$Y_missing_pattern_assign[i],]) %*% crossprod(coefs, private$.X[i,])))
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - private$Ainv %*% D
        else intercept = 0
        if(is.null(dim(coefs))){
          c(intercept, coefs)
        }else{
          mat = as.matrix(rbind(t(intercept), coefs))
          rownames(mat) = NULL
          return(mat)
        }
      }
    }
  ),
  active = list(
    XtY = function() {
      if (is.null(private$.XtY))
        # J by R matrix
        if(private$approximate){
          private$.XtY = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i)
                                                                                (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                                       private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*%
                                                                                private$.Y[i,] * private$X_for_Y_missing[i,j,]) )))
        }else{
          private$.XtY = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*%
                                                                                (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                                       private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*%
                                                                                private$.Y[i,]))))
        }
      if(private$R == 1) private$.XtY = t(private$.XtY)
      return(private$.XtY)
    },
    XtX = function() {
      if (is.null(private$.XtX))
        # FIXME: not sure how to compute XtX with missing data
        private$.XtX = cor(private$.X)
      return(private$.XtX)
    },
    XtR = function() {
      if(private$approximate){
        res = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i)
                                                                     (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                                        t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                            private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*%
                                                                     private$.residual[i,] * private$X_for_Y_missing[i,j,]) )))
      }else{
        res = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*%
                                                                     (private$missing_pattern[private$Y_missing_pattern_assign[i],] *
                                                                        t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] *
                                                                            private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*%
                                                                     private$.residual[i,]))))
      }
      if(nrow(res) != private$J) res = t(res)
      return(res)
    }
  ),
  private = list(
    X_for_Y_missing = NULL,
    Y_non_missing = NULL,
    missing_pattern = NULL,
    Y_missing_pattern_assign = NULL,
    Ainv = NULL,
    approximate = FALSE
  )
)

#' @title Summary statistics object
# Z and R
#' @importFrom R6 R6Class
#' @keywords internal
RSSData <- R6Class("RSSData",
  inherit = DenseData,
  public = list(
    initialize = function(Z, R, tol) {
      if(any(is.infinite(Z))){
        stop('Z scores contain infinite value.')
        }
      # Check NA in R
      if (any(is.na(R))) {
        stop('R matrix contains missing values.')
      }
      # Check input R.
      if (!susieR:::is_symmetric_matrix(R)) {
        stop('R is not a symmetric matrix.')
      }
      if (!(is.double(R) &
            is.matrix(R)) & !inherits(R, "CsparseMatrix"))
        stop("Input R must be a double-precision matrix, or a sparse matrix.")

      if (is.null(dim(Z)))
        Z = matrix(Z, length(Z), 1)

      if (nrow(R) != nrow(Z)) {
        stop(paste0('The dimension of correlation matrix (',nrow(R),' by ',ncol(R),
            ') does not agree with expected (',nrow(Z),' by ',nrow(Z),')'))
      }

      # replace NA in z with 0
      if (any(is.na(Z))) {
        warning('NA values in Z-scores are replaced with 0.')
        Z[is.na(Z)] = 0
      }
      private$.XtX = R
      private$J = nrow(R)
      if (is.null(dim(Z))) private$R = 1
      else private$R = ncol(Z)
      private$check_semi_pd(tol)
      private$.X = t(private$eigenvectors[, private$eigenvalues !=0]) * private$eigenvalues[private$eigenvalues != 0] ^ (0.5)
      private$.Y = (t(private$eigenvectors[, private$eigenvalues != 0]) * private$eigenvalues[private$eigenvalues != 0] ^ (-0.5)) %*% Z
      private$.Y_has_missing = private$.X_has_missing = FALSE
      private$.XtY = private$UUt %*% Z # = Z when Z is in eigen space of R
      private$.residual = private$.Y
    }
  ),
  active = list(
    XtX = function() private$.XtX,
    XtY = function() private$.XtY,
    # n_sample doesn't mean sample size here, it means the number of non zero eigenvalues
    n_sample = function() sum(private$eigenvalues > 0)
  ),
  private = list(
    UUt = NULL,
    eigenvectors = NULL,
    eigenvalues = NULL,
    check_semi_pd = function(tol) {
      eigenR = eigen(private$.XtX, symmetric = TRUE)
      eigenR$values[abs(eigenR$values) < tol] = 0
      if (any(eigenR$values < 0)) {
        eigenR$values[eigenR$values < 0] = 0
        warning('Negative eigenvalues are set to 0.')
      }
      private$.XtX = eigenR$vectors %*% (t(eigenR$vectors) * eigenR$values)
      private$csd = rep(1, length = private$J)
      private$d = diag(private$.XtX)
      private$eigenvectors = eigenR$vectors
      private$eigenvalues = eigenR$values
      private$UUt = tcrossprod(private$eigenvectors[, which(private$eigenvalues > 0)])
    }
  )
)
