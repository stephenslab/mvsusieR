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
      if (any(dim(X) == 0)) stop('X input dimension is invalid.')
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
      if(private$.Y_has_missing){
        # Length R
        # sum_i Gamma_i Sigma_i^{-1} Gamma_i b^T X[i,] 
        D = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                 t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                     private$missing_pattern[private$Y_missing_pattern_assign[i],]) %*% crossprod(b, private$.X[i,])))
        if (!is.null(private$Y_mean)) intercept = private$Y_mean - private$Ainv %*% D
        else intercept = 0
        if(is.null(dim(coefs))){
          c(intercept, coefs)
        }else{
          mat = as.matrix(rbind(t(intercept), coefs))
          rownames(mat) = NULL
          return(mat)
        }
      }else{
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
    },
    get_sumstats = function(residual_variances, residual_correlation=NULL, sbhat_only=FALSE) {
      if (is.null(residual_correlation)) {
        if (!is.matrix(residual_variances))
          stop("residual variance has to be a matrix if residual correlation is not specified")
        residual_correlation = cov2cor(residual_variances)
        residual_variances = diag(residual_variances)
      }
      # private$d is either vector or matrix
      sbhat = sqrt(do.call(rbind, lapply(1:length(private$d), function(j) residual_variances / private$d[j])))
      sbhat[which(is.nan(sbhat) | is.infinite(sbhat))] = 1E3
      is_common_sbhat = is_mat_common(sbhat)
      if (!sbhat_only) {
        bhat = self$XtY/private$d
        bhat[which(is.nan(bhat))] = 0
        if (is_common_sbhat) SVS = sbhat[1,] * t(residual_correlation * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
        else SVS = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(residual_correlation * sbhat[j,]))
      } else {
        SVS = bhat = NA
      }
      return(list(svs=SVS, sbhat=sbhat, is_common_sbhat = is_common_sbhat, bhat=bhat))
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
    Y_has_missing = function() private$.Y_has_missing
  ),
  private = list(
    .X = NULL,
    X_for_Y_missing = NULL,
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
    Y_non_missing = NULL,
    .Y_has_missing = FALSE,
    .X_has_missing = NULL
  )
)

#' @title Regression data object with missing values in Y
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseDataYMissing <- R6Class("DenseDataYMissing",
  inherit = DenseData,
  public = list(
    initialize = function(X,Y) {
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
      # store missing pattern
      private$missing_pattern = unique(private$Y_non_missing)
      private$Y_missing_pattern_assign = numeric(private$N)
      for(k in 1:nrow(private$missing_pattern)){
        idx = which(apply(private$Y_non_missing, 1, function(x) identical(x, private$missing_pattern[k,])))
        private$Y_missing_pattern_assign[idx] = k
      }
      private$.Y[Y_missing] = 0
      private$.residual = private$.Y
      private$X_for_Y_missing = array(X, dim = c(private$N, private$J, private$R, private$R))
      if(private$R > 1){
        for(i in 1:private$N){
          for(j in 1:private$J){
            private$X_for_Y_missing[i,j,,] = diag(diag(private$X_for_Y_missing[i,j,,]))
          }
        }
      }
    },
    adjust = function(residual_variance, approximation = FALSE){
      private$.residual_variance_inv = list()
      if(approximation){
        ## FIXME: order of eigenvectors? 
        ## Ignore this part for now
        n_per_condition = crossprod(private$Y_non_missing)
        # adjusted residual variance
        residual_variance = residual_variance * cov2cor(n_per_condition)
        private$.residual_variance = residual_variance
        eigenSigma = eigen(residual_variance, symmetric = TRUE)
        dinv = 1/(eigenSigma$values)
        dinv[is.infinite(dinv)] = 0
        for(k in 1:nrow(private$missing_pattern)){
          private$.residual_variance_inv[[k]] = eigenSigma$vectors %*% (dinv[private$missing_pattern[k,]] * 
                                                                          t(eigenSigma$vectors))
        }
      }else{
        private$.residual_variance = residual_variance
        for(k in 1:nrow(private$missing_pattern)){
          if(private$R == 1){
            private$.residual_variance_inv[[k]] = private$missing_pattern[k,]/residual_variance
          }else{
            eigenSigmak = eigen(t(residual_variance * private$missing_pattern[k,]) * private$missing_pattern[k,], symmetric = TRUE)
            dinv = 1/(eigenSigmak$values)
            dinv[is.infinite(dinv)] = 0
            private$.residual_variance_inv[[k]] = eigenSigmak$vectors %*% (dinv * t(eigenSigmak$vectors))
          }
        }
        # For missing Y, d is a J by R by R array
        private$d = array(0, dim=c(private$J, private$R, private$R))
        for(j in 1:private$J){
          # For variant j, sum_i private$X_for_Y_missing[i,j,,] Gamma_i Sigma_i^{-1} Gamma_i private$X_for_Y_missing[i,j,,], R by R matrix
          # when there is no missing, it is sum(x_j^2) * Sigma^{-1}
          private$d[j,,] = Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*% 
                                                (private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                                       private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*% 
                                                private$X_for_Y_missing[i,j,,]))
        }
      }
    },
    standardize = function(center, scale) {
      if(scale==TRUE){
        warning('Scaling varaibles is currently not available with missing data.')
        scale = FALSE
      }
      if(center){
        # sum_i Gamma_i Sigma_i^{-1} Gamma_i R by R matrix
        A = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                 t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                     private$missing_pattern[private$Y_missing_pattern_assign[i],]) ))
        # sum_i Gamma_i Sigma_i^{-1} Gamma_i y_i R by 1 matrix
        B = Reduce('+', lapply(1:private$N, function(i) private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                 t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                     private$missing_pattern[private$Y_missing_pattern_assign[i],]) %*% private$.Y[i,] ))
        # For variant j, sum_i X_{i,j} Gamma_i Sigma_i^{-1} Gamma_i R by R matrix
        C = array(0,dim = c(private$J, private$R, private$R)) # J by R by R array
        for(j in 1:private$J){
          C[j,,] = Reduce('+', lapply(1:private$N, function(i) private$.X[i,j] * (private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                                                                    t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                                                                        private$missing_pattern[private$Y_missing_pattern_assign[i],])) ))
        }
        private$Ainv = invert_via_chol(A)
        # center Y
        private$Y_mean = as.numeric(private$Ainv %*% B)
        private$.Y = t(t(private$.Y) - private$Y_mean)
        private$.Y[!private$Y_non_missing] = 0
        # center X
        for(i in 1:private$N){
          for(j in 1:private$J){
            private$X_for_Y_missing[i,j,,] = private$X_for_Y_missing[i,j,,] - private$Ainv %*% C[j,,]
          }
        }
        # For missing Y, d is a J by R by R array
        for(j in 1:private$J){
          # For variant j, sum_i private$X_for_Y_missing[i,j,,] Gamma_i Sigma_i^{-1} Gamma_i private$X_for_Y_missing[i,j,,], R by R matrix
          # when there is no missing, it is sum(x_j^2) * Sigma^{-1}
          private$d[j,,] = Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*% 
                                                (private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                                   t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                                       private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*% 
                                                private$X_for_Y_missing[i,j,,]))
        }
      }
    },
    compute_Xb = function(b) {
      if(is.vector(b)){
        b = matrix(b, length(b),1)
      }
      Xb = t(sapply(1:private$N, function(i){
        Reduce('+', lapply(1:private$J, function(j) private$X_for_Y_missing[i,j,,] %*% b[j,]))
      }))
      if(nrow(Xb) != private$N) Xb = t(Xb)
      return(Xb)
    },
    # Compute multivariate summary statistics in the presence of missing data
    get_sumstats = function(residual_variances, residual_correlation=NULL, sbhat_only=FALSE) {
      if (!sbhat_only) {
        # private$XtY is a J by R matrix, private$d is a J by R by R array
        bhat = t(sapply(1:private$J, function(j) solve(private$d[j,,], self$XtY[j,])))
        bhat[which(is.nan(bhat))] = 0
        # length J list, each element is an R by R matrix
        SVS = lapply(1:private$J, function(j) solve(private$d[j,,]))
      
        is_common_sbhat = is_list_common(SVS)
        if (is_common_sbhat) {
          SVS = SVS[[1]]
        }
        # FIXME: haven't figured out how to compute it for missing data case ...
        # but this will not have an impact later since we've computed SVS instead
        # that can be used for likelihood and posterior calculations
        sbhat = NA
        return(list(svs=SVS, sbhat=sbhat, is_common_sbhat = is_common_sbhat, bhat=bhat))
      }else{
        svs = bhat = NA
        # FIXME: haven't figured out how to compute it for missing data case ...
        stop("Cannot provide sbhat statistics when there is missing data.")
      }
    }
  ),
  active = list(
    XtY = function() {
      if (is.null(private$.XtY))
        # J by R matrix
        private$.XtY = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*% 
                                                                              (private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                                                                 t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                                                                     private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*% 
                                                                              private$.Y[i,]))))
      if(nrow(private$.XtY) != private$J) private$.XtY = t(private$.XtY)
      return(private$.XtY)
    },
    XtX = function() {
      if (is.null(private$.XtX))
        # FIXME: not sure how to compute XtX with missing data
        private$.XtX = cor(private$.X)
      return(private$.XtX)
    },
    XtR = function() {
      res = t(sapply(1:private$J, function(j) Reduce('+', lapply(1:private$N, function(i) private$X_for_Y_missing[i,j,,] %*% 
                                                                   (private$missing_pattern[private$Y_missing_pattern_assign[i],] * 
                                                                      t(private$.residual_variance_inv[[private$Y_missing_pattern_assign[i]]] * 
                                                                          private$missing_pattern[private$Y_missing_pattern_assign[i],])) %*% 
                                                                   private$.residual[i,]))))
      if(nrow(res) != private$J) res = t(res)
      return(res)
    }
  ),
  private = list(
    missing_pattern = NULL,
    Y_missing_pattern_assign = NULL,
    Ainv = NULL,
    .residual_variance = NULL,
    .residual_variance_inv = NULL
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
