#' @title Regression data object with missing values in Y
#' @importFrom R6 R6Class
#' @importFrom matrixStats colSds
#' @keywords internal
DenseDataYMissing <- R6Class("DenseDataYMissing",
  inherit = DenseData,
  portable = FALSE,
  public = list(
    initialize = function (X, Y, approximate = FALSE) {
        
      # Initialize with super class but postpone center and scaling to
      # later.
      super$initialize(X,Y)
      if(!.Y_has_missing) {
        warning("Y does not have any missing values in it. You should consider using DenseData class instead. Here we force set attribute Y_has_missing = TRUE")
        # To force use this class when there is no missing data in Y
        .Y_has_missing <<- TRUE
      }
      .Y_non_missing <<- !.Y_missing
      .approximate <<- approximate
      .missing_pattern <<- unique(.Y_non_missing)  # Store missing pattern.
      .Y_missing_pattern_assign <<- numeric(.N)
      for(k in 1:nrow(.missing_pattern)){
        idx = which(apply(.Y_non_missing, 1, function(x) identical(x, .missing_pattern[k,])))
        .Y_missing_pattern_assign[idx] <<- k
      }
      .Y[.Y_missing]   <<- 0
      .residual        <<- .Y
      .X_for_Y_missing <<- array(.X,dim = c(.N,.J,.R))
      for (r in 1:.R)
        .X_for_Y_missing[.Y_missing[,r],,r] <<- as.numeric(NA)
      return(invisible(self))
    },
      
    # This returns the R6 object invisibly.
    set_residual_variance = function(residual_variance=NULL, numeric = FALSE,
                                     quantities = c("residual_variance","effect_variance")){
      if('residual_variance' %in% quantities){
        if (is.null(residual_variance)) {
          if (.R > 1) {
            stop("Unspecified residual_variance is not allowed in the presence of missing data in Y")
          }
          else residual_variance = var(.Y[.Y_non_missing], na.rm = TRUE)
        }
        if (is.matrix(residual_variance)) {
          if (anyNA(diag(residual_variance)))
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
        .residual_variance <<- residual_variance
        .residual_variance_inv <<- list()
        .residual_variance_eigen <<- list()
        for(k in 1:nrow(.missing_pattern)){
          if(.R == 1){
            .residual_variance_inv[[k]] <<- .missing_pattern[k,] / .residual_variance
            if(sum(.missing_pattern[k,])>0){
              .residual_variance_eigen[[k]] <<- .residual_variance
            }else{
              .residual_variance_eigen[[k]] <<- numeric(0)
            }
          }else{
            .residual_variance_inv[[k]] <<- matrix(0, .R, .R)
            if(sum(.missing_pattern[k,])>0){
              Vk = .residual_variance[which(.missing_pattern[k,]), which(.missing_pattern[k,])]
              eigenVk = eigen(Vk, symmetric = TRUE)
              dinv = 1/(eigenVk$values)
              .residual_variance_eigen[[k]] <<- eigenVk$values
              .residual_variance_inv[[k]][which(.missing_pattern[k,]), which(.missing_pattern[k,])] <<- eigenVk$vectors %*% (dinv * t(eigenVk$vectors))
            }else{
              .residual_variance_eigen[[k]] <<- numeric(0)
            }

          }
        }
        if(is.matrix(residual_variance)){
          .residual_correlation <<- cov2cor(residual_variance)
        }else{
          .residual_correlation <<- 1
        }
      }
      if('effect_variance' %in% quantities){
        .svs_inv <<- list()
        .svs <<- list()
        for(j in 1:.J){
          # R by R matrix
          # when there is no missing, it is sum(x_j^2) * V^{-1}
          if(.approximate){
            .svs_inv[[j]] <<- Reduce("+", lapply(1:.N, function(i) t(.residual_variance_inv[[.Y_missing_pattern_assign[i]]] *
                                                                       .X_for_Y_missing[i,j,]) * .X_for_Y_missing[i,j,]))
            .svs[[j]] <<- tryCatch(invert_via_chol(.svs_inv[[j]])$inv, error = function(e){
              invert_via_chol(.svs_inv[[j]] + 1e-8 * diag(.R))$inv} )
          }else{
            if(.R == 1){
              .svs_inv[[j]] <<- Reduce("+", lapply(1:.N, function(i) ((.X_for_Y_missing[i,j,] - .Xbar[j,,])^2) *
                                                     .residual_variance_inv[[.Y_missing_pattern_assign[i]]]))
            }else{
              A1_list = list()
              A2_list = list()
              for(i in 1:.N){
                A1_list[[i]] = t(.residual_variance_inv[[.Y_missing_pattern_assign[i]]] *
                                   .X_for_Y_missing[i,j,]) * .X_for_Y_missing[i,j,]
                A2_list[[i]] = t(t(.residual_variance_inv[[.Y_missing_pattern_assign[i]]]) * .X_for_Y_missing[i,j,])
              }
              A1 = Reduce("+",A1_list)
              A2 = Reduce("+",A2_list)
              Vinvsum = Reduce("+", lapply(1:nrow(.missing_pattern), function(i) .residual_variance_inv[[i]] * sum(.Y_missing_pattern_assign == i)))
              .svs_inv[[j]] <<- A1 - crossprod(.Xbar[j,,], A2) - crossprod(A2, .Xbar[j,,]) + crossprod(.Xbar[j,,], Vinvsum %*% .Xbar[j,,])
            }
            .svs[[j]] <<- tryCatch(invert_via_chol(.svs_inv[[j]])$inv, error = function(e){
              invert_via_chol(.svs_inv[[j]] + 1e-8 * diag(.R))$inv} )
          }
        }
        .is_common_sbhat <<- is_list_common(.svs)
      }

      return(invisible(self))
    },
      
    get_coef = function(use_residual = FALSE){
      if (use_residual) XtY = self$XtR
      else XtY = self$XtY
      bhat = t(sapply(1:.J, function(j) .svs[[j]] %*% XtY[j,]))
      bhat[which(is.nan(bhat))] = 0
      if(.R == 1){
        bhat = t(bhat)
      }
      return(bhat)
    },
      
    # This method returns the R6 object invisibly.
    standardize = function (center, scale) {
        
      # precompute scale
      if(center){
        cm_x = colMeans(.X_for_Y_missing, na.rm = TRUE) # J by R
      }else{
        cm_x = matrix(0, .J, .R) # J by R
      }
      .csd <<- matrix(1, .J, .R)
      for(r in 1:.R){
        if (scale) {
          .csd[,r] <<- colSds(.X_for_Y_missing[,,r], center = cm_x[,r], na.rm = TRUE)
          .csd[.csd[,r]==0, r] <<- 1
        }
        .X_for_Y_missing[,,r] <<- t( (t(.X_for_Y_missing[,,r]) - cm_x[,r]) / .csd[,r] )
        .X_for_Y_missing[,,r][is.na(.X_for_Y_missing[,,r])] <<- 0
      }

      if(.approximate){
        .cm <<- cm_x
        if(center){
          # center Y
          if (.R == 1) .Y_mean <<- mean(.Y[.Y_non_missing])
          else .Y_mean <<- sapply(1:.R, function(r) mean(.Y[.Y_non_missing[,r],r]))
          .Y <<- t(t(.Y) - .Y_mean)
          .Y[!.Y_non_missing] <<- 0
        }
      }else{ # exact computation
        if(center){
          # sum_i V_i^{-1} R by R matrix
          Vinvsum = Reduce("+", lapply(1:nrow(.missing_pattern), function(i)
                                   .residual_variance_inv[[i]] * sum(.Y_missing_pattern_assign == i)))
          .Vinvsuminv <<- invert_via_chol(Vinvsum)$inv

          # sum_i V_i^{-1} y_i R by 1 matrix
          Ysum = Reduce("+", lapply(1:.N, function(i)
                                   .residual_variance_inv[[.Y_missing_pattern_assign[i]]] %*% .Y[i,] ))

          # center Y
          .Y_mean <<- as.numeric(.Vinvsuminv %*% Ysum)
          .Y <<- t(t(.Y) - .Y_mean)
          .Y[!.Y_non_missing] <<- 0

          # center X
          .Xbar <<- array(0, dim=c(.J, .R, .R))
          for(j in 1:.J){
            # For variant j, Vinvsuminv sum_i V_i^{-1} X_{i,j} R by R matrix
            .Xbar[j,,] <<- .Vinvsuminv %*% Reduce("+", lapply(1:.N, function(i) t(t(.residual_variance_inv[[.Y_missing_pattern_assign[i]]]) * .X_for_Y_missing[i,j,]) ))
          }
        }
      }

      return(invisible(self))
    },
      
    compute_Xb = function(b) {
      if(is.vector(b)){
        b = matrix(b, length(b),1)
      }
      if(.approximate){
        Xb = sapply(1:.R, function(r) .X_for_Y_missing[,,r] %*% b[,r])
      }else{
        Xbarb = Reduce("+", lapply(1:.J, function(j) .Xbar[j,,] %*% b[j,]))
        Xb = sapply(1:.R, function(r) .X_for_Y_missing[,,r] %*% b[,r]) - matrix(Xbarb, .N, .R,byrow = TRUE)
      }
      if (nrow(Xb) != .N)
        Xb = t(Xb)
      return(Xb)
    },
      
    rescale_coef = function (b) {
      if(.R == 1){
        .csd <<- as.vector(.csd)
        .cm <<- as.vector(.cm)
      }
      coefs = b/.csd
      if(.approximate){
        if (is.null(dim(coefs))) {
          if (!is.null(.Y_mean)){
            intercept = .Y_mean - sum(.cm * coefs)}
          else intercept = 0
          return(c(intercept, coefs))
        } else {
          if (!is.null(.Y_mean)) intercept = .Y_mean - colSums(.cm * coefs)
          else intercept = 0
          mat = as.matrix(rbind(intercept, coefs))
          rownames(mat) = NULL
          return(mat)
        }
      }else{
        # Length R
        # sum_i V_i^{-1} b^T X[i,]
        D = Reduce('+', lapply(1:.N, function(i) .residual_variance_inv[[.Y_missing_pattern_assign[i]]] %*% crossprod(coefs, .X[i,])))
        if (!is.null(.Y_mean)) intercept = .Y_mean - .Vinvsuminv %*% D
        else intercept = 0
        if(is.null(dim(coefs))){
          return(c(intercept,coefs))
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
      # J by R matrix
      if (is.null(.XtY)){
        # V_i^(-1) y_i = z_i
        VinvY = t(sapply(1:.N, function(i) .residual_variance_inv[[.Y_missing_pattern_assign[i]]] %*% .Y[i,])) # N by R
        if(.R == 1) VinvY = t(VinvY)
        if(.approximate){
          .XtY <<- t(sapply(1:.J, function(j) colSums(.X_for_Y_missing[,j,] * VinvY) ))
        }else{
          # sum_{i=1}^N (diag(X_for_Y_missing[i,j,]) - Xbar[j,,])^T z_i
          VinvYcolsum = colSums(VinvY)
          .XtY <<- t(sapply(1:.J, function(j) colSums(.X_for_Y_missing[,j,] * VinvY) - crossprod(.Xbar[j,,], VinvYcolsum) ))
        }
      }
      if (.R == 1) .XtY <<- t(.XtY)
      return(.XtY)
    },
      
    XtX = function() {
        
      # FIXME: not sure how to compute XtX with missing data.
      if (is.null(.XtX))
        .XtX <<- cor(.X)
      return(.XtX)
    },
      
    XtR = function() {
      # V_i^(-1) r_i = z_i
      VinvR = t(sapply(1:.N, function(i) .residual_variance_inv[[.Y_missing_pattern_assign[i]]] %*% .residual[i,])) # N by R
      if(.R == 1) VinvR = t(VinvR)
      if(.approximate){
        res = t(sapply(1:.J, function(j) colSums(.X_for_Y_missing[,j,] * VinvR) ))
      }else{
        # sum_{i=1}^N (diag(X_for_Y_missing[i,j,]) - Xbar[j,,])^T z_i
        res = t(sapply(1:.J, function(j) colSums(.X_for_Y_missing[,j,] * VinvR) - crossprod(.Xbar[j,,], colSums(VinvR)) ))
      }
      if(nrow(res) != .J)
        res = t(res)
      return(res)
    },
      
    Y_missing_pattern_assign = function()
      .Y_missing_pattern_assign,
      
    residual_variance_eigenvalues = function()
      .residual_variance_eigen
  ),
                             
  private = list(
    .X_for_Y_missing = NULL,
    .Y_non_missing = NULL,
    .missing_pattern = NULL,
    .Y_missing_pattern_assign = NULL,
    .Vinvsuminv = NULL,
    .approximate = FALSE,
    .Xbar = NULL,
    .residual_variance_eigen = NULL
  )
)
