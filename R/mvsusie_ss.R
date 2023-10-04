#' @rdname mvsusie
#' 
#' @param XtX A J x J matrix \eqn{X^TX} in which the columns of \eqn{X}
#'   are centered to have mean zero.
#' 
#' @param XtY A J x R matrix \eqn{X^TY} in which the columns of
#'   \eqn{X} and \eqn{Y} are centered to have mean zero.
#' 
#' @param YtY An R x R matrix \eqn{Y^TY} in which the columns of
#' \eqn{Y} are centered to have mean zero.
#' 
#' @param N The sample size.
#' 
#' @param X_colmeans A vector of length J giving the column means of
#'   \eqn{X}. If it is provided with \code{Y_colmeans}, the intercept is
#'   estimated; otherwise, the intercept is \code{NA}.
#' 
#' @param Y_colmeans A vector of length R giving the column means of
#' \eqn{Y}. If it is provided with \code{X_colmeans}, the intercept is
#'   estimated; otherwise, the intercept is \code{NA}.
#' 
#' @importFrom stats var
#' @importFrom stats cov2cor
#' @importFrom susieR susie_get_cs
#' 
#' @export
#' 
mvsusie_suff_stat = function (XtX, XtY, YtY, N, L = 10, X_colmeans = NULL,
                              Y_colmeans = NULL, prior_variance = 0.2,
                              residual_variance = NULL, prior_weights = NULL,
                              standardize = TRUE,
                              estimate_residual_variance = FALSE,
                              estimate_prior_variance = TRUE,
                              estimate_prior_method = "EM",
                              check_null_threshold = 0, prior_tol = 1e-9,
                              compute_objective = TRUE,
                              precompute_covariances = FALSE, s_init = NULL,
                              coverage = 0.95,min_abs_corr = 0.5,n_thread = 1,
                              max_iter = 100,tol = 1e-3, verbosity = 2,
                              track_fit = FALSE) {
    
  # Adjust prior weights.
  if (is.null(prior_weights))
    prior_weights = c(rep(1/ncol(XtX),ncol(XtX)))
  else
    prior_weights = prior_weights / sum(prior_weights)

  # Check and process prior variance.
  if(inherits(prior_variance,"MashInitializer"))
    prior_variance = prior_variance$clone(deep = TRUE)
  is_numeric_prior = !(is.matrix(prior_variance) ||
                       inherits(prior_variance,"MashInitializer"))
  if (!is.null(dim(YtY)) && nrow(YtY) > 1 && is_numeric_prior)
    stop("Please specify prior variance for the multivariate response Y")
  if (standardize && !is_numeric_prior) {
      
    # Scale prior variance; see
    # https://github.com/stephenslab/mvsusieR/blob/master/
    # inst/prototypes/prior_matrices_scale.ipynb
    sigma = sqrt(diag(YtY)/(N-1))
    sigma = sigma/sqrt(N)
    
    # Make sigma numerically more robust against extreme values.
    if (estimate_prior_variance)
      sigma = sigma/max(sigma)
    if (is.matrix(prior_variance))
      prior_variance = scale_covariance(prior_variance,sigma)
    else
      prior_variance$scale_prior_variance(sigma)
  }
  
  # Set data object.
  data = SSData$new(XtX,XtY,YtY,N,X_colmeans,Y_colmeans)
  
  # Include residual variance in data.
  data$set_residual_variance(residual_variance,
                             numeric = is_numeric_prior,
                             quantities = "residual_variance")
  data$standardize(standardize)

  # Include residual variance in the data object.
  data$set_residual_variance(quantities = "effect_variance")
  
  # Fit the susie model.
  s = mvsusie_core(data,s_init,L,prior_variance,prior_weights,
                   estimate_residual_variance,estimate_prior_variance,
                   estimate_prior_method,check_null_threshold,
                   precompute_covariances,compute_objective,n_thread,
                   max_iter,tol,prior_tol,track_fit,verbosity)
  
  # Compute CSs and PIPs.
  if (!is.null(coverage) && !is.null(min_abs_corr))
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = cov2cor(XtX),
                          min_abs_corr = min_abs_corr)

  # Set row and column names of the outputs, and fix dimensions of
  # outputs if needed.
  s$coef   = drop(s$coef)
  s$fitted = drop(s$fitted)
  s$residual_variance = data$residual_variance
  if (is.null(colnames(XtY)))
    if(is.null(dim(XtY))){
      s$condition_names = 'cond1'
    }else{
      s$condition_names = paste0("cond",1:ncol(XtY))
    }
  else
    s$condition_names = colnames(XtY)
  if (is.null(colnames(XtX)))
    s$variable_names = paste0("var",1:ncol(XtX))
  else
    s$variable_names = colnames(XtX)
  if (length(s$condition_names) == 1) {
    names(s$coef)   = c("(Intercept)",s$variable_names)
    names(s$fitted) = s$variable_names
    rownames(s$b1)  = paste0("l",1:L)
    colnames(s$b2)  = s$variable_names
    rownames(s$b2)  = paste0("l",1:L)
    colnames(s$b2)  = s$variable_names
  } else {
    rownames(s$coef)   = c("(Intercept)",s$variable_names)
    colnames(s$coef)   = s$condition_names
    colnames(s$fitted) = s$condition_names
    rownames(s$fitted) = s$variable_names
    dimnames(s$b1) = list(single_effect = paste0("l",1:L),
                          variable      = s$variable_names,
                          condition     = s$condition_names)
    dimnames(s$b2) = list(single_effect = paste0("l",1:L),
                          variable      = s$variable_names,
                          condition     = s$condition_names)
    rownames(s$residual_variance) = s$condition_names
    colnames(s$residual_variance) = s$condition_names
  }
  names(s$pip)       = s$variable_names
  colnames(s$alpha)  = s$variable_names
  names(s$intercept) = s$condition_names
  rownames(s$alpha)  = paste0("l",1:L)
  class(s) = "mvsusie"
  return(s)
}
