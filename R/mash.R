#' @title Compute a list of canonical covariance matrices
#' @param R an integer indicating the number of conditions
#' @param singletons a logical value indicating whether the singleton matrices are computed
#' @param hetgrid a vector of numbers between -1 and 1, each representing the off-diagonal elements of matrices with 1s on the diagonal.
#' If 0 is included, the identity matrix will be returned which corresponds to assuming effects are independent across conditions.
#' IF NULL, these matrices are not computed.
#' @return a list of canonical covariance matrices
#' @details This function computes canonical covariance matrices to be provided to mash
#' @examples
#'  mvsusieR:::create_cov_canonical(3)
#'  mvsusieR:::create_cov_canonical(3, singletons=FALSE)
#'  mvsusieR:::create_cov_canonical(3, hetgrid=NULL)
#' 
#' @keywords internal
#' 
create_cov_canonical <- function(R, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1)){
  mats <- list()
  s_idx <- 0
  nms <- vector()
  ###Singleton matrices
  if(singletons) {
    for(i in 1:R) {
      mats[[i]] <- matrix(0, nrow=R, ncol=R)
      mats[[i]][i, i] <- 1
      nms[i] = paste0('singleton_', i)
    }
    s_idx <- R
  }
  ###Heterogeneity matrices
  if(!is.null(hetgrid)) {
    for(j in 1:length(hetgrid)) {
      mats[[s_idx+j]] <- matrix(1, nrow=R, ncol=R)
      mats[[s_idx+j]][lower.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      mats[[s_idx+j]][upper.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      nms[s_idx+j] = paste0('shared_', j)
    }
  }
  names(mats) = nms
  return(mats)
}

#' @title Create mash prior object
#' @param fitted_g from mashr::mash
#' @param mixture_prior a list of (weights = vector(), matrices = list()) where  matrices is a list of prior matrices and have same length as weights.
#' @param sample_data a list of (X=X,Y=Y,residual_variance=residual_variance,center=T,scale=T) to allow for automatically determine canonical priors with equal weights
#' @param null_weight whether or not to add a weight for null in single effect models. By default it takes the null weight from fitted_g
#' if available. Use `null_weight = 0` to override the behavior.
#' @param use_grid expand mixture by grid values as in MASH (not necessary when prior scalar is estimated)
#' @param weights_tol filter out priors with weights smaller than weights_tol
#' @param max_mixture_len only keep the top priors by weight so that the list of mixture prior is of max_mixture_len.
#' Use `max_mixture_len=-1` to include all input weights after weights_tol filtering. Default is set to use all input prior matrices.
#' @param include_indices postprocess input prior to only include conditions from this indices
#' @param ... other parameters, for mvsusieR:::create_cov_canonical
#' @return mash prior object for use with mvsusie() function
#' @details ...
#'
#' @importFrom stats cov2cor
#' 
#' @export
#' 
create_mash_prior = function(fitted_g = NULL, mixture_prior = NULL, sample_data = NULL,
                             null_weight = NULL, use_grid = FALSE, weights_tol = 1e-10, max_mixture_len = -1,
                             include_indices = NULL, ...) {
  if (sum(is.null(fitted_g), is.null(mixture_prior), is.null(sample_data)) != 2)
    stop("Require one and only one of fitted_g, mixture_prior and sample_data to be not NULL.")
  if (!is.null(fitted_g)) {
    # fitted_g: list(pi=pi_s, Ulist=Ulist, grid=grid, usepointmass=usepointmass)
    for (item in c('pi', 'Ulist', 'grid', 'usepointmass')) {
      if (!(item %in% names(fitted_g))) stop(paste("Cannot find", item, "in fitted_g input"))
    }
    if (fitted_g$usepointmass) {
      prior_weights = fitted_g$pi[-1]
      if (is.null(null_weight)) null_weight = fitted_g$pi[1]
    } else {
      prior_weights = mash$fitted_g$pi
    }
    return(MashInitializer$new(fitted_g$Ulist, fitted_g$grid,
                               prior_weights=prior_weights, null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
  if (!is.null(mixture_prior)) {
    for (item in c('matrices')) {
      if (!(item %in% names(mixture_prior))) stop(paste("Cannot find", item, "in mixture_prior input"))
    }
    if (is.null(mixture_prior$weights)) mixture_prior$weights = rep(1/length(mixture_prior$matrices), length(mixture_prior$matrices))
    if (is.null(null_weight)) null_weight = 0
    return(MashInitializer$new(NULL, NULL, xUlist=mixture_prior$matrices, prior_weights=mixture_prior$weights,
                               null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
  if (!is.null(sample_data)) {
    for (item in c('X', 'Y', 'residual_variance')) {
      if (!(item %in% names(sample_data))) stop(paste("Cannot find", item, "in sample_data input"))
    }
    if (is.null(sample_data$center)) {
      #message("Assuming intercept is fitted (otherwise please set 'sample_data$center=F')", stderr())
      sample_data$center = TRUE
    }
    if (is.null(sample_data$scale)) {
      #message("Assuming X is not yet scaled and will scale X (otherwise please set 'sample_data$scale=F')", stderr())
      sample_data$scale = TRUE
    }
    # compute canonical covariances
    Ulist = create_cov_canonical(ncol(sample_data$Y), ...)
    # compute grid
    if (use_grid) {
      d = DenseData$new(sample_data$X, sample_data$Y)
      d$standardize(sample_data$center, sample_data$scale)
      res = d$get_sumstats(diag(sample_data$residual_variance), cov2cor(sample_data$residual_variance))
      ## Use sqrt(3) giving a coarser grid than mash default in exchange for less walltime
      grid = mashr:::autoselect_grid(list(Bhat=res$bhat, Shat=res$sbhat), sqrt(3))
    } else {
      grid = 1
    }
    comp_len = length(grid) * length(Ulist)
    if (max_mixture_len<comp_len && max_mixture_len>0)
      warning(paste0('Automatically generated uniform mixture prior is of length ', comp_len, ' and is greater than currently specified max_mixture_len ', max_mixture_len, ". Please set max_mixture_len=-1 to allow using all of them (although computational speed will suffer)."))
    if (is.null(null_weight)) null_weight = 0
    return(MashInitializer$new(Ulist, grid,
                               prior_weights=NULL, null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
}
