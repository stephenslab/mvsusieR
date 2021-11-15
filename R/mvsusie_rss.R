#' @title SUm of Single Effect (SuSiE) Regression using Summary Statistics Z and R
#' @param Z a J by R matrix of z scores
#' @param R a J by J LD matrix
#' @param L maximum number of non-zero effects
#' @param prior_variance Can be 1) a vector of length L, or a scalar, for scaled prior variance when Y is univariate (equivalent to `susieR::susie`); 2) a matrix for simple Multivariate regression or 3) a MASH fit that contains an array of prior covariance matrices and their weights
#' @param residual_variance the residual variance (defaults to 1)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param estimate_prior_variance indicates whether to estimate prior (currently only works for univariate Y and for multivariate Y when prior is a single matrix)
#' @param estimate_prior_method the method used for estimating prior variance: "optim", "uniroot" and "em" for univariate Y, "optim" and "simple" for multivariate Y.
#' @param check_null_threshold when prior variance is estimated, compare the estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example, you can set it to 0.1 to nudge the estimate towards zero. When used with "EM" method
#' setting \code{check_null_threshold=NA} will skip the check and instead relying solely on EM to update this parameter.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param s_init a previous model fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param n_thread maximum number of threads to use for parallel computation (only applicable to mixture prior)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbosity set to 0 for no message output, 1 for a concise progress bar massage output and 2 for one line of message at the end of each iteration.
#' @param track_fit add an attribute \code{trace} to output that saves some current quantities of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{b1}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{b2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{KL}{an L vector of KL divergence}
#' \item{lbf}{an L vector of logBF}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{niter}{number of iterations took for convergence}
#' \item{convergence}{convergence status}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{walltime}{records runtime of the fitting algorithm}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' R = t(X) %*% X
#' z = susieR:::calc_z(X,y)
#' res = mvsusie_rss(z,R,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_cs 
#' @export
mvsusie_rss = function(Z,R,L=10,
                      prior_variance=50,
                      residual_variance=NULL,
                      prior_weights=NULL,
                      estimate_prior_variance=TRUE,
                      estimate_prior_method='EM',
                      check_null_threshold=0, prior_tol=1e-9,
                      compute_objective=TRUE,
                      precompute_covariances = FALSE,
                      s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                      n_thread=1, max_iter=100,tol=1e-3,
                      verbosity=2,track_fit=FALSE) {
  is_numeric_prior = !(is.matrix(prior_variance) || inherits(prior_variance, 'MashInitializer'))
  if (!is.null(dim(Z)) && ncol(Z) > 1 && is_numeric_prior) stop("Please specify prior variance for the multivariate z-scores")
  
  if (any(is.na(Z))) {
    warning('NA values in Z-scores are replaced with 0.')
    Z[is.na(Z)] = 0
  }
  is_numeric_matrix(R,'R')
  
  if(is.null(dim(Z))){
    YtY = 1
  }else{
    YtY = diag(ncol(Z))
  }
  
  s = mvsusie_suff_stat(XtX=R, XtY=Z, YtY=YtY, N=2, L=L,
                       prior_variance=prior_variance,
                       residual_variance=residual_variance,
                       prior_weights=prior_weights,
                       standardize=FALSE,
                       estimate_residual_variance=FALSE,
                       estimate_prior_variance=estimate_prior_variance,
                       estimate_prior_method=estimate_prior_method,
                       check_null_threshold=check_null_threshold, prior_tol=prior_tol,
                       compute_objective=compute_objective,
                       precompute_covariances = precompute_covariances,
                       s_init = s_init,coverage=coverage,min_abs_corr=min_abs_corr,
                       n_thread=n_thread, max_iter=max_iter,tol=tol,
                       verbosity=verbosity,track_fit=track_fit)
  return(s)
}
