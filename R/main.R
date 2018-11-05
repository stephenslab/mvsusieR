#' @title SUm of Single Effect (SuSiE) Regression of Y on X
#' @details Performs Bayesian multiple linear regression of Y on X.
#' That is, this function
#' fits the regression model Y= sum_l Xb_l + e, where elements of e are iid N(0,residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=var(Y)*scaled_prior_variance).
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times). The prior variance on each non-zero element of b is set to be var(Y)*scaled_prior_variance.
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `scaled_prior_variance` specifies the prior on the coefficients of X *after* standardization (if performed).
#' If you do not standardize you may need
#' to think more carefully about specifying
#' `scaled_prior_variance`. Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' Any column of X that has zero variance is not standardized, but left as is.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE). The latter is generally not recommended.
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not fully tested and assessed)
#' @param s_init a previous susie fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param compute_univariate_zscore if true, outputs z-score from per variable univariate regression
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Xr}{an n vector of fitted values, equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res =susie(X,y,L=10)
#' coef(res)
#' plot(y,predict(res))
#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
susie = function(X,Y,L=10,scaled_prior_variance=0.2,residual_variance=NULL,
                 prior_weights=NULL, null_weight=NULL,
                 standardize=TRUE,intercept=TRUE,
                 estimate_residual_variance=TRUE,
                 estimate_prior_variance = FALSE,
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 compute_univariate_zscore = FALSE,
                 max_iter=100,tol=1e-3,
                 verbose=FALSE,track_fit=FALSE) {
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X)*(1-null_weight), ncol(X)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    X = cbind(X,0)
  }
  # FIXME: input check and initialization

  data = DenseData$new(intercept, standardize)
  if (estimate_prior_variance) prior = NULL
  else prior =  scaled_prior_variance * as.numeric(var(Y))
  BR_model = BayesBaysianRegression$new(prior)
  SER_model = SingleEffectRegression$new(BR_model, data$get_n_effect(), prior_weights)
  SuSiE_model = SuSiE$new(SER_model, L, residual_variance, estimate_residual_variance, max_iter, tol, track_pip, track_lbf)
  SuSiE_model.fit(data)
  reporter = SuSiEReporter$new(SuSiE_model)
  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    reporter.comp_cs(coverage=coverage, X=X, min_abs_corr=min_abs_corr)
    reporter.comp_pip()
  }
  ## report z-scores from univariate regression
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0) {
      X = X[,1:(ncol(X)-1)]
    }
    reporter$annotate('z', calc_z(X,Y,centered=intercept))
  }
  return(reporter)
}