#' @title SUm of Single Effect (SuSiE) Regression of multivariate Y on X
#' @details Performs Bayesian multiple linear regression of Y on X.
#' That is, this function
#' fits the regression model Y= sum_l Xb_l + e, where elements of e are iid N(0,residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero.
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L maximum number of non-zero effects
#' @param prior_variance Can be 1) a vector of length L, or a scalar, for scaled prior variance when Y is univariate (equivalent to `susieR::susie`); 2) a matrix for simple Multivariate regression or 3) a MASH fit that contains an array of prior covariance matrices and their weights
#' @param residual_variance the residual variance (defaults to sample variance of Y)
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
#' @param estimate_prior_method the method used for estimating prior variance.
#' @param s_init a previous susie fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param compute_univariate_zscore if true, outputs z-score from per variable univariate regression
#' @param precompute_covariances if true, precomputes various covariance quantities to speed up computations at the cost of increased memory usage
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
#' res = msusie(X,y,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_pip susie_get_cs
#' @export
msusie = function(X,Y,L=10,
                 prior_variance=0.2,
                 residual_variance=NULL,
                 prior_weights=NULL, null_weight=NULL,
                 standardize=TRUE,intercept=TRUE,
                 estimate_residual_variance=TRUE,
                 estimate_prior_variance=FALSE,
                 estimate_prior_method='optim',
                 compute_objective=FALSE,
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 compute_univariate_zscore = FALSE,
                 precompute_covariances = FALSE,
                 max_iter=100,tol=1e-3,
                 verbose=TRUE,track_fit=FALSE) {
  # FIXME: this main function code needs to be examined and cleaned up
  # to make it more elegant and robust
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if (any(is.na(X))) {
    stop("Input X must not contain missing values.")
  }
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (is.null(prior_weights)) prior_weights = c(rep(1/ncol(X), ncol(X)))
  else prior_weights = prior_weights / sum(prior_weights)
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    prior_weights = c(prior_weights * (1-null_weight), null_weight)
    X = cbind(X,0)
  }

  ptm = proc.time()
  ## BEGIN mmbr code
  data = DenseData$new(X, Y, intercept, standardize)
  # FIXME: this is because of issue #5
  if (data$X_has_missing) stop("Missing data in input matrix X is not allowed at this point.")
  if (is.null(residual_variance)) {
    #if (dim(Y)[2] > 1) residual_variance = diag(apply(Y, 2, function(x) var(x, na.rm=T)))
    # FIXME: either need better initialization method, or just quit on error for unspecified residual variance
    if (dim(Y)[2] > 1) residual_variance = cov(Y, use = "pairwise.complete.obs")
    else residual_variance = var(Y, na.rm=T)
    if (is.numeric(prior_variance) && !is.matrix(prior_variance)) residual_variance = as.numeric(residual_variance)
    residual_variance[which(is.na(residual_variance))] = 0
  }
  if (is.matrix(residual_variance)) mashr:::check_positive_definite(residual_variance)

  # for now the type of prior_variance controls the type of regression 
  if (is.numeric(prior_variance)) {
    if (!(is.null(dim(Y)) || dim(Y)[2] == 1) && !is.matrix(prior_variance))
      stop("prior variance cannot be a number when Y is a multivariate variable.")
    if (is.matrix(prior_variance)) {
      base = BayesianMultivariateRegression
    } else {
      base = BayesianSimpleRegression
      # Here prior variance is scaled prior variance
      prior_variance = prior_variance * residual_variance
    }
  } else {
    # FIXME: check prior_variance is valid MASH object
    if (prior_variance$n_condition != ncol(Y)) stop("Dimension mismatch between input prior covariance and response data.")
    base = MashRegression
    if ((data$Y_has_missing && !is_diag_mat(residual_variance)) || precompute_covariances) prior_variance$precompute_cov_matrices(data, residual_variance)
  }
  if (!estimate_prior_variance) estimate_prior_method = NULL
  # Below are the core computations
  SER_model = SingleEffectModel(base)$new(data$n_effect, residual_variance, prior_variance)
  SuSiE_model = SuSiE$new(SER_model, L, estimate_residual_variance, compute_objective, max_iter, tol, track_pip=track_fit, track_lbf=track_fit)
  if (!is.null(s_init)) SuSiE_model$init_coef(s_init$coef_index, s_init$coef_value, ncol(X), ncol(Y))
  SuSiE_model$fit(data, prior_weights, estimate_prior_method, verbose)
  s = report_susie_model(data, SuSiE_model, estimate_prior_variance)
  ## END new mmbr code
  s$walltime = proc.time() - ptm 

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s)
  }
  ## report z-scores from univariate regression
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0) {
      X = X[,1:(ncol(X)-1)]
    }
    s$z = susieR:::calc_z(X,Y,center=intercept,scale=standardize)
  }
  return(s)
}