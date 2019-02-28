#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_compute_weight = function(value, weight, log = TRUE) {
    if (!log) value = log(value)
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(alpha = as.vector(w_weighted / weighted_sum_w), log_total = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @importFrom abind abind
#' @keywords internal
report_susie_model = function(d, m) {
    if (is.null(dim(m$posterior_b1[[1]]))) {
      mu = t(do.call(cbind, m$posterior_b1))
      mu2 = t(do.call(cbind, m$posterior_b2))
      b = colSums(t(m$pip)*mu)
    } else {
      mu =  aperm(abind::abind(m$posterior_b1,along=3), c(3,1,2))
      mu2 = aperm(abind::abind(m$posterior_b2,along=3), c(3,1,2))
      b = do.call(cbind, lapply(1:dim(mu)[3], function(i) colSums(t(m$pip) * mu[,,i])))
    }
    if (is.null(m$mixture_posterior_weights)) mixture_weights = NA
    else mixture_weights = aperm(abind::abind(m$mixture_posterior_weights,along=3), c(3,1,2))
    s = list(
        alpha = t(m$pip),
        mu = mu,
        mu2 = mu2,
        KL = m$kl,
        lbf = m$lbf,
        sigma2 = m$residual_variance,
        V = m$prior_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$get_niter(),
        fitted = d$fitted,
        coef = d$rescale_coef(b),
        null_index = -9,
        mixture_weights = mixture_weights 
        )
    s$intercept = s$coef[1]
    class(s) = 'susie'
    return(s)
}

#' @title Compute condition specific posterior inclusion probability
#' @param m M&M model
#' @param prior_obj prior mixture object
#' @export
mmbr_get_pip_per_condition = function(m, prior_obj) {
  condition_indicator = do.call(rbind, lapply(1:length(prior_obj$prior_covariance$xUlist), function(i) as.integer(diag(prior_obj$prior_covariance$xUlist[[i]]) != 0)))
  condition_pip = array(0, dim=dim(m$mu))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r])) {
        condition_pip[,,r] = condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    }
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(do.call(cbind, lapply(1:dim(condition_pip)[3], function(r) apply(condition_pip[,,r], 2, function(x) 1-prod(1-x)))))
}

#' @title A null progressbar, because currently `progressbar_enabled` feature does not work for `progress_bar`
#' @importFrom R6 R6Class
#' @keywords internal
null_progress_bar = R6Class('null_progress_bar', public = list(tick = function(...) {}))

#' @title check if all elements are the same in matrix of J by R, J >> R
#' @keywords internal
is_mat_common = function(mat) {
  all((t(mat) - mat[1,]) == 0)
}

#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param s number of effect variables per condition if greater than 1; otherwise percentage of effect variables per condition
#' @param center_scale FALSE by default
#' @export
mmbr_sim1 = function(n=200,p=500,r=2,s=4,center_scale=FALSE) {
  X = matrix(rnorm(n*p,0,1),n,p)
  if (s>1) {
    beta = matrix(0, p, r)
    for (i in 1:r) beta[sample(1:p,s), i] = 1
  } else {
    beta = matrix(runif(p*r)>s, p, r)
  }
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  if (center_scale) {
    X = scale(X)
    y = y - apply(y,2,mean)
  }
  scaled_prior_variance = 0.2
  return(list(X=X,y=y, d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y), b=beta))
}

#' @title Get lfsr for one condition
#' @importFrom stats pnorm
#' @keywords internal
mmbr_get_one_lfsr = function(mu, mu2, alpha) {
    pos_prob = pnorm(0,mean=t(mu),sd=sqrt(mu2-mu^2))
    neg_prob = 1 - pos_prob
    pmax(0, 1 - rowSums(alpha * t(pmax(pos_prob,neg_prob))))
}

#' @title Local false sign rate (lfsr) for credible sets
#' @details This computes the lfsr of CS identified for each condition.
#' @param m a mmbr fit, the output of `mmbr::susie()`
#' @return a L by R matrix of lfsr
#' @export
mmbr_get_lfsr = function(m) {
    do.call(cbind, lapply(1:dim(m$mu)[3], function(i) mmbr_get_one_lfsr(m$mu[,,i], m$mu2[,,i], m$alpha)))
}