#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_compute_weight = function(value, weight, log = TRUE) {
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(alpha = as.vector(w_weighted / weighted_sum_w), log_total = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @keywords internal
report_susie_model = function(d, m) {
    s = list(
        alpha = t(m$pip),
        mu = t(do.call(cbind, m$posterior_b1)),
        mu2 = t(do.call(cbind, m$posterior_b2)),
        KL = m$kl,
        lbf = m$lbf,
        sigma2 = m$residual_variance,
        V = m$prior_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$get_niter(),
        fitted = d$fitted,
        coef = m$coef(d),
        null_index = -9
        )
    s$intercept = s$coef[1]
    class(s) = 'susie'
    return(s)
}

#' @title A null progressbar, because currently `progressbar_enabled` feature does not work for `progress_bar`
#' @importFrom R6 R6Class
#' @keywords internal
null_progress_bar = R6Class('null_progress_bar', public = list(tick = function(...) {}))