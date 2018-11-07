#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_comp_weight = function(value, weight, log = TRUE) {
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(alpha = w_weighted / weighted_sum_w, log_total = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @keywords internal
report_susie_model = function(d, m) {
    s = list()
    s$fitted = m$fitted
    s$coef = m$coef(d)
    s$intercept = s$coef[1]
    return(s)
}