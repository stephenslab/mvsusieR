#' @title Compute condition specific posterior inclusion probability.
#' @description This is only relevant when canonical priors are used
#' @param m M&M model
#' @param prior_obj prior mixture object
#' @return P by R matrix of PIP per condition
#' @keywords internal
mvsusie_get_pip_per_condition = function(m, prior_obj) {
  condition_pip = mvsusie_get_alpha_per_condition(m, prior_obj)
  return(do.call(cbind, lapply(1:dim(condition_pip)[3], function(r) apply(condition_pip[,,r], 2, function(x) 1-prod(1-x)))))
}

#' @title Compute condition specific posterior inclusion probability per effect
#' @keywords internal
mvsusie_get_alpha_per_condition = function(m, prior_obj) {
  condition_indicator = do.call(rbind, lapply(1:length(prior_obj$prior_variance$xUlist), function(i) as.integer(diag(prior_obj$prior_variance$xUlist[[i]]) != 0)))
  condition_pip = array(0, dim=dim(m$b1))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r])) {
        condition_pip[,,r] = condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    }
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(condition_pip)
}
