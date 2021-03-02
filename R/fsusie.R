fsusie = function(X, Y, L=2, K=3, prior_weights=NULL,
                 compute_objective=TRUE,
                 coverage = 0.95, min_abs_corr = 0.5,
                 max_iter=100,tol=1e-3,
                 verbosity=2,track_fit=FALSE) {
  if (is.null(prior_weights)) prior_weights = c(rep(1/ncol(X), ncol(X)))
  else prior_weights = prior_weights / sum(prior_weights)
  # set data object
  # FIXME: need to design object for DWT data
  Y = grove::DWT(Y)
  Y = cbind(Y$C,Y$D)
  data = DenseData$new(X, Y)
  data$set_residual_variance(residual_variance=diag(ncol(Y)), quantities = "residual_variance")
  #
  s = mmbr_core(data, NULL, L, NULL, prior_weights,
            F, F, NULL, 0, F, compute_objective, 1, max_iter, tol, track_fit, verbosity)
  # CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$null_index = -9
    s$sets = susie_get_cs(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s, prior_tol=prior_tol)
    s$null_index = NULL
  }
  return(s)
}