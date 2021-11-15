#' @title Core mvsusie code
#' @importFrom susieR susie_get_pip 
#' @keywords internal
mvsusie_core = function(data, s_init, L, prior_variance, prior_weights,
            estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
            precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbosity) {
  start_time = proc.time()
  # for now the type of prior_variance controls the type of regression
  if (is.numeric(prior_variance)) {
    n_thread = NULL
    if (data$n_condition > 1 && !is.matrix(prior_variance))
      stop(paste("prior variance cannot be a number for multivariate analysis with", data$n_condition, "response variables."))
    if (is.matrix(prior_variance)) {
      base = BayesianMultivariateRegression
    } else {
      base = BayesianSimpleRegression
      prior_variance = prior_variance * data$residual_variance
    }
  } else {
    if (inherits(prior_variance, 'MashInitializer')) {
      if (prior_variance$n_condition != data$n_condition) stop("Dimension mismatch between input prior covariance and response variable data.")
      base = MashRegression
      if (verbosity > 1) {
        message("Initializing prior object ...")
        message(paste("Number of components in the mixture prior:", prior_variance$n_component))
      }
      estimate_prior_scalar = estimate_prior_variance && !is.null(estimate_prior_method) && estimate_prior_method == 'EM'
      if (estimate_prior_scalar) {
        if (precompute_covariances) {
          warning("precompute_covariances option is disabled when prior variances are to be updated.")
        }
        prior_variance$compute_prior_inv()
      } else {
        if (!precompute_covariances) {
          warning("precompute_covariances option is set to FALSE by default to save memory usage with MASH prior. The computation can be a lot slower as a result. It is recommended that you try setting it to TRUE, see if there is a memory usage issue and only switch back if it is a problem.")
        } else {
          prior_variance$precompute_cov_matrices(data)
        }
      }
    } else if (is.null(prior_variance)) {
      # FIXME: a temporary interface for MS for methylation data
      base = NIGMGRegression
      prior_variance = data$n_condition 
    } else {
      stop("prior_variance should be a scalar for univariate response, or a matrix or MashInitializer object for multivariate response.")
    }
  }
  if (!estimate_prior_variance) estimate_prior_method = NULL
  if (verbosity > 1 && "pryr" %in% installed.packages()) {
    message(paste("Memory used by data object", round(pryr::object_size(data)/1024^3, 3), "GB"))
    message(paste("Memory used by prior object", round(pryr::object_size(prior_variance)/1024^3, 3), "GB"))
  }
  # Below are the core computations
  SER_model = SingleEffectModel(base)$new(data$n_effect, prior_variance)
  if (!is.null(n_thread)) {
    if (n_thread<1) stop("number of threads cannot be smaller than 1")
    SER_model$set_thread(floor(n_thread))
  }
  SuSiE_model = SuSiE$new(SER_model, L, estimate_residual_variance, compute_objective, max_iter, tol, track_pip=track_fit, track_lbf=track_fit, track_prior=track_fit)
  if (!is.null(s_init)) SuSiE_model$init_from(s_init)
  SuSiE_model$fit(data, prior_weights, estimate_prior_method, check_null_threshold, verbosity)
  s = report_susie_model(data, SuSiE_model, estimate_prior_variance)
  s$pip = susie_get_pip(s, prior_tol=prior_tol)
  ## clean up prior object
  if (inherits(prior_variance, "MashInitializer")) prior_variance$remove_precomputed()
  s$walltime = proc.time() - start_time
  return(s)
}
