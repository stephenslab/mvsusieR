# Register S3 methods for susieR's internal generics and cache dependency
# internals from mashr and susieR.
#
# susieR defines S3 generics (compute_kl, check_convergence, etc.) that are
# not exported. mvsusieR provides methods for these generics (dispatched on
# mv_individual, mv_ss, and mvsusie classes). We use registerS3method() so
# that R dispatches to our methods when susieR calls these generics internally.
#
# mashr and susieR internal functions are cached as package-level bindings
# so we can call them directly without mashr::: or susieR:::.

# Dependency internals -- populated by .onLoad()
# mashr:
calc_lik_rcpp <- NULL
calc_sermix_rcpp <- NULL
compute_null_loglik_from_matrix <- NULL
compute_alt_loglik_from_matrix_and_pi <- NULL
expand_cov <- NULL
check_covmat_basics <- NULL
issemidef <- NULL
check_positive_definite <- NULL
# susieR:
get_var_y <- NULL
initialize_susie_model <- NULL
initialize_fitted <- NULL

.onLoad <- function(libname, pkgname) {
  susie_ns <- asNamespace("susieR")
  pkg_ns   <- asNamespace(pkgname)
  mashr_ns <- asNamespace("mashr")

  # Cache mashr internal functions
  for (fn in c("calc_lik_rcpp", "calc_sermix_rcpp",
               "compute_null_loglik_from_matrix",
               "compute_alt_loglik_from_matrix_and_pi",
               "expand_cov", "check_covmat_basics",
               "issemidef", "check_positive_definite")) {
    assign(fn, get(fn, envir = mashr_ns), envir = pkg_ns)
  }

  # Cache susieR internal generics
  for (fn in c("get_var_y", "initialize_susie_model",
               "initialize_fitted")) {
    assign(fn, get(fn, envir = susie_ns), envir = pkg_ns)
  }

  # Register S3 methods for both mv_individual and mv_ss classes
  mv_generics <- c(
    "ibss_initialize",
    "SER_posterior_e_loglik",
    "calculate_posterior_moments",
    "check_convergence",
    "cleanup_model",
    "compute_kl",
    "compute_residuals",
    "compute_ser_statistics",
    "em_update_prior_variance",
    "get_cs",
    "get_fitted",
    "get_intercept",
    "get_objective",
    "get_scale_factors",
    "get_var_y",
    "get_variable_names",
    "get_zscore",
    "initialize_fitted",
    "initialize_susie_model",
    "loglik",
    "neg_loglik",
    "track_ibss_fit",
    "trim_null_effects",
    "update_derived_quantities",
    "update_fitted_values",
    "update_model_variance",
    "update_variance_components",
    "validate_prior"
  )
  for (g in mv_generics) {
    for (cls in c("mv_individual", "mv_ss")) {
      method_fn <- get(paste0(g, ".", cls), envir = pkg_ns)
      registerS3method(g, cls, method_fn, envir = susie_ns)
    }
  }

  # Register S3 methods for the mvsusie model class
  mvsusie_generics <- c(
    "get_alpha_l",
    "get_posterior_mean_l",
    "get_posterior_mean_sum",
    "get_posterior_moments_l",
    "get_prior_variance_l",
    "set_prior_variance_l"
  )
  for (g in mvsusie_generics) {
    method_fn <- get(paste0(g, ".mvsusie"), envir = pkg_ns)
    registerS3method(g, "mvsusie", method_fn, envir = susie_ns)
  }
}
