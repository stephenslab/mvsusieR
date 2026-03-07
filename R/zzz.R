# Register S3 methods for susieR's internal generics.
#
# susieR defines S3 generics (compute_kl, check_convergence, etc.) that are
# not exported — they are internal to susieR's IBSS framework. mvsusieR
# provides methods for these generics (dispatched on mv_individual, mv_ss,
# and mvsusie classes). We use registerS3method() in .onLoad() so that R
# knows to dispatch to our methods when susieR calls these generics internally.
#
# This approach:
#   - Does NOT require susieR to export the generics
#   - Does NOT require mvsusieR to use susieR:::
#   - Survives devtools::document() (since it's R code, not NAMESPACE directives)
#
.onLoad <- function(libname, pkgname) {
  susie_ns <- asNamespace("susieR")
  pkg_ns   <- asNamespace(pkgname)

  # Generics with methods for both mv_individual and mv_ss classes
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

  # Generics with methods for the mvsusie model class
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
