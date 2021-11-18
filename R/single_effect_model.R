# Single effect model object; it is meant to be dynamically inherited
# from any regression model.
#
#' @importFrom R6 R6Class
SingleEffectModel <- function(base)
  R6Class("SingleEffectModel",
  inherit = base,
  public = list(
    initialize = function (J, prior_variance) {
      super$initialize(J,prior_variance)
      private$.pip = rep(0,J)
  },
      
        fit = function(d, prior_weights=NULL, estimate_prior_variance_method=NULL, check_null_threshold=0, save_var=FALSE) {
            if (is.null(prior_weights)) prior_weights = rep(1/private$J, private$J)
            super$fit(d, use_residual=TRUE, prior_weights=prior_weights, estimate_prior_variance_method=estimate_prior_variance_method, save_var=save_var, check_null_threshold=ifelse(is.na(check_null_threshold), 0, check_null_threshold))
            ws = compute_softmax(private$.lbf, prior_weights, log = TRUE)
            private$.pip = ws$weights
            private$lbf_single_effect = ws$log_sum
            if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method == "EM") {
                V = private$estimate_prior_variance_em(private$.pip)
                # when check_null_threshold = NA we skip this check with zero estimate
                # see details in https://github.com/stephenslab/mvsusieR/issues/26
                if (!is.na(check_null_threshold)) {
                    if (private$loglik(0,private$cache$b,private$cache$s,prior_weights) + check_null_threshold >= 
                        private$loglik(V,private$cache$b,private$cache$s,prior_weights)) {
                        V=0
                        # set the corresponding posterior also to zero
                        private$.posterior_b1 = private$.posterior_b1 * 0
                        private$.posterior_b2 = private$.posterior_b2 * 0
                        private$.lbf = private$.lbf * 0
                    }
                }
                private$prior_variance_scalar = V
                private$cache = NULL
            }
        },
        predict = function(d) {
            d$compute_Xb(self$posterior_b1)
        },
        compute_kl = function(d) {
            # compute KL divergence
            pp_eloglik = private$compute_expected_loglik_partial(d)
            private$.kl = pp_eloglik - private$lbf_single_effect
        }
    ),
          
    private = list(
        .pip = NULL, # posterior inclusion probability, alpha
        lbf_single_effect = NULL, # logBF for SER model: sum of logBF of all J effects
        .kl = NULL,
        .vbxxb = NULL,
        .bxxb = NULL,
        
        # This is the expected loglik minus the loglik_null, that is, N(R|B,V) - N(R|0,V)
        compute_expected_loglik_partial = function(d) {
            if (is.matrix(d$residual_variance)) {
                private$compute_expected_loglik_partial_multivariate(d)
            } else {
                private$compute_expected_loglik_partial_univariate(d)
            }
        },
        
        compute_expected_loglik_partial_univariate = function(d) {
            if(d$Y_has_missing){
              Xb = d$compute_Xb(self$posterior_b1)
              resid_var_inv = unlist(d$residual_variance_inv)[d$Y_missing_pattern_assign]
              E1 = sum(d$residual * Xb * resid_var_inv)
              private$.vbxxb = sum(self$posterior_b2*unlist(d$svs_inv))
              return(E1 - (private$.vbxxb / 2))
            }else{
              return(- (0.5/d$residual_variance) * (- 2*sum(self$posterior_b1*d$XtR) + sum(d$X2_sum*as.vector(self$posterior_b2))))
            }
        },
        
    # Posterior variance covariance matrix, weighted by PIP.
    compute_expected_loglik_partial_multivariate = function (d) {
      if (length(dim(private$.posterior_b2)) == 3)
        pb2 = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * private$.posterior_b2[,,j])
      else
        pb2 = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * matrix(private$.posterior_b2[j,]))
      if (d$Y_has_missing) {
        Xb = d$compute_Xb(self$posterior_b1)
        E1 = sum(sapply(1:d$n_sample, function(i) crossprod(d$residual[i,], d$residual_variance_inv[[d$Y_missing_pattern_assign[i]]] %*% Xb[i,])))
        private$.vbxxb = sum(sapply(1:length(pb2), function(j) tr(d$svs_inv[[j]] %*% pb2[[j]])))
        return(E1 - (private$.vbxxb/2))
      } else {
        E1 = crossprod(self$posterior_b1, d$XtR)
        E1 = tr(d$residual_variance_inv %*% (E1 + t(E1)))
        private$.vbxxb = sum(d$X2_sum * sapply(1:length(pb2), function(j) tr(d$residual_variance_inv %*% pb2[[j]])))
              private$.bxxb = Reduce("+", lapply(1:length(pb2), function(j) d$X2_sum[j] * pb2[[j]]))
              return((E1 - private$.vbxxb) / 2)
            }
        }
    ),
          
  active = list(
      
    # Allow for initialization of coefficients.
    mu = function (v) {
      if (missing(v))
        private$.posterior_b1
      else
        private$.posterior_b1 = v
    },
      
    # Posterior first moment, alpha * posterior_b1_reg.
    posterior_b1 = function() private$.pip * private$.posterior_b1,
      
    # Posterior second moment, alpha * posterior_b2_reg.
    posterior_b2 = function() {
      if (length(dim(private$.posterior_b2)) == 3)
        b2 = t(apply(private$.posterior_b2, 3, diag))
      else
        b2 = private$.posterior_b2
      return(private$.pip * b2)
    },
      
    pip = function (v) {
      if (missing(v))
        private$.pip
      else
        private$.pip = v
    },

    lbf   = function() private$lbf_single_effect,
    kl    = function() private$.kl,
    vbxxb = function() private$.vbxxb,
    bxxb  = function() private$.bxxb
  )
)
