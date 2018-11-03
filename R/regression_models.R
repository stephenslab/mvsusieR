# Base class for regression model
BaseBayesianRegression <- R6Class("BaseBayesianRegression",
  public = list(
    lbf = NULL, # log Bayes factor
    posterior_b1 = NULL, # posterior first moment
    posterior_b2 = NULL, # posterior second moment
    loglik = NULL, # log-likelihood
    mloglik = NULL, # marginal log-likelihood
    initialize = function() {
      self$set_prior()
    },
    set_prior = function() private$exit(), # set prior
    fit = function(d,j,r=NULL) {
      # d: data object
      # j: the j-th effect of data
      # r: the r-th condition of data
      private$comp_posterior(d,j)
      private$comp_lbf()
    }
  ),
  private = list(
    prior = NULL, # prior on effect size
    comp_lbf = function() private$q(), # compute Bayes factor
    comp_loglik = function() private$q(), # compute loglik
    comp_posterior = function(d,j) {
       # compute posterior
       private$q()
    }
    comp_mloglik = function () private$exit(), # compute marginal loglik
    comp_prior = function() private$exit(), # compute prior in EM updates
    exit = function() {
        warning("Not implemented.")
    }
  )
)