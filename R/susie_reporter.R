#' @title SuSiE reporter object
#' @importFrom R6 R6Class
#' @keywords internal
SuSiEReporter <- R6Class("SuSiEReporter",
  public = list(
    cs = NULL,
    purity = NULL,
    cs_index = NULL,
    pip = NULL,
    annotations = NULL,
    initialize = function(susie_model) {
        private$model = susie_model
        annotations = list()
    },
    comp_cs = function(coverage = 0.95, min_abs_corr = 0.5, X = NULL, Xcorr = NULL) {
        self$cs = NULL
    },
    comp_pip = function() {
        self$pip = NULL
    },
    annotate = function(name, value) {
        self$annotations[[name]] = value
    }
  ), 
  private = list(
    model = NULL
  )
)