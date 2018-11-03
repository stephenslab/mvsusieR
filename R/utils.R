#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_comp_weight = function(value, weight, log = FALSE) {
    return(value/sum(value * weight))
}