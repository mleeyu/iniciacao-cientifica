#' @export
value_at_risk_norm <- function(risk_level, sigma) {
    return(.Call("R2Cpp_value_at_risk_norm", risk_level, sigma))
}

#' @export
value_at_risk_std <- function(risk_level, sigma, nu) {
    return(.Call("R2Cpp_value_at_risk_std", risk_level, sigma, nu))
}
