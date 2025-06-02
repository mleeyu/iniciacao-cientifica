#' @export
expected_shortfall_norm <- function(risk_level, sigma, value_at_risk) {
    return(.Call("R2Cpp_expected_shortfall_norm", risk_level, sigma, value_at_risk))
}

#' @export
expected_shortfall_std <- function(risk_level, sigma, value_at_risk, nu) {
    return(.Call("R2Cpp_expected_shortfall_std", risk_level, sigma, value_at_risk, nu))
}
