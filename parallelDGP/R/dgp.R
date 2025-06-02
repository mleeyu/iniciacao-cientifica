#' @export
dgp_norm <- function(returns, n, rolling_window, risk_level) {
    return(.Call("R2Cpp_dgp_norm", returns, n, rolling_window, risk_level))
}

#' @export
dgp_std <- function(returns, n, rolling_window, risk_level) {
    return(.Call("R2Cpp_dgp_std", returns, n, rolling_window, risk_level))
}
