#' @export
simulate_norm <- function(parameters, n, n_bis, risk_level, seed) {
    return(.Call("R2Cpp_simulate_norm", parameters, n, n_bis, risk_level, seed))
}

#' @export
simulate_std <- function(parameters, n, n_bis, risk_level, seed) {
    return(.Call("R2Cpp_simulate_std", parameters, n, n_bis, risk_level, seed))
}
