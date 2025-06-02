#' @export
estimate_norm <- function(returns) {
    return(.Call("R2Cpp_estimate_norm", returns))
}

#' @export
estimate_std <- function(returns) {
    return(.Call("R2Cpp_estimate_std", returns))
}
