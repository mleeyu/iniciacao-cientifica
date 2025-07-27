#' @export
simulate_norm <- function(params, n, n_bis, risk_level, seed) {
    return(.Call("R2Cpp_simulate_norm", params, n, n_bis, risk_level, seed))
}

#' @export
simulate_std <- function(params, n, n_bis, risk_level, seed) {
    return(.Call("R2Cpp_simulate_std", params, n, n_bis, risk_level, seed))
}

#' @export
simulate_dist <- function(dist, params, n, n_bis, risk_level, seed) {
    switch(dist,
        norm = {
            return(.Call("R2Cpp_simulate_norm", params, n, n_bis, risk_level, seed))
        },
        std = {
            return(.Call("R2Cpp_simulate_std", params, n, n_bis, risk_level, seed))
        },
        {
            stop("Invalid dist.")
        }
    )
}
