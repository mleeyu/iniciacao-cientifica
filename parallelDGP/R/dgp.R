#' @export
dgp_norm <- function(returns, n_oos, rws, risk_level) {
    return(.Call("R2Cpp_dgp_norm", returns, n_oos, rws, risk_level))
}

#' @export
dgp_std <- function(returns, n_oos, rws, risk_level) {
    return(.Call("R2Cpp_dgp_std", returns, n_oos, rws, risk_level))
}

#' @export
dgp_dist <- function(dist, returns, n_oos, rws, risk_level) {
    switch(dist,
        norm = {
            return(.Call("R2Cpp_dgp_norm", returns, n_oos, rws, risk_level))
        },
        std = {
            return(.Call("R2Cpp_dgp_std", returns, n_oos, rws, risk_level))
        },
        {
            stop("Invalid dist.")
        }
    )
}
