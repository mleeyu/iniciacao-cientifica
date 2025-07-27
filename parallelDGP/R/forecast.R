#' @export
forecast_past <- function(params, returns) {
    return(.Call("R2Cpp_forecast_past", params, returns))
}

#' @export
forecast_one_step_ahead <- function(params, returns) {
    return(.Call("R2Cpp_forecast_one_step_ahead", params, returns))
}
