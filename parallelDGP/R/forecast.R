#' @export
forecast_past <- function(parameters, returns) {
    return(.Call("R2Cpp_forecast_past", parameters, returns))
}

#' @export
forecast_one_step_ahead <- function(parameters, returns) {
    return(.Call("R2Cpp_forecast_one_step_ahead", parameters, returns))
}
