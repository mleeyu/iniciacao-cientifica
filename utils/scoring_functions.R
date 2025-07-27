al <- function(returns, value_at_risk, expected_shortfall, risk_level) {

  g2 <- function(x) return(-1/x)
  gcal2 <- function(x) return(-log(-x))
  a <- function(x) return(1 - log(1 - risk_level))

  hits <- ifelse(returns <= value_at_risk, 1, 0)

  score <- mean(
    g2(expected_shortfall) * (expected_shortfall - value_at_risk + hits * (value_at_risk - returns) / risk_level)
    - gcal2(expected_shortfall) + a(returns)
  )

  return(score)
}

nz <- function(returns, value_at_risk, expected_shortfall, risk_level) {

  g2 <- function(x) return(0.5 * (-x)^-0.5)
  gcal2 <- function(x) return(-(-x)^0.5)

  hits <- ifelse(returns <= value_at_risk, 1, 0)

  score <- mean(
    g2(expected_shortfall) * (expected_shortfall - value_at_risk + hits * (value_at_risk - returns) / risk_level)
    - gcal2(expected_shortfall)
  )

  return(score)
}

fzg <- function(returns, value_at_risk, expected_shortfall, risk_level) {

  g1 <- function(x) return(x)
  g2 <- function(x) return(exp(x) / (1 + exp(x)))
  gcal2 <- function(x) return(log(1 + exp(x)))
  a <- function(x) return(log(2))

  hits <- ifelse(returns <= value_at_risk, 1, 0)

  score <- mean(
    (hits - risk_level) * g1(value_at_risk) - hits * g1(returns)
    + g2(expected_shortfall) * (expected_shortfall - value_at_risk + hits * (value_at_risk - returns) / risk_level)
    - gcal2(expected_shortfall) + a(returns)
  )

  return(score)
}

as <- function(returns, value_at_risk, expected_shortfall, risk_level) {

  W <- 1
  n <- length(returns)
  while (sum(W * value_at_risk < expected_shortfall) == n) W <- W + 1

  g1 <- function(x) return(-0.5 * W * x^2)
  g2 <- function(x) return(risk_level * x)
  gcal2 <- function(x) return(0.5 * risk_level * x^2)

  hits <- ifelse(returns <= value_at_risk, 1, 0)

  score <- mean(
    (hits - risk_level) * g1(value_at_risk) - hits * g1(returns)
    + g2(expected_shortfall) * (expected_shortfall - value_at_risk + hits * (value_at_risk - returns) / risk_level)
    - gcal2(expected_shortfall)
  )

  return(score)
}

scoring_functions <- function(returns, value_at_risk, expected_shortfall, risk_level) {
    return(c(
        al(returns, value_at_risk, expected_shortfall, risk_level),
        nz(returns, value_at_risk, expected_shortfall, risk_level),
        fzg(returns, value_at_risk, expected_shortfall, risk_level),
        as(returns, value_at_risk, expected_shortfall, risk_level)
    ))
}

scoring_functions_names <- function() {
    return(c("AL", "NZ", "FZG", "AS"))
}
