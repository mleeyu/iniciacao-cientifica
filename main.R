start <- Sys.time()

# Load library and wrapper functions
library(parallelDGP)
source("backtests.R")
source("scores.R")

# Setup
n_monte_carlo <- 10000
risk_level <- 0.05
n <- 1000
burn_in <- 500
rolling_window <- c(500, 1000, 2500)
n_rw <- length(rolling_window)
max_rw <- max(rolling_window)
norm_parameters <- c(0.01, 0.1, 0.85)
std7_parameters <- c(norm_parameters, 7.0)
power7_parameters <- c(norm_parameters, 7.0)
power15_parameters <- c(norm_parameters, 15.0)
power25_parameters <- c(norm_parameters, 25.0)

# Create output object
distributions <- c("norm", "std7", "power7", "power15", "power25")
output <- list()
for (i in distributions) {
    output$sem_outlier$backtests[[i]] <- matrix(
        0, nrow = length(backtests_names()), ncol = n_rw + 1,
        dimnames = list(backtests_names(), c(rolling_window, "True")))
    output$sem_outlier$scores[[i]] <- matrix(
        0, nrow = length(scores_names()), ncol = n_rw + 1,
        dimnames = list(scores_names(), c(rolling_window, "True")))
    output$com_outlier$backtests[[i]] <- matrix(
        0, nrow = length(backtests_names()), ncol = n_rw + 1,
        dimnames = list(backtests_names(), c(rolling_window, "True")))
    output$com_outlier$scores[[i]] <- matrix(
        0, nrow = length(scores_names()), ncol = n_rw + 1,
        dimnames = list(scores_names(), c(rolling_window, "True")))
}

# Simulations
for (i in 1:n_monte_carlo) {
    print(paste("Simulation", i, "out of", n_monte_carlo))
    set.seed(i)
    outlier_index <- sample(1:n, size = 1) + max_rw

    # each simulation matrix has columns: returns, returns with outlier, sigma, var, es
    simulation <- list()
    for (j in distributions) {
        simulation[[j]] <- matrix(0, max_rw + n, 5)
    }
    simulation$norm[, c(1, 3, 4, 5)] <- simulate_norm(
        parameters = norm_parameters, n = max_rw + n, n_bis = burn_in,
        risk_level = risk_level, seed = i)
    simulation$std7[, c(1, 3, 4, 5)] <- simulate_std(
        parameters = std7_parameters, n = max_rw + n, n_bis = burn_in,
        risk_level = risk_level, seed = i)
    simulation$power7[, c(1, 3, 4, 5)] <- simulate_std(
        parameters = power7_parameters, n = max_rw + n, n_bis = burn_in,
        risk_level = risk_level, seed = i)
    simulation$power15[, c(1, 3, 4, 5)] <- simulate_std(
        parameters = power15_parameters, n = max_rw + n, n_bis = burn_in,
        risk_level = risk_level, seed = i)
    simulation$power25[, c(1, 3, 4, 5)] <- simulate_std(
        parameters = power25_parameters, n = max_rw + n, n_bis = burn_in,
        risk_level = risk_level, seed = i)
    for (j in distributions) {
        simulation[[j]][, 2] <- simulation[[j]][, 1]
        simulation[[j]][, 2][outlier_index] <- simulation[[j]][, 2][outlier_index] + ifelse(simulation[[j]][, 2][outlier_index] > 0, 5, -5)
    }

    # dgp functions returns matrix with columns:
    # sigma_hat_rolling_window[1], var_hat_rolling_window[1], es_hat_rolling_window[1],
    # ...
    # sigma_hat_rolling_window[n_rw], var_hat_rolling_window[n_rw], es_hat_rolling_window[n_rw]
    dgp <- list()
    dgp$sem_outlier$norm <- dgp_norm(
        returns = simulation$norm[, 1], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$sem_outlier$std7 <- dgp_std(
        returns = simulation$std7[, 1], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$sem_outlier$power7 <- dgp_norm(
        returns = simulation$power7[, 1], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$sem_outlier$power15 <- dgp_norm(
        returns = simulation$power15[, 1], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$sem_outlier$power25 <- dgp_norm(
        returns = simulation$power25[, 1], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$com_outlier$norm <- dgp_norm(
        returns = simulation$norm[, 2], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$com_outlier$std7 <- dgp_std(
        returns = simulation$std7[, 2], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$com_outlier$power7 <- dgp_norm(
        returns = simulation$power7[, 2], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$com_outlier$power15 <- dgp_norm(
        returns = simulation$power15[, 2], n = n,
        rolling_window = rolling_window, risk_level = risk_level)
    dgp$com_outlier$power25 <- dgp_norm(
        returns = simulation$power25[, 2], n = n,
        rolling_window = rolling_window, risk_level = risk_level)

    for (j in distributions) {
        for (k in 1:n_rw) {
            output$sem_outlier$backtests[[j]][, k] <- output$sem_outlier$backtests[[j]][, k] +
                backtests(
                    returns = tail(simulation[[j]][, 1], n),
                    sigma = dgp$sem_outlier[[j]][, 3*k - 2],
                    value_at_risk = dgp$sem_outlier[[j]][, 3*k - 1],
                    expected_shortfall = dgp$sem_outlier[[j]][, 3*k],
                    risk_level = risk_level)
            output$sem_outlier$scores[[j]][, k] <- output$sem_outlier$scores[[j]][, k] +
                scores(
                    returns = tail(simulation[[j]][, 1], n),
                    value_at_risk = dgp$sem_outlier[[j]][, 3*k - 1],
                    expected_shortfall = dgp$sem_outlier[[j]][, 3*k], risk_level = risk_level)
            output$com_outlier$backtests[[j]][, k] <- output$com_outlier$backtests[[j]][, k] +
                backtests(
                    returns = tail(simulation[[j]][, 2], n),
                    sigma = dgp$com_outlier[[j]][, 3*k - 2],
                    value_at_risk = dgp$com_outlier[[j]][, 3*k - 1],
                    expected_shortfall = dgp$com_outlier[[j]][, 3*k],
                    risk_level = risk_level)
            output$com_outlier$scores[[j]][, k] <- output$com_outlier$scores[[j]][, k] +
                scores(
                    returns = tail(simulation[[j]][, 2], n),
                    value_at_risk = dgp$com_outlier[[j]][, 3*k - 1],
                    expected_shortfall = dgp$com_outlier[[j]][, 3*k], risk_level = risk_level)
        }
        output$sem_outlier$backtests[[j]][, n_rw + 1] <- output$sem_outlier$backtests[[j]][, n_rw + 1] +
            backtests(
                returns = tail(simulation[[j]][, 1], n),
                sigma = tail(simulation[[j]][, 3], n),
                value_at_risk = tail(simulation[[j]][, 4], n),
                expected_shortfall = tail(simulation[[j]][, 5], n),
                risk_level = risk_level)
        output$sem_outlier$scores[[j]][, n_rw + 1] <- output$sem_outlier$scores[[j]][, n_rw + 1] +
            scores(
                returns = tail(simulation[[j]][, 1], n),
                value_at_risk = tail(simulation[[j]][, 4], n),
                expected_shortfall = tail(simulation[[j]][, 5], n), risk_level = risk_level)
        output$com_outlier$backtests[[j]][, n_rw + 1] <- output$com_outlier$backtests[[j]][, n_rw + 1] +
            backtests(
                returns = tail(simulation[[j]][, 2], n),
                sigma = tail(simulation[[j]][, 3], n),
                value_at_risk = tail(simulation[[j]][, 4], n),
                expected_shortfall = tail(simulation[[j]][, 5], n),
                risk_level = risk_level)
        output$com_outlier$scores[[j]][, n_rw + 1] <- output$com_outlier$scores[[j]][, n_rw + 1] +
            scores(
                returns = tail(simulation[[j]][, 2], n),
                value_at_risk = tail(simulation[[j]][, 4], n),
                expected_shortfall = tail(simulation[[j]][, 5], n), risk_level = risk_level)
    }
}
for (i in distributions) {
    output$sem_outlier$backtests[[i]] <- output$sem_outlier$backtests[[i]] / n_monte_carlo
    output$sem_outlier$scores[[i]] <- output$sem_outlier$scores[[i]] / n_monte_carlo
    output$com_outlier$backtests[[i]] <- output$com_outlier$backtests[[i]] / n_monte_carlo
    output$com_outlier$scores[[i]] <- output$com_outlier$scores[[i]] / n_monte_carlo
}

# Save output object
save(output, file = "./data/output.RData")

warnings()
print(Sys.time() - start)
