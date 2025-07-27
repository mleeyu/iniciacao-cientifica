start <- Sys.time()

library(parallelDGP)
source("./utils/calibration_tests.R")
source("./utils/scoring_functions.R")




# Setup
n_replicates <- 10000
significance_level <- 0.05
risk_level <- 0.05

n_bis <- 500              # burn-in-sample size
rws <- c(500, 1000, 2500) # rolling window sizes
n_rws <- length(rws)      # number of rolling window sizes
n_ins <- max(rws)         # in-sample size
n_oos <- 1000             # out-of-sample size

# sim = "dist" means it simulates GARCH(1, 1)-dist when simulating:
#   returns (without|with) outlier, sigmas, value at risks and expected shortfalls
# dgp = "dist" means it assumes GARCH(1, 1)-dist when generating:
#   sigmas hat, value at risks hat and expected shortfalls hat for each rolling window size
dists <- list(
    size_norm   = list(sim = "norm", dgp = "norm", params = c(0.01, 0.1, 0.85)),
    size_std7   = list(sim = "std",  dgp = "std",  params = c(0.01, 0.1, 0.85, 7.0)),
    power_std7  = list(sim = "std",  dgp = "norm", params = c(0.01, 0.1, 0.85, 7.0)),
    power_std15 = list(sim = "std",  dgp = "norm", params = c(0.01, 0.1, 0.85, 15.0)),
    power_std25 = list(sim = "std",  dgp = "norm", params = c(0.01, 0.1, 0.85, 25.0))
)

# Create objects to be saved
tests <- list()
scores <- list()
for (i in names(dists)) {
    for (j in c("without_outlier", "with_outlier")) {
        tests[[i]][[j]] <- matrix(0,
                                  nrow = length(calibration_tests_names()),
                                  ncol = n_rws + 1,
                                  dimnames = list(calibration_tests_names(), c(rws, "True")))
        scores[[i]][[j]] <- matrix(0,
                                   nrow = length(scoring_functions_names()),
                                   ncol = n_rws + 1,
                                   dimnames = list(scoring_functions_names(), c(rws, "True")))
    }
}
filename <- "./data/mc_sims_alpha5_noos1000.RData"




# Monte Carlo replicates
for (i in 1:n_replicates) {
    print(paste("Monte Carlo replicate", i, "out of", n_replicates))

    # sims[["dist"]][["(without|with)_outlier"]] is a matrix with (n_ins + n_oos) rows and 4 columns
    # sims[["dist"]][["(without|with)_outlier"]][, 1] = returns (without|with) outlier
    # sims[["dist"]][["(without|with)_outlier"]][, 2] = sigmas
    # sims[["dist"]][["(without|with)_outlier"]][, 3] = value at risks
    # sims[["dist"]][["(without|with)_outlier"]][, 4] = expected shortfalls
    sims <- list()
    for (j in names(dists)) {
        sims[[j]][["without_outlier"]] <- simulate_dist(dist = dists[[j]][["sim"]],
                                                       params = dists[[j]][["params"]],
                                                       n = n_ins + n_oos,
                                                       n_bis = n_bis,
                                                       risk_level = risk_level,
                                                       seed = i)

        set.seed(i)
        outlier_index <- sample(1:n_oos, size = 1) + n_ins
        sims[[j]][["with_outlier"]] <- sims[[j]][["without_outlier"]]
        sims[[j]][["with_outlier"]][outlier_index, 1] <- sims[[j]][["with_outlier"]][outlier_index, 1] +
            sign(sims[[j]][["with_outlier"]][outlier_index, 1]) * 5
    }

    # dgps[["dist"]][["(without|with)_outlier"]] is a matrix with n_oos rows and 3*n_rws columns
    # dgps[["dist"]][["(without|with)_outlier"]][, 1] = sigmas hat using rws[1]
    # dgps[["dist"]][["(without|with)_outlier"]][, 2] = value at risks hat using rws[1]
    # dgps[["dist"]][["(without|with)_outlier"]][, 3] = expected shortfalls hat using rws[1]
    # ...
    # dgps[["dist"]][["(without|with)_outlier"]][, 3*n_rws - 2] = sigmas hat using rws[n_rws]
    # dgps[["dist"]][["(without|with)_outlier"]][, 3*n_rws - 1] = value at risks hat using rws[n_rws]
    # dgps[["dist"]][["(without|with)_outlier"]][, 3*n_rws] = expected shortfall hat using rws[n_rws]
    dgps <- list()
    for (j in names(dists)) {
        for (k in c("without_outlier", "with_outlier")) {
            dgps[[j]][[k]] <- dgp_dist(dist = dists[[j]][["dgp"]],
                                       returns = sims[[j]][[k]][, 1],
                                       n_oos = n_oos,
                                       rws = rws,
                                       risk_level = risk_level)
        }
    }

    for (j in names(dists)) {
        for (k in c("without_outlier", "with_outlier")) {
            # For each rolling window size column
            for (l in 1:n_rws) {
                tests[[j]][[k]][, l] <- tests[[j]][[k]][, l] +
                    (calibration_tests(returns = tail(sims[[j]][[k]][, 1], n_oos),
                                       sigma = dgps[[j]][[k]][, 3*l - 2],
                                       value_at_risk = dgps[[j]][[k]][, 3*l - 1],
                                       expected_shortfall = dgps[[j]][[k]][, 3*l],
                                       risk_level = risk_level) < significance_level)
                scores[[j]][[k]][, l] <- scores[[j]][[k]][, l] +
                    scoring_functions(returns = tail(sims[[j]][[k]][, 1], n_oos),
                                      value_at_risk = dgps[[j]][[k]][, 3*l - 1],
                                      expected_shortfall = dgps[[j]][[k]][, 3*l],
                                      risk_level = risk_level)
            }
            # For the "True" column
            l <- n_rws + 1
            tests[[j]][[k]][, l] <- tests[[j]][[k]][, l] +
                (calibration_tests(returns = tail(sims[[j]][[k]][, 1], n_oos),
                                   sigma = tail(sims[[j]][[k]][, 2], n_oos),
                                   value_at_risk = tail(sims[[j]][[k]][, 3], n_oos),
                                   expected_shortfall = tail(sims[[j]][[k]][, 4], n_oos),
                                   risk_level = risk_level) < significance_level)
            scores[[j]][[k]][, l] <- scores[[j]][[k]][, l] +
                scoring_functions(returns = tail(sims[[j]][[k]][, 1], n_oos),
                                  value_at_risk = tail(sims[[j]][[k]][, 3], n_oos),
                                  expected_shortfall = tail(sims[[j]][[k]][, 4], n_oos),
                                  risk_level = risk_level)
        }
    }
}
for (i in names(dists)) {
    for (j in c("without_outlier", "with_outlier")) {
        tests[[i]][[j]] <- tests[[i]][[j]] / n_replicates
        scores[[i]][[j]] <- scores[[i]][[j]] / n_replicates
    }
}




# Save objects
save(tests, scores, file = filename)

warnings()
print(Sys.time() - start)
