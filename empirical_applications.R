start <- Sys.time()

library(parallelDGP)
library(wavelets)
source("./utils/calibration_tests.R")
source("./utils/scoring_functions.R")
source("./utils/outlier_detector.R")




# Filter and save stocks data
# stocks <- readxl::read_excel(
#         "~/Downloads/retornos_diarios_b3.xlsx",
#         skip = 3,
#         col_types = c("date", rep("numeric", 512)),
#         na = c("", "-", NA)
#     ) |>
#     dplyr::filter(Data > "2010-01-01" & Data < "2025-01-01") |>
#     dplyr::filter(!dplyr::if_all(dplyr::where(is.numeric), is.na)) |>
#     dplyr::select(dplyr::where(~ !any(is.na(.x)))) |>
#     dplyr::rename_with(~ stringr::str_remove(.x, "^(?s).*prov\n"), -Data)
# save(stocks, file = "./data/stocks.RData")

# Load stocks data
load("./data/stocks.RData")
stock_names <- colnames(stocks)[-1]
n_stocks <- length(stock_names)




# Setup
risk_level <- 0.05
rw <- 2500    # rolling window size
n_oos <- 1000 # out-of-sample size

# Create objects to be saved
hits <- matrix(0,
               nrow = 1,
               ncol = n_stocks,
               dimnames = list("hits", stock_names))
tests <- matrix(0,
                nrow = length(calibration_tests_names()),
                ncol = n_stocks,
                dimnames = list(calibration_tests_names(), stock_names))
scores <- matrix(0,
                 nrow = length(scoring_functions_names()),
                 ncol = n_stocks,
                 dimnames = list(scoring_functions_names(), stock_names))
outlier_indexes <- setNames(vector("list", length(stock_names)), stock_names)
filename <- "./data/emp_apps_alpha5_noos1000.RData"




# Calculate number of hits, p-values of calibration tests, scoring functions and
# detect outlier indexes (in out-of-sample period)
for (i in stock_names) {
    print(paste(i, "stock"))

    # dgp is a matrix with n_oos rows and 3 columns
    # dgp[, 1] = sigmas hat using rw
    # dgp[, 2] = value at risks hat using rw
    # dgp[, 3] = expected shortfalls hat using rw
    dgp <- dgp_std(returns = tail(stocks[, i], rw + n_oos),
                    n_oos = n_oos,
                    rws = rw,
                    risk_level = risk_level)

    # number of hits
    hits[1, i] <- sum(tail(stocks[, i], n_oos) < dgp[, 2])

    # p-values of calibration tests
    tests[, i] <- calibration_tests(returns = tail(stocks[, i], n_oos),
                                    sigma = dgp[, 1],
                                    value_at_risk = dgp[, 2],
                                    expected_shortfall = dgp[, 3],
                                    risk_level = risk_level)

    # scoring functions
    scores[, i] <- scoring_functions(returns = tail(stocks[, i], n_oos),
                                     value_at_risk = dgp[, 2],
                                     expected_shortfall = dgp[, 3],
                                     risk_level = risk_level)

    # detect outlier indexes (in out-of-sample period)
    nu_hat <- estimate_std(tail(stocks[, i], rw + n_oos))[4]
    k <- threshold_std(df = nu_hat)
    standardized_residuals <- tail(stocks[, i], n_oos) / dgp[, 1]
    outlier_indexes[[i]] <- outlier_detector(standardized_residuals, k)$outlier_indexes
}




# Save objects
save(hits, tests, scores, outlier_indexes, file = filename)

warnings()
print(Sys.time() - start)
