# Monte Carlo Simulations
table_body <- function(x, y) {
    n_rows <- nrow(x)
    n_cols <- ncol(x)
    tmp <- matrix(rbind(tests[[i]][[1]], tests[[i]][[2]]),
                  nrow = n_rows, ncol = 2*n_cols)
    out <- matrix("", nrow = n_rows, ncol = 2*n_cols)
    rownames(out) <- rownames(x)

    for (i in 1:n_rows) {
        for (j in 1:(2*n_cols)) {
            out[i, j] <- sprintf("%.3f", round(tmp[i, j], 3))
        }
    }
    return(out)
}

load("../data/sims_risklevel5_noos1000.RData")
dists <- names(tests)

dists
i <- dists[1]
xtable::xtable(table_body(tests[[i]][[1]], tests[[i]][[2]]),
               only.contents=TRUE,
               include.rownames=TRUE,
               include.colnames=FALSE,
               hline.after = NULL)

dists
i <- dists[1]
xtable::xtable(table_body(scores[[i]][[1]], scores[[i]][[2]]),
               only.contents=TRUE,
               include.rownames=TRUE,
               include.colnames=FALSE,
               hline.after = NULL)




# Empirical Applications
load("../data/apps_risklevel5_noos1000.RData")

stock_names <- sort(colnames(tests))
n_rows <- length(stock_names)
n_cols <- 1 + 11 + 1 # hits (1) + tests (11) + outliers (1)
out <- matrix("", nrow = n_rows, ncol = n_cols)
rownames(out) <- stock_names

for (i in stock_names) {
    n_hits <- as.character(hits[, i])
    p_values <- sapply(tests[, i], function(x) sprintf("%.3f", round(x, 3)))
    n_outliers <- as.character(length(outlier_indexes[[i]]))
    out_row <- c(n_hits, p_values, n_outliers)
    names(out_row) <- NULL
    out[i, ] <- out_row
}
xtable::xtable(out,
               only.contents=TRUE,
               include.rownames=TRUE,
               include.colnames=FALSE,
               hline.after = NULL)
