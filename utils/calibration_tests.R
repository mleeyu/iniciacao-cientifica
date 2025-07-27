dq <- function(returns, value_at_risk, risk_level, L = 4, M = 1) {
    n <- length(returns)
    hits <- ifelse(returns < value_at_risk, 1 - risk_level, - risk_level)
    min_t <- max(L, M) + 1

    y <- hits[min_t:n]

    lagged_hits <- matrix(nrow = n - (min_t - 1), ncol = L)
    for (i in 1:L) {
        lagged_hits[, i] <- hits[(min_t - i):(n - i)]
    }

    lagged_value_at_risk <- value_at_risk[min_t:n]

    lagged_squared_returns <- matrix(nrow = n - (min_t - 1), ncol = M)
    for (i in 1:M) {
        lagged_squared_returns[, i] <- returns[(min_t - i):(n - i)]^2
    }

    min_nrow <- min(nrow(lagged_hits), length(lagged_value_at_risk), nrow(lagged_squared_returns))
    lagged_hits <- head(lagged_hits, min_nrow)
    lagged_value_at_risk <- head(lagged_value_at_risk, min_nrow)
    lagged_squared_returns <- head(lagged_squared_returns, min_nrow)
    x <- cbind(1, lagged_hits, lagged_value_at_risk, lagged_squared_returns)

    test_statistic <- t(y) %*% x %*% MASS::ginv(t(x) %*% x) %*% t(x) %*% y / (risk_level * (1 - risk_level))
    p_value <- 1 - pchisq(test_statistic, df = ncol(x))

    return(list(test_statistic = test_statistic, p_value = p_value))
}

vqr <- function(returns, value_at_risk, risk_level) {
    fit <- suppressWarnings(
        summary(quantreg::rq(returns ~ value_at_risk, tau = risk_level, method = "fn"),
                method = "fn", se = "nid", cov = TRUE)
    )

    a1 <- fit$coefficients[1]
    a2 <- fit$coefficients[2]
    M <- matrix(nrow = 2, ncol = 1)
    M[1, 1] <- a1
    M[2, 1] <- a2 - 1

    aa <- fit$cov[1, 1]
    bb <- fit$cov[1, 2]
    cc <- fit$cov[2, 1]
    dd <- fit$cov[2, 2]
    icov <- matrix(nrow = 2, ncol = 2)
    icov[2, 1] <- 1 / (bb - aa * dd / cc)
    icov[2, 2] <- 1 / (dd - cc * bb / aa)
    icov[1, 1] <- - icov[2, 1] * dd / cc
    icov[1, 2] <- - icov[2, 2] * bb / aa

    statistic <- (t(M)) %*% icov %*% M

    if (is.na(statistic)) {
        fit <- suppressWarnings(
            summary(quantreg::rq(returns ~ value_at_risk, tau = risk_level, method = "fn"),
                    method = "fn", se = "boot", cov = TRUE)
        )

        a1 <- fit$coefficients[1]
        a2 <- fit$coefficients[2]
        M <- matrix(nrow = 2, ncol = 1)
        M[1, 1] <- a1
        M[2, 1] <- a2 - 1

        aa <- fit$cov[1, 1]
        bb <- fit$cov[1, 2]
        cc <- fit$cov[2, 1]
        dd <- fit$cov[2, 2]
        icov <- matrix(nrow = 2, ncol = 2)
        icov[2, 1] <- 1 / (bb - aa * dd / cc)
        icov[2, 2] <- 1 / (dd - cc * bb / aa)
        icov[1, 1] <- - icov[2, 1] * dd / cc
        icov[1, 2] <- - icov[2, 2] * bb / aa

        statistic <- (t(M)) %*% icov %*% M
    }

    p_value <- 1 - pchisq(statistic[1, 1], df = 2)

    return(p_value)
}

calibration_tests <- function(returns, sigma, value_at_risk, expected_shortfall, risk_level) {
    UC_CC <- GAS::BacktestVaR(data = returns, VaR = value_at_risk, alpha = risk_level)
    DQ <- dq(returns = returns, value_at_risk = value_at_risk, risk_level = risk_level, L = 4, M = 1)
    VQR <- vqr(returns = returns, value_at_risk = value_at_risk, risk_level = risk_level)

    ER <- esback::er_backtest(r = returns, q = value_at_risk, e = expected_shortfall, s = sigma, B = 1000)
    gCoC_sCoC <- esback::cc_backtest(r = returns, q = value_at_risk, e = expected_shortfall, s = sigma, alpha = risk_level)
    aESR <- esback::esr_backtest(r = returns, q = value_at_risk, e = expected_shortfall, alpha = risk_level, version = 2)
    sESR <- esback::esr_backtest(r = returns, q = value_at_risk, e = expected_shortfall, alpha = risk_level, version = 1)
    osiESR_tsiESR <- esback::esr_backtest(r = returns, q = value_at_risk, e = expected_shortfall, alpha = risk_level, version = 3)

    p_values <- c(
        unname(UC_CC$LRuc[2]),
        unname(UC_CC$LRcc[2]),
        DQ$p_value,
        VQR,
        ER$pvalue_onesided_standardized,
        gCoC_sCoC$pvalue_twosided_general,
        gCoC_sCoC$pvalue_twosided_simple,
        aESR$pvalue_twosided_asymptotic,
        sESR$pvalue_twosided_asymptotic,
        osiESR_tsiESR$pvalue_onesided_asymptotic,
        osiESR_tsiESR$pvalue_twosided_asymptotic
    )
    return(p_values)
}

calibration_tests_names <- function() {
    return(c("UC", "CC", "DQ (L=4, M=1)", "VQR", "ER", "gCoC", "sCoC", "aESR", "sESR", "osiESR", "tsiESR"))
}
