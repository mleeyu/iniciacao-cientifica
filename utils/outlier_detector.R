# Wavelet-based detection of outliers in financial time series
# Aurea Grané and Helena Veiga (2010)
# https://doi.org/10.1016/j.csda.2009.12.010

# install.packages("wavelets")
# library(wavelets)

# Threshold for Student-t distribution with variance = 1
# param n_mc Number of iterations in Monte Carlo simulation
# param n_sample Size of sample in each iteration
# param alpha Risk level
# param seed Seed for reproducibility
# param df Degrees of freedom of Student-t distribution
threshold_std <- function(n_mc = 20000, n_sample = 1000, alpha = 0.05, seed = 1, df) {
    set.seed(seed)
    s <- sqrt((df - 2) / df)
    max_std_sample <- rep(0, n_mc)
    for (i in 1:n_mc) {
        std_sample <- rt(n = n_sample, df = df) * s
        dwt_std_sample <- dwt(std_sample, filter = "haar", n.levels = 1, boundary = "reflection")
        max_std_sample[i] <- max(abs(dwt_std_sample@W$W1))
    }
    k <- quantile(max_std_sample, prob = 1 - alpha)
    return(k)
}

# Additive Level Outliers (ALOs) Detector (only for isolated ALOs)
# param x Standardized residuals
# param k Threshold
outlier_detector <- function(x, k) {
    n <- length(x)
    ss <- c() # indexes of d_max

    # Step 1: Apply DWT
    dwt_x <- dwt(x, filter = "haar", n.levels = 1, boundary = "reflection")
    d1 <- dwt_x@W$W1 # first level detail coefficients
    # a1 <- dwt_x@V$V1 # first level approx coefficients

    # Step 2: Threshold k already given

    # check if there is a element in d1 (in abs) greater than k
    while (length(d1[abs(d1) > k]) > 0) {
        # Step 3: Find d_max and s
        d_max <- max(abs(d1[abs(d1) > k])) # should I save d_max with abs or not?
        s <- which(abs(d1) == d_max)[1]
        ss <- c(ss, s)

        # Step 4 and 5: Apply IDWT with d1[s] <- 0
        dwt_x@W$W1[s] <- 0 # d1 tilde
        new_x <- idwt(dwt_x)

        # Step 6: Repeat steps 1 to 5 until every element in d1 (in abs) <= k
        dwt_x <- dwt(new_x, filter = "haar", n.levels = 1, boundary = "reflection")
        d1 <- dwt_x@W$W1
    }

    # Step 6: Sort ss
    ss <- sort(ss)

    # Step 7: Locate outlier indexes
    outlier_indexes <- c()
      for (s in ss) {
          x_bar <- (sum(x) - x[2*s] - x[2*s - 1]) / (n - 2)
          if (abs(x[2*s] - x_bar) > abs(x[2*s - 1] - x_bar)) {
              outlier_indexes <- c(outlier_indexes, 2*s)
          } else {
              outlier_indexes <- c(outlier_indexes, 2*s - 1)
          }
      }

    return(list(ss = ss, outlier_indexes = outlier_indexes))
}

# Additive Level Outliers (ALOs) Corrector (only for isolated ALOs)
# param x Returns
# param ss Indexes d_max used during outlier detection
outlier_corrector <- function(x, ss) {
    dwt_x <- dwt(x, filter = "haar", n.levels = 1, boundary = "reflection")
    dwt_x@W$W1[ss] <- 0
    new_x <- idwt(dwt_x)
    return(new_x)
}
