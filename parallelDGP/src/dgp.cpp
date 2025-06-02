#include <Rcpp.h>
#include <RcppEigen.h>

#include <omp.h>

#include "simulate.hpp"
#include "estimate.hpp"
#include "forecast.hpp"
#include "value_at_risk.hpp"
#include "expected_shortfall.hpp"

Eigen::MatrixXd dgp_norm(
    const Eigen::Ref<const Eigen::VectorXd>& returns,
    int n,
    const Eigen::Ref<const Eigen::VectorXi>& rolling_window,
    double risk_level
)
{
    int n_returns = returns.size();
    int n_rw = rolling_window.size();

    // out columns:
    // sigma_hat_rolling_window(0), var_hat_rolling_window(0), es_hat_rolling_window(0),
    // ...,
    // sigma_hat_rolling_window(n_rw - 1), var_hat_rolling_window(n_rw - 1), es_hat_rolling_window(n_rw - 1)
    Eigen::MatrixXd out(n, 3 * n_rw);

    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic, 1)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n_rw; ++j) {
            int start_pos = n_returns - n - rolling_window(j);
            const Eigen::Ref<const Eigen::VectorXd> returns_rw = returns.segment(start_pos + i, rolling_window(j));
            Eigen::VectorXd estimate = estimate_norm(returns_rw);
            out(i, 3 * j) = forecast_one_step_ahead(estimate, returns_rw);
            out(i, 3 * j + 1) = value_at_risk_norm(risk_level, out(i, 3 * j));
            out(i, 3 * j + 2) = expected_shortfall_norm(risk_level, out(i, 3 * j), out(i, 3 * j + 1));
        }
    }

    return out;
}

extern "C" SEXP R2Cpp_dgp_norm(SEXP R2Cpp_returns, SEXP R2Cpp_n, SEXP R2Cpp_rolling_window, SEXP R2Cpp_risk_level) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());
    int n = Rcpp::as<int>(R2Cpp_n);
    Rcpp::IntegerVector Rcpp_rolling_window = Rcpp::as<Rcpp::IntegerVector>(R2Cpp_rolling_window);
    Eigen::Map<const Eigen::VectorXi> rolling_window(Rcpp_rolling_window.begin(), Rcpp_rolling_window.size());
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);

    Eigen::MatrixXd out = dgp_norm(returns, n, rolling_window, risk_level);
    return Rcpp::wrap(out);
}

Eigen::MatrixXd dgp_std(
    const Eigen::Ref<const Eigen::VectorXd>& returns,
    int n,
    const Eigen::Ref<const Eigen::VectorXi>& rolling_window,
    double risk_level
)
{
    int n_returns = returns.size();
    int n_rw = rolling_window.size();

    // out columns:
    // sigma_hat_rolling_window(0), var_hat_rolling_window(0), es_hat_rolling_window(0),
    // ...,
    // sigma_hat_rolling_window(n_rw - 1), var_hat_rolling_window(n_rw - 1), es_hat_rolling_window(n_rw - 1)
    Eigen::MatrixXd out(n, 3 * n_rw);

    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic, 1)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n_rw; ++j) {
            int start_pos = n_returns - n - rolling_window(j);
            const Eigen::Ref<const Eigen::VectorXd> returns_rw = returns.segment(start_pos + i, rolling_window(j));
            Eigen::VectorXd estimate = estimate_std(returns_rw);
            out(i, 3 * j) = forecast_one_step_ahead(estimate, returns_rw);
            out(i, 3 * j + 1) = value_at_risk_std(risk_level, out(i, 3 * j), estimate(3));
            out(i, 3 * j + 2) = expected_shortfall_std(risk_level, out(i, 3 * j), out(i, 3 * j + 1), estimate(3));
        }
    }

    return out;
}

extern "C" SEXP R2Cpp_dgp_std(SEXP R2Cpp_returns, SEXP R2Cpp_n, SEXP R2Cpp_rolling_window, SEXP R2Cpp_risk_level) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());
    int n = Rcpp::as<int>(R2Cpp_n);
    Rcpp::IntegerVector Rcpp_rolling_window = Rcpp::as<Rcpp::IntegerVector>(R2Cpp_rolling_window);
    Eigen::Map<const Eigen::VectorXi> rolling_window(Rcpp_rolling_window.begin(), Rcpp_rolling_window.size());
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);

    Eigen::MatrixXd out = dgp_std(returns, n, rolling_window, risk_level);
    return Rcpp::wrap(out);
}
