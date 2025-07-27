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
    int n_oos,
    const Eigen::Ref<const Eigen::VectorXi>& rws,
    double risk_level
)
{
    int n_returns = returns.size();
    int n_rws = rws.size();

    // out columns:
    // sigma_hat_rws(0), var_hat_rws(0), es_hat_rws(0),
    // ...,
    // sigma_hat_rws(n_rws - 1), var_hat_rws(n_rws - 1), es_hat_rws(n_rws - 1)
    Eigen::MatrixXd out(n_oos, 3 * n_rws);

    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic, 1)
    for (int i = 0; i < n_oos; ++i) {
        for (int j = 0; j < n_rws; ++j) {
            int start_pos = n_returns - n_oos - rws(j);
            const Eigen::Ref<const Eigen::VectorXd> returns_rw = returns.segment(start_pos + i, rws(j));
            Eigen::VectorXd estimate = estimate_norm(returns_rw);
            out(i, 3 * j) = forecast_one_step_ahead(estimate, returns_rw);
            out(i, 3 * j + 1) = value_at_risk_norm(risk_level, out(i, 3 * j));
            out(i, 3 * j + 2) = expected_shortfall_norm(risk_level, out(i, 3 * j), out(i, 3 * j + 1));
        }
    }

    return out;
}

extern "C" SEXP R2Cpp_dgp_norm(SEXP R2Cpp_returns, SEXP R2Cpp_n_oos, SEXP R2Cpp_rws, SEXP R2Cpp_risk_level) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());
    int n_oos = Rcpp::as<int>(R2Cpp_n_oos);
    Rcpp::IntegerVector Rcpp_rws = Rcpp::as<Rcpp::IntegerVector>(R2Cpp_rws);
    Eigen::Map<const Eigen::VectorXi> rws(Rcpp_rws.begin(), Rcpp_rws.size());
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);

    Eigen::MatrixXd out = dgp_norm(returns, n_oos, rws, risk_level);
    return Rcpp::wrap(out);
}

Eigen::MatrixXd dgp_std(
    const Eigen::Ref<const Eigen::VectorXd>& returns,
    int n_oos,
    const Eigen::Ref<const Eigen::VectorXi>& rws,
    double risk_level
)
{
    int n_returns = returns.size();
    int n_rws = rws.size();

    // out columns:
    // sigma_hat_rws(0), var_hat_rws(0), es_hat_rws(0),
    // ...,
    // sigma_hat_rws(n_rws - 1), var_hat_rws(n_rws - 1), es_hat_rws(n_rws - 1)
    Eigen::MatrixXd out(n_oos, 3 * n_rws);

    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic, 1)
    for (int i = 0; i < n_oos; ++i) {
        for (int j = 0; j < n_rws; ++j) {
            int start_pos = n_returns - n_oos - rws(j);
            const Eigen::Ref<const Eigen::VectorXd> returns_rw = returns.segment(start_pos + i, rws(j));
            Eigen::VectorXd estimate = estimate_std(returns_rw);
            out(i, 3 * j) = forecast_one_step_ahead(estimate, returns_rw);
            out(i, 3 * j + 1) = value_at_risk_std(risk_level, out(i, 3 * j), estimate(3));
            out(i, 3 * j + 2) = expected_shortfall_std(risk_level, out(i, 3 * j), out(i, 3 * j + 1), estimate(3));
        }
    }

    return out;
}

extern "C" SEXP R2Cpp_dgp_std(SEXP R2Cpp_returns, SEXP R2Cpp_n_oos, SEXP R2Cpp_rws, SEXP R2Cpp_risk_level) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());
    int n_oos = Rcpp::as<int>(R2Cpp_n_oos);
    Rcpp::IntegerVector Rcpp_rws = Rcpp::as<Rcpp::IntegerVector>(R2Cpp_rws);
    Eigen::Map<const Eigen::VectorXi> rws(Rcpp_rws.begin(), Rcpp_rws.size());
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);

    Eigen::MatrixXd out = dgp_std(returns, n_oos, rws, risk_level);
    return Rcpp::wrap(out);
}
