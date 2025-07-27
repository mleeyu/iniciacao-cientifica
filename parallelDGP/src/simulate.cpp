#include <Rcpp.h>
#include <RcppEigen.h>

#include <random>

#include "value_at_risk.hpp"
#include "expected_shortfall.hpp"

Eigen::MatrixXd simulate_norm(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    int n,
    int n_bis,
    double risk_level,
    int seed
)
{
    int n_total = n + n_bis;

    double omega = params(0);
    double alpha = params(1);
    double beta = params(2);

    // out columns: returns, sigma, var, es
    Eigen::MatrixXd out(n_total, 4);
    out(0, 1) = std::sqrt(omega / (1 - alpha - beta));
    out(0, 0) = 0;

    std::mt19937 rng(seed);
    std::normal_distribution<double> rnorm{0.0, 1.0};
    for (int i = 1; i < n_total; ++i) {
        out(i, 1) = std::sqrt(omega + alpha * std::pow(out(i - 1, 0), 2) + beta * std::pow(out(i - 1, 1), 2));
        out(i, 0) = out(i, 1) * rnorm(rng);
        out(i, 2) = value_at_risk_norm(risk_level, out(i, 1));
        out(i, 3) = expected_shortfall_norm(risk_level, out(i, 1), out(i, 2));
    }

    return out.bottomRows(n_total - n_bis);
}

extern "C" SEXP R2Cpp_simulate_norm(SEXP R2Cpp_params, SEXP R2Cpp_n, SEXP R2Cpp_n_bis, SEXP R2Cpp_risk_level, SEXP R2Cpp_seed) {
    Rcpp::NumericVector Rcpp_params = Rcpp::as<Rcpp::NumericVector>(R2Cpp_params);
    Eigen::Map<const Eigen::VectorXd> params(Rcpp_params.begin(), Rcpp_params.size());
    int n = Rcpp::as<int>(R2Cpp_n);
    int n_bis = Rcpp::as<int>(R2Cpp_n_bis);
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    int seed = Rcpp::as<int>(R2Cpp_seed);

    Eigen::MatrixXd out = simulate_norm(params, n, n_bis, risk_level, seed);
    return Rcpp::wrap(out);
}

Eigen::MatrixXd simulate_std(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    int n,
    int n_bis,
    double risk_level,
    int seed
)
{
    int n_total = n + n_bis;

    double omega = params(0);
    double alpha = params(1);
    double beta = params(2);
    double nu = params(3);
    double s = std::sqrt((nu - 2) / nu);

    // out columns: returns, sigma, var, es
    Eigen::MatrixXd out(n_total, 4);
    out(0, 1) = std::sqrt(omega / (1 - alpha - beta));
    out(0, 0) = 0;

    std::mt19937 rng(seed);
    std::student_t_distribution<double> rt{nu};
    for (int i = 1; i < n_total; ++i) {
        out(i, 1) = std::sqrt(omega + alpha * std::pow(out(i - 1, 0), 2) + beta * std::pow(out(i - 1, 1), 2));
        out(i, 0) = out(i, 1) * rt(rng) * s;
        out(i, 2) = value_at_risk_std(risk_level, out(i, 1), nu);
        out(i, 3) = expected_shortfall_std(risk_level, out(i, 1), out(i, 2), nu);
    }

    return out.bottomRows(n_total - n_bis);
}

extern "C" SEXP R2Cpp_simulate_std(SEXP R2Cpp_params, SEXP R2Cpp_n, SEXP R2Cpp_n_bis, SEXP R2Cpp_risk_level, SEXP R2Cpp_seed) {
    Rcpp::NumericVector Rcpp_params = Rcpp::as<Rcpp::NumericVector>(R2Cpp_params);
    Eigen::Map<const Eigen::VectorXd> params(Rcpp_params.begin(), Rcpp_params.size());
    int n = Rcpp::as<int>(R2Cpp_n);
    int n_bis = Rcpp::as<int>(R2Cpp_n_bis);
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    int seed = Rcpp::as<int>(R2Cpp_seed);

    Eigen::MatrixXd out = simulate_std(params, n, n_bis, risk_level, seed);
    return Rcpp::wrap(out);
}
