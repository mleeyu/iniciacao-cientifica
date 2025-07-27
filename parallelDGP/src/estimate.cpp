#include <Rcpp.h>
#include <RcppEigen.h>

#include <RcppGSL.h>
#include <gsl/gsl_statistics.h>
#include <nloptrAPI.h>

#include "forecast.hpp"

double log_likelihood_norm(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
)
{
    int n = returns.size();

    Eigen::VectorXd sigma = forecast_past(params, returns);

    double sigma2;
    double sum_log_likelihood = 0;

    for (int i = 0; i < n; ++i) {
        sigma2 = std::pow(sigma(i), 2);
        sum_log_likelihood += - std::log(sigma2) - std::pow(returns(i), 2) / sigma2;
    }

    return sum_log_likelihood;
}

double objective_norm(unsigned n, const double *x, double *grad, void *data)
{
    Eigen::VectorXd& returns = *(Eigen::VectorXd*) data;

    Eigen::Map<const Eigen::VectorXd> params(x, n);

    double epsilon = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);

    if (grad) {
        Eigen::VectorXd x_forward;
        Eigen::VectorXd x_backward;
        double f_forward;
        double f_backward;

        for (unsigned i = 0; i < n; ++i) {
            x_forward = params;
            x_backward = params;

            x_forward[i] += epsilon;
            x_backward[i] -= epsilon;

            f_forward = - log_likelihood_norm(x_forward, returns);
            f_backward = - log_likelihood_norm(x_backward, returns);

            grad[i] = (f_forward - f_backward) / (2 * epsilon);
        }
    }

    return - log_likelihood_norm(params, returns);
}

double constraint_norm(unsigned n, const double *x, double *grad, void *data)
{
    if (grad) {
        grad[0] = 0;
        grad[1] = 1;
        grad[2] = 1;
    }

    return x[1] + x[2] - 1; // alpha + beta - 1 < 0
}

Eigen::VectorXd estimate_norm(
    const Eigen::Ref<const Eigen::VectorXd>& returns
)
{
    int d = 3; // Dimension
    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, d); // L-BFGS

    double lower_bounds[d] = { std::numeric_limits<double>::epsilon(), 0.0, 0.0 };
    double upper_bounds[d] = { 1.0                                   , 1.0, 1.0 };
    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    nlopt_set_min_objective(opt, objective_norm, (void*) &returns);
    nlopt_add_inequality_constraint(opt, constraint_norm, NULL, 1e-8);
    nlopt_set_xtol_rel(opt, 1e-8);
    nlopt_set_ftol_rel(opt, 1e-8);
    nlopt_set_maxeval(opt, 250);

    // Same as rugarch in R, for GARCH(p, q):
    // initial guess = (var(returns) / 1000, 0.05 / p, 0.9 / q)
    double x[d] = { gsl_stats_variance(returns.data(), returns.stride(), returns.size()) / 1000, 0.05, 0.9 }; // Initial guess
    double f_min; // Smallest f(x)

    nlopt_result res = nlopt_optimize(opt, x, &f_min); // Optimize
    nlopt_destroy(opt);

    return Eigen::Map<Eigen::VectorXd>(x, d);
}

extern "C" SEXP R2Cpp_estimate_norm(SEXP R2Cpp_returns) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());

    Eigen::VectorXd estimate = estimate_norm(returns);
    return Rcpp::wrap(estimate);
}

double log_likelihood_std(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
)
{
    int n = returns.size();
    double nu = params(3);

    Eigen::VectorXd sigma = forecast_past(params, returns);

    double sum_log_likelihood = 0;
    double sigma2;

    for (int i = 0; i < n; ++i) {
        sigma2 = std::pow(sigma(i), 2);
        sum_log_likelihood += - 0.5 * std::log(sigma2) - ((nu + 1) / 2) * std::log(1 + std::pow(returns(i), 2) / (sigma2 * (nu - 2)));
    }
    sum_log_likelihood += n * std::log(std::tgamma((nu + 1) / 2) / (std::sqrt(nu - 2) * std::tgamma(nu / 2)));

    return sum_log_likelihood;
}

double objective_std(unsigned n, const double *x, double *grad, void *data)
{
    Eigen::VectorXd& returns = *(Eigen::VectorXd*) data;

    Eigen::Map<const Eigen::VectorXd> params(x, n);

    double epsilon = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);

    if (grad) {
        Eigen::VectorXd x_forward;
        Eigen::VectorXd x_backward;
        double f_forward;
        double f_backward;

        for (unsigned i = 0; i < n; ++i) {
            x_forward = params;
            x_backward = params;

            x_forward[i] += epsilon;
            x_backward[i] -= epsilon;

            f_forward = - log_likelihood_std(x_forward, returns);
            f_backward = - log_likelihood_std(x_backward, returns);

            grad[i] = (f_forward - f_backward) / (2 * epsilon);
        }
    }

    return - log_likelihood_std(params, returns);
}

double constraint_std(unsigned n, const double *x, double *grad, void *data) {
    if (grad) {
        grad[0] = 0;
        grad[1] = 1;
        grad[2] = 1;
        grad[3] = 0;
    }

    return x[1] + x[2] - 1; // alpha + beta - 1 <= 0
}

Eigen::VectorXd estimate_std(
    const Eigen::Ref<const Eigen::VectorXd>& returns
) 
{
    int d = 4; // Dimension
    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, d); // L-BFGS

    double lower_bounds[d] = { std::numeric_limits<double>::epsilon(), 0.0, 0.0,   2.1 };
    double upper_bounds[d] = { 1.0                                   , 1.0, 1.0, 100.0 };
    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    nlopt_set_min_objective(opt, objective_std, (void*) &returns);
    nlopt_add_inequality_constraint(opt, constraint_std, NULL, 1e-8);
    nlopt_set_xtol_rel(opt, 1e-8);
    nlopt_set_ftol_rel(opt, 1e-8);
    nlopt_set_maxeval(opt, 250);

    // Same as rugarch in R, for GARCH(p, q):
    // initial guess = (var(returns) / 1000, 0.05 / p, 0.9 / q, 4.0)
    double x[d] = { gsl_stats_variance(returns.data(), returns.stride(), returns.size()) / 1000, 0.05, 0.9, 4.0 }; // Initial guess
    double f_min; // Smallest f(x)

    nlopt_result res = nlopt_optimize(opt, x, &f_min); // Optimize
    nlopt_destroy(opt);

    return Eigen::Map<Eigen::VectorXd>(x, d);
}

extern "C" SEXP R2Cpp_estimate_std(SEXP R2Cpp_returns) {
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());

    Eigen::VectorXd estimate = estimate_std(returns);
    return Rcpp::wrap(estimate);
}
