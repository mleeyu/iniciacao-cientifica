#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::VectorXd forecast_past(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
)
{
    double omega = params(0);
    double alpha = params(1);
    double beta = params(2);

    int n = returns.size();
    Eigen::VectorXd sigma(n);

    sigma(0) = std::sqrt(omega / (1 - alpha - beta));

    for (int i = 1; i < n; ++i) {
        sigma(i) = std::sqrt(omega + alpha * std::pow(returns(i - 1), 2) + beta * std::pow(sigma(i - 1), 2));
    }

    return sigma;
}

extern "C" SEXP R2Cpp_forecast_past(SEXP R2Cpp_params, SEXP R2Cpp_returns) {
    Rcpp::NumericVector Rcpp_params = Rcpp::as<Rcpp::NumericVector>(R2Cpp_params);
    Eigen::Map<const Eigen::VectorXd> params(Rcpp_params.begin(), Rcpp_params.size());
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());

    Eigen::VectorXd sigma = forecast_past(params, returns);
    return Rcpp::wrap(sigma);
}

double forecast_one_step_ahead(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
)
{
    double omega = params(0);
    double alpha = params(1);
    double beta = params(2);

    int n = returns.size();
    double sigma = 0;
    double past_sigma = omega / (1 - alpha - beta);

    for (int i = 1; i < n; ++i) {
        sigma = std::sqrt(omega + alpha * std::pow(returns(i - 1), 2) + beta * std::pow(past_sigma, 2));
        past_sigma = sigma;
    }
  
    return std::sqrt(params(0) + params(1) * std::pow(returns(n - 1), 2) + params(2) * std::pow(sigma, 2));
}

extern "C" SEXP R2Cpp_forecast_one_step_ahead(SEXP R2Cpp_params, SEXP R2Cpp_returns) {
    Rcpp::NumericVector Rcpp_params = Rcpp::as<Rcpp::NumericVector>(R2Cpp_params);
    Eigen::Map<const Eigen::VectorXd> params(Rcpp_params.begin(), Rcpp_params.size());
    Rcpp::NumericVector Rcpp_returns = Rcpp::as<Rcpp::NumericVector>(R2Cpp_returns);
    Eigen::Map<const Eigen::VectorXd> returns(Rcpp_returns.begin(), Rcpp_returns.size());

    double sigma = forecast_one_step_ahead(params, returns);
    return Rcpp::wrap(sigma);
}
