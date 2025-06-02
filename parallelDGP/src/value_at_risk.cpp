#include <Rcpp.h>
#include <RcppEigen.h>

#include <RcppGSL.h>
#include <gsl/gsl_cdf.h>

double value_at_risk_norm(double risk_level, double sigma) {
    return sigma * gsl_cdf_gaussian_Pinv(risk_level, 1);
}

extern "C" SEXP R2Cpp_value_at_risk_norm(SEXP R2Cpp_risk_level, SEXP R2Cpp_sigma) {
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    double sigma = Rcpp::as<double>(R2Cpp_sigma);

    double value_at_risk = value_at_risk_norm(risk_level, sigma);
    return Rcpp::wrap(value_at_risk);
}

double value_at_risk_std(double risk_level, double sigma, double nu) {
    double s = std::sqrt((nu - 2) / nu);
    return sigma * gsl_cdf_tdist_Pinv(risk_level, nu) * s;
}

extern "C" SEXP R2Cpp_value_at_risk_std(SEXP R2Cpp_risk_level, SEXP R2Cpp_sigma, SEXP R2Cpp_nu) {
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    double sigma = Rcpp::as<double>(R2Cpp_sigma);
    double nu = Rcpp::as<double>(R2Cpp_nu);

    double value_at_risk = value_at_risk_std(risk_level, sigma, nu);
    return Rcpp::wrap(value_at_risk);
}
