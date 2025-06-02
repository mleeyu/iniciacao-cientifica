#include <Rcpp.h>
#include <RcppEigen.h>

#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

double xdnorm(double x, void* params) {
    double sigma = *((double*) params);
    return x * gsl_ran_gaussian_pdf(x, sigma);
}

double expected_shortfall_norm(double risk_level, double sigma, double value_at_risk) {
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
    gsl_function integrand;
    double abs_error;
    integrand.function = &xdnorm;
    integrand.params = &sigma;
    double expected_shortfall;
    gsl_integration_qagil(&integrand, value_at_risk, 1e-8, 1e-8, 1000, workspace, &expected_shortfall, &abs_error);
    gsl_integration_workspace_free(workspace);
    return expected_shortfall / risk_level;
}

extern "C" SEXP R2Cpp_expected_shortfall_norm(SEXP R2Cpp_risk_level, SEXP R2Cpp_sigma, SEXP R2Cpp_value_at_risk) {
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    double sigma = Rcpp::as<double>(R2Cpp_sigma);
    double value_at_risk = Rcpp::as<double>(R2Cpp_value_at_risk);

    double expected_shortfall = expected_shortfall_norm(risk_level, sigma, value_at_risk);
    return Rcpp::wrap(expected_shortfall);
}

double xdstd(double x, void* params) {
    double* parameters = (double*) params;
    double sigma = parameters[0];
    double nu = parameters[1];
    double s = std::sqrt((nu - 2) / nu);
    return x * gsl_ran_tdist_pdf(x / (s * sigma), nu) / (s * sigma);
}

double expected_shortfall_std(double risk_level, double sigma, double value_at_risk, double nu) {
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
    gsl_function integrand;
    double abs_error;
    integrand.function = &xdstd;
    double params[2] = { sigma, nu };
    integrand.params = params;
    double expected_shortfall;
    gsl_integration_qagil(&integrand, value_at_risk, 1e-8, 1e-8, 1000, workspace, &expected_shortfall, &abs_error);
    gsl_integration_workspace_free(workspace);
    return expected_shortfall / risk_level;
}

extern "C" SEXP R2Cpp_expected_shortfall_std(SEXP R2Cpp_risk_level, SEXP R2Cpp_sigma, SEXP R2Cpp_value_at_risk, SEXP R2Cpp_nu) {
    double risk_level = Rcpp::as<double>(R2Cpp_risk_level);
    double sigma = Rcpp::as<double>(R2Cpp_sigma);
    double value_at_risk = Rcpp::as<double>(R2Cpp_value_at_risk);
    double nu = Rcpp::as<double>(R2Cpp_nu);

    double expected_shortfall = expected_shortfall_std(risk_level, sigma, value_at_risk, nu);
    return Rcpp::wrap(expected_shortfall);
}
