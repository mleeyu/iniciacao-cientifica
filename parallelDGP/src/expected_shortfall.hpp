#pragma once

#include <Rcpp.h>

double xdnorm(double x, void* params);
double expected_shortfall_norm(double risk_level, double sigma, double value_at_risk);

double xdstd(double x, void* params);
double expected_shortfall_std(double risk_level, double sigma, double value_at_risk, double nu);
