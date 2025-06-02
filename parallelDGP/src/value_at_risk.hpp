#pragma once

#include <Rcpp.h>

double value_at_risk_norm(double risk_level, double sigma);

double value_at_risk_std(double risk_level, double sigma, double nu);
