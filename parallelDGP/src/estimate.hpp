#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

double log_likelihood_norm(
    const Eigen::Ref<const Eigen::VectorXd>& parameters,
    const Eigen::Ref<const Eigen::VectorXd>& returns
);
double objective_norm(unsigned n, const double *x, double *grad, void *data);
double constraint_norm(unsigned n, const double *x, double *grad, void *data);
Eigen::VectorXd estimate_norm(
    const Eigen::Ref<const Eigen::VectorXd>& returns
);

double log_likelihood_std(
    const Eigen::Ref<const Eigen::VectorXd>& parameters,
    const Eigen::Ref<const Eigen::VectorXd>& returns
);
double objective_std(unsigned n, const double *x, double *grad, void *data);
double constraint_std(unsigned n, const double *x, double *grad, void *data);
Eigen::VectorXd estimate_std(
    const Eigen::Ref<const Eigen::VectorXd>& returns
);
