#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::MatrixXd simulate_norm(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    int n,
    int n_bis,
    double risk_level,
    int seed
);

Eigen::MatrixXd simulate_std(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    int n_oos,
    int n_ins,
    int n_bis,
    double risk_level,
    bool outlier
);
