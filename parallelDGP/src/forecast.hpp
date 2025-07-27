#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::VectorXd forecast_past(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
);

double forecast_one_step_ahead(
    const Eigen::Ref<const Eigen::VectorXd>& params,
    const Eigen::Ref<const Eigen::VectorXd>& returns
);
