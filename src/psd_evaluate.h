#ifndef INCLUDE_PSD_EVALUATE
#define INCLUDE_PSD_EVALUATE

#include <RcppEigen.h>

// Function declaration.
double rcpp_psd_error(const Eigen::MatrixXd& G,
                      const Eigen::MatrixXd& P,
                      const Eigen::MatrixXd& F);

#endif
