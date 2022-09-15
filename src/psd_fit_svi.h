#ifndef INCLUDE_PSD_FIT_SVI
#define INCLUDE_PSD_FIT_SVI

#include <RcppEigen.h>

// Function declaration.
Eigen::MatrixXd rcpp_update_pp_svi(const Eigen::MatrixXd& G,
                                   const Eigen::MatrixXd& PP, const Eigen::MatrixXd& ZP,
                                   const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                                   const Eigen::MatrixXd& ALPHA,
                                   const size_t& J, const double& rho);

Eigen::MatrixXd rcpp_update_zp(const Eigen::MatrixXd& PP);

Rcpp::List rcpp_update_zf(const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb);

#endif
