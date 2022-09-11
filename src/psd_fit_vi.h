#ifndef INCLUDE_PSD_FIT_VI
#define INCLUDE_PSD_FIT_VI

#include <RcppEigen.h>

// Function declaration.
Eigen::MatrixXd rcpp_update_pp(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                               const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                               const Eigen::MatrixXd& ALPHA);

Rcpp::List rcpp_update_ff(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                          const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                          const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb);

double rcpp_marginal_likelihood_e1(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                                   const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF);

#endif
