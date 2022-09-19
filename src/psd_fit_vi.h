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

Eigen::MatrixXd rcpp_update_zp(const Eigen::MatrixXd& PP);

Rcpp::List rcpp_update_zf(const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb);

double rcpp_marginal_likelihood(const Eigen::MatrixXd& G,
                                const Eigen::MatrixXd& ZP, const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                                const Eigen::MatrixXd& PP, const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb,
                                const Eigen::MatrixXd& ALPHA, const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb);

double marginal_likelihood_e1(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                              const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF);

double marginal_likelihood_e2(const Eigen::MatrixXd& ZP,
                              const Eigen::MatrixXd& PP,
                              const Eigen::MatrixXd& ALPHA);

double marginal_likelihood_e3(const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                              const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb,
                              const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb);

#endif
