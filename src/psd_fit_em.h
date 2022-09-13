#ifndef INCLUDE_PSD_FIT_EM
#define INCLUDE_PSD_FIT_EM

#include <RcppEigen.h>

// Function declaration.
Eigen::MatrixXd rcpp_update_p_em(const Eigen::MatrixXd& G,
                                 const Eigen::MatrixXd& P,
                                 const Eigen::MatrixXd& F);

Eigen::MatrixXd rcpp_update_f_em(const Eigen::MatrixXd& G,
                                 const Eigen::MatrixXd& P,
                                 const Eigen::MatrixXd& F);

double rcpp_psd_loss(const Eigen::MatrixXd& G,
                     const Eigen::MatrixXd& P,
                     const Eigen::MatrixXd& F);

double binomial1(size_t i, size_t j,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F);

double binomial2(size_t i, size_t j,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F);

double fraction1(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F);

double fraction2(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F);

#endif
