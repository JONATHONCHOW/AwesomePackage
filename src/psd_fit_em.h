#ifndef INCLUDE_PSD_FIT_EM
#define INCLUDE_PSD_FIT_EM

#include <RcppEigen.h>

// Function declaration.
Rcpp::List rcpp_psd_fit_em(const Eigen::MatrixXd& P,
                           const Eigen::MatrixXd& F,
                           const Eigen::MatrixXd& G,
                           const double& epsilon,
                           const int& maxiter);

double sum1(size_t i, size_t j,
            const Eigen::MatrixXd& P, const Eigen::MatrixXd& F);

double sum2(size_t i, size_t j,
            const Eigen::MatrixXd& P, const Eigen::MatrixXd& F);

double fraction1(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P, const Eigen::MatrixXd& F);

double fraction2(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P, const Eigen::MatrixXd& F);

double psd_loss(const Eigen::MatrixXd& P,
                const Eigen::MatrixXd& F,
                const Eigen::MatrixXd& G);

#endif
