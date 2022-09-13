#ifndef INCLUDE_PSD_FIT_SQP
#define INCLUDE_PSD_FIT_SQP

#include <RcppEigen.h>

// Function declaration.
Eigen::MatrixXd rcpp_update_p_sqp(const Eigen::MatrixXd& G,
                                  const Eigen::MatrixXd& P,
                                  const Eigen::MatrixXd& F,
                                  const double& zero);

Eigen::MatrixXd rcpp_update_f_sqp(const Eigen::MatrixXd& G,
                                  const Eigen::MatrixXd& P,
                                  const Eigen::MatrixXd& F,
                                  const double& zero);

Eigen::MatrixXd eqqp_qr(const Eigen::MatrixXd& H, const Eigen::MatrixXd& f,
                        const Eigen::MatrixXd& Ae, const Eigen::MatrixXd& be);

Eigen::MatrixXd ineqqp_activeset(const Eigen::MatrixXd& X,
                                 const Eigen::MatrixXd& H, const Eigen::MatrixXd& f,
                                 const Eigen::MatrixXd& Ae, const Eigen::MatrixXd& be,
                                 const Eigen::MatrixXd& Ai, const Eigen::MatrixXd& bi,
                                 const double& epsilon = 1e-9, const size_t& maxiter = 50);

Eigen::MatrixXd partialP(size_t i, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G);

Eigen::MatrixXd partialF(size_t j, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G);

Eigen::MatrixXd hessianP(size_t i, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G);

Eigen::MatrixXd hessianF(size_t j, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G);

#endif
