// [[Rcpp::depends(RcppEigen)]]
#include "psd_evaluate.h"
#include "psd_fit_em.h"
using namespace Rcpp;

// Compute the deviance residuals under the binomial model averaged over all entries as the prediction error.
// [[Rcpp::export]]
double rcpp_psd_error(const Eigen::MatrixXd& G,
                      const Eigen::MatrixXd& P,
                      const Eigen::MatrixXd& F)
{
  size_t I = G.rows();
  size_t J = G.cols();
  size_t K = P.cols();
  Eigen::MatrixXd GG = 2 * P * F;
  double res = 0;
  for (size_t i = 0; i < I; i++)
  {
    for (size_t j = 0; j < J; j++)
    {
      if (G(i, j) == 0)
      {
        res += 2 * (log(2) - log(2 - GG(i, j)));
      }
      if (G(i, j) == 1)
      {
        res += -log(GG(i, j)) - log((2 - GG(i, j)));
      }
      if (G(i, j) == 2)
      {
        res += 2 * (log(2) - log(GG(i, j)));
      }
    }
  }
  res = res / (I * J);
  return res;
}

// Update F in validation set.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_f_val(const Eigen::MatrixXd& G,
                                  const Eigen::MatrixXd& P,
                                  const Eigen::MatrixXd& F,
                                  const double& zero)
{
  size_t I = P.rows();
  size_t K = P.cols();
  size_t J = F.cols();
  Eigen::MatrixXd f = F;
  for (size_t k = 0; k < K; k++)
  {
    for (size_t j = 0; j < J; j++)
    {
      double temp1 = 0, temp2 = 0;
      for (size_t i = 0; i < I; i++)
      {
        temp1 += G(i, j) * fraction1(i, j, k, P, F);
        temp2 += (2 - G(i, j)) * fraction2(i, j, k, P, F);
      }
      f(k, j) = temp1 / (temp1 + temp2);
      if (f(k, j) < zero)
      {
        f(k, j) = zero;
      }
      if (f(k, j) > 1 - zero)
      {
        f(k, j) = 1 - zero;
      }
    }
  }
  return f;
}
