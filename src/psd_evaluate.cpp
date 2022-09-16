// [[Rcpp::depends(RcppEigen)]]
#include "psd_evaluate.h"
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
