// [[Rcpp::depends(RcppEigen)]]
#include "psd_fit_em.h"
using namespace Rcpp;

// Use EM algorithm to fit PSD model by Rcpp.
// [[Rcpp::export]]
Rcpp::List rcpp_psd_fit_em(const Eigen::MatrixXd& P,
                           const Eigen::MatrixXd& F,
                           const Eigen::MatrixXd& G,
                           const double& epsilon,
                           const size_t& maxiter)
{
  size_t I = P.rows();
  size_t K = P.cols();
  size_t J = F.cols();
  Eigen::MatrixXd p = P, f = F, g = G;
  std::vector<double> L_list;
  L_list.reserve(maxiter);
  double pre_L = 0, now_L = 0;
  now_L = psd_loss(p, f, g);
  L_list.push_back(now_L);
  size_t iter = 0;
  do
  {
    iter++;
    pre_L = now_L;
    now_L = 0;
    // Update P.
    Eigen::MatrixXd temp_p = p;
    for (size_t i = 0; i < I; i++)
    {
      for (size_t k = 0; k < K; k++)
      {
        double temp1 = 0, temp2 = 0;
        for (size_t j = 0; j < J; j++)
        {
          temp1 += g(i, j) * fraction1(i, j, k, p, f);
          temp2 += (2 - g(i, j)) * fraction2(i, j, k, p, f);
        }
        temp_p(i, k) = (temp1 + temp2) / 2 / J;
      }
    }
    // Update F.
    Eigen::MatrixXd temp_f = f;
    for (size_t k = 0; k < K; k++)
    {
      for (size_t j = 0; j < J; j++)
      {
        double temp1 = 0, temp2 = 0;
        for (size_t i = 0; i < I; i++)
        {
          temp1 += g(i, j) * fraction1(i, j, k, p, f);
          temp2 += (2 - g(i, j)) * fraction2(i, j, k, p, f);
        }
        temp_f(k, j) = temp1 / (temp1 + temp2);
      }
    }
    p = temp_p;
    f = temp_f;

    now_L = psd_loss(p, f, g);
    L_list.push_back(now_L);
  } while (fabs(pre_L - now_L) > epsilon && iter < maxiter);

  return Rcpp::List::create(Rcpp::Named("P") = p,
                            Rcpp::Named("F") = f,
                            Rcpp::Named("Loss") = L_list,
                            Rcpp::Named("Iterations") = iter);
}

// Compute the first sum.
double sum1(size_t i, size_t j,
            const Eigen::MatrixXd& P, const Eigen::MatrixXd& F)
{
  double temp = 0;
  for (size_t r = 0; r < P.cols(); r++)
  {
    temp += P(i, r) * F(r, j);
  }
  return temp;
}

// Compute the second sum.
double sum2(size_t i, size_t j,
            const Eigen::MatrixXd& P, const Eigen::MatrixXd& F)
{
  double temp = 0;
  for (size_t r = 0; r < P.cols(); r++)
  {
    temp += P(i, r) * (1 - F(r, j));
  }
  return temp;
}

// Compute the first fraction.
double fraction1(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P, const Eigen::MatrixXd& F)
{
  return P(i, k) * F(k, j) / sum1(i, j, P, F);
}

// Compute the second fraction.
double fraction2(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P, const Eigen::MatrixXd& F)
{
  return P(i, k) * (1 - F(k, j)) / sum2(i, j, P, F);
}

// Compute loss function.
double psd_loss(const Eigen::MatrixXd& P,
                const Eigen::MatrixXd& F,
                const Eigen::MatrixXd& G)
{
  double loss = 0;
  for (size_t i = 0; i < P.rows(); i++)
  {
    for (size_t j = 0; j < F.cols(); j++)
    {
      loss += G(i, j) * log(sum1(i, j, P, F)) + (2 - G(i, j)) * log(sum2(i, j, P, F));
    }
  }
  return loss;
}
