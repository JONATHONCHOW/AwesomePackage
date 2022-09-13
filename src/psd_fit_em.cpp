// [[Rcpp::depends(RcppEigen)]]
#include "psd_fit_em.h"
using namespace Rcpp;

// Update P.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_p_em(const Eigen::MatrixXd& G,
                                 const Eigen::MatrixXd& P,
                                 const Eigen::MatrixXd& F)
{
  size_t I = P.rows();
  size_t K = P.cols();
  size_t J = F.cols();
  Eigen::MatrixXd p = P;
  for (size_t i = 0; i < I; i++)
  {
    for (size_t k = 0; k < K; k++)
    {
      double temp1 = 0, temp2 = 0;
      for (size_t j = 0; j < J; j++)
      {
        temp1 += G(i, j) * fraction1(i, j, k, P, F);
        temp2 += (2 - G(i, j)) * fraction2(i, j, k, P, F);
      }
      p(i, k) = (temp1 + temp2) / 2 / J;
    }
  }
  return p;
}

// Update F.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_f_em(const Eigen::MatrixXd& G,
                                 const Eigen::MatrixXd& P,
                                 const Eigen::MatrixXd& F)
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
    }
  }
  return f;
}

// Compute loss function.
// [[Rcpp::export]]
double rcpp_psd_loss(const Eigen::MatrixXd& G,
                     const Eigen::MatrixXd& P,
                     const Eigen::MatrixXd& F)
{
  size_t I = P.rows();
  size_t J = F.cols();
  double loss = 0;
  for (size_t i = 0; i < I; i++)
  {
    for (size_t j = 0; j < J; j++)
    {
      loss += G(i, j) * log(binomial1(i, j, P, F)) + (2 - G(i, j)) * log(binomial2(i, j, P, F));
    }
  }
  loss = loss / (I * J);
  return loss;
}

// Compute the first binomial distribution parameter.
double binomial1(size_t i, size_t j,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F)
{
  double temp = 0;
  for (size_t r = 0; r < P.cols(); r++)
  {
    temp += P(i, r) * F(r, j);
  }
  return temp;
}

// Compute the second binomial distribution parameter.
double binomial2(size_t i, size_t j,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F)
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
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F)
{
  return P(i, k) * F(k, j) / binomial1(i, j, P, F);
}

// Compute the second fraction.
double fraction2(size_t i, size_t j, size_t k,
                 const Eigen::MatrixXd& P,
                 const Eigen::MatrixXd& F)
{
  return P(i, k) * (1 - F(k, j)) / binomial2(i, j, P, F);
}
