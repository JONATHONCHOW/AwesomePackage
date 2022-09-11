// [[Rcpp::depends(RcppEigen)]]
#include "psd_fit_vi.h"
using namespace Rcpp;

// Update PP.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_pp(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                               const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                               const Eigen::MatrixXd& ALPHA)
{
  size_t I = G.rows();
  size_t J = G.cols();
  size_t K = ZP.cols();
  Eigen::MatrixXd PP = Eigen::MatrixXd::Zero(I, K);
  for (size_t i = 0; i < I; i++)
  {
    for (size_t j = 0; j < J; j++)
    {
      // Normalized for k.
      double normZa = 0, normZb = 0;
      for (size_t k = 0; k < K; k++)
      {
        normZa += ZaF(k, j) * ZP(i, k);
        normZb += ZbF(k, j) * ZP(i, k);
      }
      for (size_t k = 0; k < K; k++)
      {
        PP(i, k) += (((2 - G(i, j)) * ZbF(k, j) / normZb) + (G(i, j) * ZaF(k, j) / normZa)) * ZP(i, k);
      }
    }
    for (size_t k = 0; k < K; k++)
    {
      PP(i, k) = ALPHA(k, 0) + PP(i, k);
      PP(i, k) = ALPHA(k, 0) + PP(i, k);
    }
  }
  return PP;
}

// Update FFa and FFb.
// [[Rcpp::export]]
Rcpp::List rcpp_update_ff(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                          const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                          const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb)
{
  size_t I = G.rows();
  size_t J = G.cols();
  size_t K = ZP.cols();
  Eigen::MatrixXd FFa = Eigen::MatrixXd::Zero(K, J);
  Eigen::MatrixXd FFb = Eigen::MatrixXd::Zero(K, J);
  for (size_t j = 0; j < J; j++)
  {
    Eigen::MatrixXd tempa = Eigen::MatrixXd::Zero(K, 1);
    Eigen::MatrixXd tempb = Eigen::MatrixXd::Zero(K, 1);
    for (size_t i = 0; i < I; i++)
    {
      // Normalized for k.
      double normZa = 0, normZb = 0;
      for (size_t k = 0; k < K; k++)
      {
        normZa += ZaF(k, j) * ZP(i, k);
        normZb += ZbF(k, j) * ZP(i, k);
      }
      for (size_t k = 0; k < K; k++)
      {
        tempa(k, 0) += G(i, j) * ZP(i, k) / normZa;
        tempb(k, 0) += (2 - G(i, j)) * ZP(i, k) / normZb;
      }
    }
    for (size_t k = 0; k < K; k++)
    {
      FFa(k, j) = BETAa(k, j) + ZaF(k, j) * tempa(k, 0);
      FFb(k, j) = BETAb(k, j) + ZbF(k, j) * tempb(k, 0);
    }
  }
  return Rcpp::List::create(Rcpp::Named("FFa") = FFa,
                            Rcpp::Named("FFb") = FFb);
}

// Compute E1.
// [[Rcpp::export]]
double rcpp_marginal_likelihood_e1(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
                                   const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF)
{
  size_t I = G.rows();
  size_t J = G.cols();
  size_t K = ZP.cols();
  double E1 = 0;
  for (size_t j = 0; j < J; j++)
  {
    for (size_t i = 0; i < I; i++)
    {
      double tempa = 0, tempb = 0;
      if (G(i, j) == 0)
      {
        for (size_t k = 0; k < K; k++)
        {
          tempb += ZbF(k, j) * ZP(i, k);
        }
        tempa = tempb;
      }
      else if (G(i, j) == 1)
      {
        for (size_t k = 0; k < K; k++)
        {
          tempa += ZaF(k, j) * ZP(i, k);
          tempb += ZbF(k, j) * ZP(i, k);
        }
      }
      else if (G(i, j) == 2)
      {
        for (size_t k = 0; k < K; k++)
        {
          tempa += ZaF(k, j) * ZP(i, k);
        }
        tempb = tempa;
      }
      E1 += log(tempa) + log(tempb);
    }
  }
  return E1;
}
