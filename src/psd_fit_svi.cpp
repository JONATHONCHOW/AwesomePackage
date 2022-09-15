// [[Rcpp::depends(RcppEigen)]]
#include "psd_fit_svi.h"
using namespace Rcpp;

// Update PP.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_pp_svi(const Eigen::MatrixXd& G,
                                   const Eigen::MatrixXd& PP, const Eigen::MatrixXd& ZP,
                                   const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                                   const Eigen::MatrixXd& ALPHA,
                                   const size_t& J, const double& rho)
{
  size_t I = G.rows();
  size_t K = ZP.cols();
  Eigen::MatrixXd pp = Eigen::MatrixXd::Zero(I, K);
  for (size_t i = 0; i < I; i++)
  {
    // Normalized for k.
    double normZa = 0, normZb = 0;
    for (size_t k = 0; k < K; k++)
    {
      normZa += ZaF(k, 0) * ZP(i, k);
      normZb += ZbF(k, 0) * ZP(i, k);
    }
    for (size_t k = 0; k < K; k++)
    {
      pp(i, k) = (1 - rho) * PP(i, k) + rho * (ALPHA(k, 0) + J * (((2 - G(i, 0)) * ZbF(k, 0) / normZb) + (G(i, 0) * ZaF(k, 0) / normZa)) * ZP(i, k));
    }
  }
  return pp;
}

// Update ZP.
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_zp(const Eigen::MatrixXd& PP)
{
  size_t I = PP.rows();
  size_t K = PP.cols();
  Eigen::MatrixXd ZP = Eigen::MatrixXd::Zero(I, K);
  for (size_t i = 0; i < I; i++)
  {
    double temp = 0;
    for (size_t k = 0; k < K; k++)
    {
      temp += PP(i, k);
    }
    for (size_t k = 0; k < K; k++)
    {
      ZP(i, k) = exp(R::digamma(PP(i, k)) - R::digamma(temp));
    }
  }
  return ZP;
}

// Update ZaF and ZbF.
// [[Rcpp::export]]
Rcpp::List rcpp_update_zf(const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb)
{
  size_t J = FFa.cols();
  size_t K = FFa.rows();
  Eigen::MatrixXd ZaF = Eigen::MatrixXd::Zero(K, J);
  Eigen::MatrixXd ZbF = Eigen::MatrixXd::Zero(K, J);
  for (size_t j = 0; j < J; j++)
  {
    for (size_t k = 0; k < K; k++)
    {
      ZaF(k, j) = exp(R::digamma(FFa(k, j)) - R::digamma(FFa(k, j) + FFb(k, j)));
      ZbF(k, j) = exp(R::digamma(FFb(k, j)) - R::digamma(FFa(k, j) + FFb(k, j)));
    }
  }
  return Rcpp::List::create(Rcpp::Named("ZaF") = ZaF,
                            Rcpp::Named("ZbF") = ZbF);
}
