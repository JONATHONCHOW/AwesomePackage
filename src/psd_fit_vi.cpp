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

// Compute marginal likelihood.
// [[Rcpp::export]]
double rcpp_marginal_likelihood(const Eigen::MatrixXd& G,
                                const Eigen::MatrixXd& ZP, const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                                const Eigen::MatrixXd& PP, const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb,
                                const Eigen::MatrixXd& ALPHA, const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb)
{
  size_t I = G.rows();
  size_t J = G.cols();
  double E1 = 0, E2 = 0, E3 = 0, Etotal = 0;
  E1 = marginal_likelihood_e1(G, ZP, ZaF, ZbF);
  E2 = marginal_likelihood_e2(ZP, PP, ALPHA);
  E3 = marginal_likelihood_e3(ZaF, ZbF, FFa, FFb, BETAa, BETAb);
  Etotal = (E1 + E2 + E3) / (I * J);
  return Etotal;
}

// Compute E1.
double marginal_likelihood_e1(const Eigen::MatrixXd& G, const Eigen::MatrixXd& ZP,
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

// Compute E2.
double marginal_likelihood_e2(const Eigen::MatrixXd& ZP,
                              const Eigen::MatrixXd& PP,
                              const Eigen::MatrixXd& ALPHA)
{
  size_t I = ZP.rows();
  size_t K = ZP.cols();
  double E2 = 0;
  for (size_t i = 0; i < I; i++)
  {
    double temp1 = 0, temp2 = 0, temp3 = 0;
    for (size_t k = 0; k < K; k++)
    {
      temp1 += lgamma(PP(i, k)) - lgamma(ALPHA(k, 0)) - (PP(i, k) - ALPHA(k, 0)) * log(ZP(i, k));
      temp2 += PP(i, k);
      temp3 += ALPHA(k, 0);
    }
    E2 += temp1 - lgamma(temp2) + lgamma(temp3);
  }
  return E2;
}

// Compute E3.
double marginal_likelihood_e3(const Eigen::MatrixXd& ZaF, const Eigen::MatrixXd& ZbF,
                              const Eigen::MatrixXd& FFa, const Eigen::MatrixXd& FFb,
                              const Eigen::MatrixXd& BETAa, const Eigen::MatrixXd& BETAb)
{
  size_t K = ZaF.rows();
  size_t J = ZaF.cols();
  double E3 = 0;
  double temp1 = 0, temp2 = 0, temp3 = 0;
  for (size_t j = 0; j < J; j++)
  {
    for (size_t k = 0; k < K; k++)
    {
      temp1 = lgamma(FFa(k, j)) - lgamma(BETAa(k, j)) - (FFa(k, j) - BETAa(k, j)) * log(ZaF(k, j));
      temp2 = lgamma(FFb(k, j)) - lgamma(BETAb(k, j)) - (FFb(k, j) - BETAb(k, j)) * log(ZbF(k, j));
      temp3 = -lgamma(FFa(k, j) + FFb(k, j)) + lgamma(BETAa(k, j) + BETAb(k, j));
      E3 += temp1 + temp2 + temp3;
    }
  }
  return E3;
}
