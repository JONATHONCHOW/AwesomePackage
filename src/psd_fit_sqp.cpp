// [[Rcpp::depends(RcppEigen)]]
#include "psd_fit_em.h"
#include "psd_fit_sqp.h"
using namespace Rcpp;

// Use SQP algorithm to fit PSD model by Rcpp.
// [[Rcpp::export]]
Rcpp::List rcpp_psd_fit_sqp(const Eigen::MatrixXd& P,
                            const Eigen::MatrixXd& F,
                            const Eigen::MatrixXd& G,
                            const double& epsilon,
                            const size_t& maxiter,
                            const double& zero)
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
    for (size_t i = 0; i < I; i++)
    {
      Eigen::MatrixXd deltaP(K, 1);
      deltaP = Eigen::MatrixXd::Zero(K, 1);
      Eigen::MatrixXd parP(K, 1);
      parP = partialP(i, p, f, g);
      Eigen::MatrixXd hesP(K, K);
      hesP = hessianP(i, p, f, g);
      Eigen::MatrixXd Ae(1, K);
      Ae = Eigen::MatrixXd::Ones(1, K);
      Eigen::MatrixXd be(1, 1);
      be << 0;
      Eigen::MatrixXd Ai(K, K);
      Ai << Eigen::MatrixXd::Identity(K, K);
      Eigen::MatrixXd bi(K, 1);
      bi << -p.row(i).transpose();
      double sc = hesP.norm();
      deltaP = ineqqp_activeset(deltaP, hesP / sc, parP / sc, Ae, be, Ai, bi);
      for (size_t k = 0; k < K; k++)
      {
        p(i, k) += deltaP(k, 0);
      }
    }
    for (size_t i = 0; i < I; i++)
    {
      for (size_t k = 0; k < K; k++)
      {
        if (p(i, k) < zero)
        {
          p(i, k) = zero;
        }
        if (p(i, k) > 1 - zero)
        {
          p(i, k) = 1 - zero;
        }
      }
    }
    // Update F.
    for (size_t j = 0; j < J; j++)
    {
      Eigen::MatrixXd deltaF(K, 1);
      deltaF = Eigen::MatrixXd::Zero(K, 1);
      Eigen::MatrixXd parF(K, 1);
      parF = partialF(j, p, f, g);
      Eigen::MatrixXd hesF(K, K);
      hesF = hessianF(j, p, f, g);
      Eigen::MatrixXd Ae(0, 0);
      Eigen::MatrixXd be(0, 0);
      Eigen::MatrixXd Ai(2 * K, K);
      Ai << Eigen::MatrixXd::Identity(K, K),
            -Eigen::MatrixXd::Identity(K, K);
      Eigen::MatrixXd bi(2 * K, 1);
      bi << -f.col(j),
            f.col(j) - Eigen::MatrixXd::Ones(K, 1);
      double sc = hesF.norm();
      deltaF = ineqqp_activeset(deltaF, hesF / sc, parF / sc, Ae, be, Ai, bi);
      for (size_t k = 0; k < K; k++)
      {
        f(k, j) += deltaF(k, 0);
      }
    }
    for (size_t j = 0; j < J; j++)
    {
      for (size_t k = 0; k < K; k++)
      {
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

    now_L = psd_loss(p, f, g);
    L_list.push_back(now_L);
  } while (fabs(pre_L - now_L) > epsilon && iter <= maxiter);

  return Rcpp::List::create(Rcpp::Named("P") = p,
                            Rcpp::Named("F") = f,
                            Rcpp::Named("Loss") = L_list,
                            Rcpp::Named("Iterations") = iter + 1);
}

// Solve equality-constrained quadratic programs using Lagrange duel and QR decomposition.
Eigen::MatrixXd eqqp_qr(const Eigen::MatrixXd& H, const Eigen::MatrixXd& f,
                        const Eigen::MatrixXd& Ae, const Eigen::MatrixXd& be)
{
  Eigen::MatrixXd result;
  size_t n = H.rows();
  size_t ne = Ae.rows();

  // Compute KKT matrix.
  Eigen::MatrixXd K(n + ne, n + ne);
  Eigen::MatrixXd q(n + ne, 1);
  K << H, -Ae.transpose(),
       Ae, Eigen::MatrixXd::Zero(ne, ne);
  q << -f,
       be;

  // Compute x.
  result = K.colPivHouseholderQr().solve(q);
  return result;
}

// Solve inequality-constrained quadratic programs using active-set methods.
Eigen::MatrixXd ineqqp_activeset(const Eigen::MatrixXd& X,
                                 const Eigen::MatrixXd& H, const Eigen::MatrixXd& f,
                                 const Eigen::MatrixXd& Ae, const Eigen::MatrixXd& be,
                                 const Eigen::MatrixXd& Ai, const Eigen::MatrixXd& bi,
                                 const double& epsilon, const size_t& maxiter)
{
  Eigen::MatrixXd x = X;
  size_t n = H.rows();
  size_t ne = Ae.rows();
  size_t ni = Ai.rows();

  // Init active set.
  Eigen::MatrixXd Sk = Eigen::MatrixXd::Zero(ni, 1);
  std::vector<double> asIndex;
  Eigen::MatrixXd Aee(1, 1);
  Eigen::MatrixXd Aix(ni, 1);
  Aix = Ai * x;
  for (size_t i = 0; i < ni; i++)
  {
    if (Aix(i, 0) <= bi(i, 0) + epsilon)
    {
      Sk(i, 0) = 1;
    }
  }

  // Main loop.
  size_t iter = 0;
  while (iter < maxiter)
  {
    size_t asSize = 0;
    for (size_t i = 0; i < ni; i++)
    {
      if (Sk(i, 0) == 1)
      {
        asSize += 1;
      }
    }
    if (ne != 0)
    {
      Aee.resize(ne + asSize, n);
      Aee.block(0, 0, ne, n) = Ae;
      size_t aeeIndex = 0;
      for (size_t i = 0; i < ni; i++)
      {
        if (Sk(i, 0) == 1)
        {
          Aee.block(ne + aeeIndex, 0, 1, n) = Ai.row(i);
          asIndex.push_back(i);
          aeeIndex += 1;
        }
      }
    }
    else {
      Aee.resize(asSize, n);
      size_t aeeIndex = 0;
      for (size_t i = 0; i < ni; i++)
      {
        if (Sk(i, 0) == 1)
        {
          Aee.block(aeeIndex, 0, 1, n) = Ai.row(i);
          asIndex.push_back(i);
          aeeIndex += 1;
        }
      }
    }

    // Solve subproblem.
    size_t nee = Aee.rows();
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(n, 1);
    Eigen::MatrixXd lambda = Eigen::MatrixXd::Zero(nee, 1);
    Eigen::MatrixXd gk;
    gk = H * x + f;
    Eigen::MatrixXd bee = Eigen::MatrixXd::Zero(nee, 1);
    Eigen::MatrixXd temp(n + nee, 1);
    temp = eqqp_qr(H, gk, Aee, bee);
    d = temp.block(0, 0, n, 1);
    lambda = temp.block(n, 0, nee, 1);
    // Calculate the length of d.
    double dNorm = 0;
    for (size_t i = 0; i < n; i++)
    {
      dNorm += d(i) * d(i);
    }
    if (dNorm < epsilon)
    {
      // Calculate the stop conditions.
      double minlambda = 0;
      size_t minlambdaIndex;
      for (size_t i = ne; i < nee; i++)
      {
        if (lambda(i, 0) < minlambda)
        {
          minlambda = lambda(i, 0);
          minlambdaIndex = i;
        }
      }
      if (minlambda >= 0)
      {
        break;
      }
      // Eliminate constraints from the active set.
      else
      {
        size_t removed_cons = asIndex[minlambdaIndex - ne];
        Sk(removed_cons, 0) = 0;
      }
    }
    // Caculate step length.
    else
    {
      Eigen::MatrixXd ad(ni, 1);
      ad = Ai * d;
      Eigen::MatrixXd ax(ni, 1);
      ax = Ai * x;
      double minalpha = 1;
      size_t minalphaIndex;
      for (size_t i = 0; i < ni; i++)
      {
        if (Sk(i) == 0 && ad(i, 0) < 0)
        {
          double alpha = (bi(i, 0) - ax(i, 0)) / ad(i, 0);
          if (alpha < minalpha)
          {
            minalpha = alpha;
            minalphaIndex = i;
          }
        }
      }
      if (minalpha >= 1)
      {
        x = x + d;
      }
      else
      {
        x = x + minalpha * d;
        // Add the constraint corresponding to the smallest alpha to the active set.
        Sk(minalphaIndex) = 1;
      }
    }

    iter += 1;
    asIndex.clear();
  }

  return x;
}

Eigen::MatrixXd partialP(size_t i, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G)
{
  size_t K = P.cols();
  size_t J = F.cols();
  Eigen::MatrixXd parP(K, 1);
  for (size_t k = 0; k < K; k++)
  {
    double temp1 = 0, temp2 = 0;
    for (size_t j = 0; j < J; j++)
    {
      temp1 += G(i, j) * F(k, j) / sum1(i, j, P, F);
      temp2 += (2 - G(i, j)) * (1 - F(k, j)) / sum2(i, j, P, F);
    }
    parP(k, 0) = temp1 + temp2;
  }
  return parP;
}

Eigen::MatrixXd partialF(size_t j, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G)
{
  size_t K = P.cols();
  size_t I = P.rows();
  Eigen::MatrixXd parF(K, 1);
  for (size_t k = 0; k < K; k++)
  {
    double temp1 = 0, temp2 = 0;
    for (size_t i = 0; i < I; i++)
    {
      temp1 += G(i, j) * P(i, k) / sum1(i, j, P, F);
      temp2 += -(2 - G(i, j)) * P(i, k) / sum2(i, j, P, F);
    }
    parF(k, 0) = temp1 + temp2;
  }
  return parF;
}

Eigen::MatrixXd hessianP(size_t i, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G)
{
  size_t K = P.cols();
  size_t J = F.cols();
  Eigen::MatrixXd hesP(K, K);
  for (size_t k = 0; k < K; k++)
  {
    for (size_t l = 0; l < k + 1; l++)
    {
      double temp1 = 0, temp2 = 0;
      for (size_t j = 0; j < J; j++)
      {
        temp1 += -G(i, j) * F(k, j) * F(l, j) / sum1(i, j, P, F) / sum1(i, j, P, F);
        temp2 += -(2 - G(i, j)) * (1 - F(k, j)) * (1 - F(l, j)) / sum2(i, j, P, F) / sum2(i, j, P, F);
      }
      hesP(k, l) = temp1 + temp2;
      hesP(l, k) = temp1 + temp2;
    }
  }
  return hesP;
}

Eigen::MatrixXd hessianF(size_t j, const Eigen::MatrixXd& P,
                         const Eigen::MatrixXd& F, const Eigen::MatrixXd& G)
{
  size_t K = P.cols();
  size_t I = P.rows();
  Eigen::MatrixXd hesF(K, K);
  for (size_t k = 0; k < K; k++)
  {
    for (size_t l = 0; l < k + 1; l++)
    {
      double temp1 = 0, temp2 = 0;
      for (size_t i = 0; i < I; i++)
      {
        temp1 += -G(i, j) * P(i, k) * P(i, l) / sum1(i, j, P, F) / sum1(i, j, P, F);
        temp2 += -(2 - G(i, j)) * P(i, k) * P(i, l) / sum2(i, j, P, F) / sum2(i, j, P, F);
      }
      hesF(k, l) = temp1 + temp2;
      hesF(l, k) = temp1 + temp2;
    }
  }
  return hesF;
}
