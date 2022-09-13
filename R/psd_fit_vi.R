#' @title Use VI algorithm to fit PSD model
#'
#' @importFrom Rcpp evalCpp
#' @importFrom progress progress_bar
#'
#' @description Fit PSD model with VI algorithm, and use the loss
#'     function as a stopping criterion.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param K An integer 2 or greater giving the matrix rank.
#' @param epsilon Convergence criterion.
#' @param maxiter The maximum number of iterations.
#'
#' @return A List with the following parameters:
#' \describe{
#' \item{\code{P}}{The population scale matrix of the individuals.}
#' \item{\code{F}}{The gene scale matrix of the populations.}
#' \item{\code{Loss}}{A vector represents the value of the loss function which records once for 10 iterations.}
#' \item{\code{Iterations}}{An integer represents the number of iterations.}}
#'
#' @export
#'
#' @examples
#' G <- matrix(c(0,0,1, 0,2,1, 1,0,1, 0,1,0, 1,0,0), 3, 5)
#' psd_fit_vi(G, 2, 1e-5, 10)
psd_fit_vi <- function (G, K, epsilon = 1e-5, maxiter = 500)
{
  I <- nrow(G)
  J <- ncol(G)
  # Add progress bar.
  pb <- progress_bar$new(
    format = '[:bar] :current/:total (:elapsed)',
    total = maxiter, clear = FALSE, width = 80
  )
  # Init.
  ALPHA <- matrix(rep(1/K, K) ,K, 1)
  BETAa <- matrix(1, K, J)
  BETAb <- matrix(1, K, J)
  PP <- matrix(1, I, K) + 0.1 * rand(I, K)
  ZP <- update_zp(PP)
  FFa <- matrix(1, K, J) + 0.1 * rand(K, J)
  FFb <- 10 * matrix(1, K, J) + 0.1 * rand(K, J)
  ZF <- update_zf(FFa, FFb)
  ZaF <- ZF$ZaF
  ZbF <- ZF$ZbF
  pre_L <- 0
  now_L <- marginal_likelihood(G, ZP, ZaF, ZbF, PP, FFa, FFb, ALPHA, BETAa, BETAb)
  L_list <- now_L
  # Loop.
  iter <- 0
  repeat
  {
    pb$tick()
    iter <- iter + 1
    PP <- rcpp_update_pp(G, ZP, ZaF, ZbF, ALPHA)
    FF <- rcpp_update_ff(G, ZP, ZaF, ZbF, BETAa, BETAb)
    FFa <- FF$FFa
    FFb <- FF$FFb
    ZP <- update_zp(PP)
    ZF <- update_zf(FFa, FFb)
    ZaF <- ZF$ZaF
    ZbF <- ZF$ZbF
    if (iter %% 10 == 0)
    {
      pre_L <- now_L
      now_L <- marginal_likelihood(G, ZP, ZaF, ZbF, PP, FFa, FFb, ALPHA, BETAa, BETAb)
      L_list <- append(L_list, now_L)
    }
    if (! (abs(pre_L - now_L) > epsilon && iter < maxiter) ) {break}
  }
  # Posterior.
  P <- matrix(0,I,K)
  for (i in 1:I)
  {
    P[i,] <- PP[i,] / sum(PP[i,])
  }
  F = FFa / (FFa + FFb)
  return(list(P=P, F=F, Loss = L_list, Iterations = iter))
}

# Update ZP.
update_zp <- function (PP)
{
  I <- nrow(PP)
  K <- ncol(PP)
  ZP <- matrix(nrow = I, ncol = K)
  for (i in 1:I)
  {
    temp <- digamma(sum(PP[i,]))
    ZP[i,] <- exp(digamma(PP[i,]) - temp)
  }
  return(ZP)
}

# Update ZaF and ZbF.
update_zf <- function (FFa, FFb)
{
  ZaF <- exp(digamma(FFa) - digamma(FFa + FFb))
  ZbF <- exp(digamma(FFb) - digamma(FFa + FFb))
  return(list(ZaF=ZaF, ZbF=ZbF))
}

# Compute marginal likelihood.
marginal_likelihood <- function (G, ZP, ZaF, ZbF,
                                 PP, FFa, FFb,
                                 ALPHA, BETAa, BETAb)
{
  I <- nrow(G)
  J <- ncol(G)
  E1 <- rcpp_marginal_likelihood_e1(G, ZP, ZaF, ZbF)
  E2 <- 0
  for (i in 1:I)
  {
    E2 <- E2 + sum(lgamma(PP[i,]) - lgamma(ALPHA) - (PP[i,]-ALPHA)*log(ZP[i,])) - lgamma(sum(PP[i,])) + lgamma(sum(ALPHA))
  }
  E3 <- sum(lgamma(FFa) - lgamma(BETAa) - (FFa-BETAa)*log(ZaF)
            + lgamma(FFb) - lgamma(BETAb) - (FFb-BETAb)*log(ZbF)
            - lgamma(FFa+FFb) + lgamma(BETAa+BETAb))
  Etotal <- (E1 + E2 + E3) / (I * J)
  return(Etotal)
}
