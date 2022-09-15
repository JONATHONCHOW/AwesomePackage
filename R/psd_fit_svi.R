#' @title Use SVI algorithm to fit PSD model
#'
#' @importFrom Rcpp evalCpp
#' @importFrom progress progress_bar
#'
#' @description Fit PSD model with SVI algorithm, and use the loss
#'     function as a stopping criterion.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param K An integer 2 or greater giving the matrix rank.
#' @param epsilon Convergence criterion.
#' @param maxiter The maximum number of iterations.
#' @param subiter The number of iterations in the sampling section.
#' @param val_J Sample proportion of SNPs in validation set.
#' @param val_I Sample proportion of individuals in validation set.
#' @param val_iter The number of iterations between each validation set sampling.
#' @param tau A parameter of the descending direction of SVI algorithm.
#' @param kappa A parameter of the descending direction of SVI algorithm.
#'
#' @return A List with the following parameters:
#' \describe{
#' \item{\code{P}}{The population scale matrix of the individuals.}
#' \item{\code{Loss}}{A vector represents the value of the loss function which records once for 10 iterations.}
#' \item{\code{Iterations}}{An integer represents the number of iterations.}}
#'
#' @export
#'
#' @examples
#' # Refer to Articles in AwesomePackage.
psd_fit_svi <- function (G, K,
                         epsilon = 1e-5, maxiter = 5e+5, subiter = 100,
                         val_J = 5e-2, val_I = 1e-1, val_iter = 1e+4,
                         tau = 1, kappa = 0.5)
{
  ind_J <- sample(2, ncol(G), replace = TRUE, prob = c(1 - val_J, val_J))
  ind_I <- sample(2, nrow(G), replace = TRUE, prob = c(1 - val_I, val_I))
  G_train <- G[, ind_J == 1]
  G_val <- G[ind_I == 2, ind_J == 2]
  I <- nrow(G_train)
  J <- ncol(G_train)
  I_val <- nrow(G_val)
  J_val <- ncol(G_val)
  # Add progress bar.
  pb <- progress_bar$new(
    format = '[:bar] :current/:total (:elapsed)',
    total = maxiter, clear = FALSE, width = 80
  )
  # Init.
  ALPHA <- matrix(rep(1/K, K) ,K, 1)
  PP <- matrix(1, I, K) + 0.1 * rand(I, K)
  ZP <- rcpp_update_zp(PP)
  P <- matrix(0,I,K)
  for (i in 1:I)
  {
    P[i,] <- PP[i,] / sum(PP[i,])
  }
  # Validate init.
  ALPHA_val <- matrix(rep(1/K, K) ,K, 1)
  PP_val <- PP[ind_I == 2,]
  ZP_val <- ZP[ind_I == 2,]
  P_val <- ZP[ind_I == 2,]
  BETAa_val <- matrix(1, K, J_val)
  BETAb_val <- matrix(1, K, J_val)
  FFa_val <- matrix(1, K, J_val) + 0.1 * rand(K, J_val)
  FFb_val <- 10 * matrix(1, K, J_val) + 0.1 * rand(K, J_val)
  ZF_val <- rcpp_update_zf(FFa_val, FFb_val)
  ZaF_val <- ZF_val$ZaF
  ZbF_val <- ZF_val$ZbF
  F_val = FFa_val / (FFa_val + FFb_val)
  pre_L <- 0
  now_L <- rcpp_psd_loss(G_val, P_val, F_val)
  L_list <- now_L
  # Loop.
  iter <- 0
  repeat
  {
    pb$tick()
    iter <- iter + 1
    # Subloop.
    G_sample <- as.matrix(G_train[, sample(J, 1)])
    ffzf <- update_ffzf_svi(G_sample, ZP, subiter)
    FFa <- ffzf$FFa
    FFb <- ffzf$FFb
    ZaF <- ffzf$ZaF
    ZbF <- ffzf$ZbF
    rho <- (tau + iter)^(-kappa)
    PP <- rcpp_update_pp_svi(G_sample, PP, ZP, ZaF, ZbF, ALPHA, J, rho)
    ZP <- rcpp_update_zp(PP)
    # for (i in 1:I)
    # {
    #   P[i,] <- PP[i,] / sum(PP[i,])
    # }
    if (iter %% val_iter == 0)
    {
      for (i in 1:I)
      {
        P[i,] <- PP[i,] / sum(PP[i,])
      }
      ZP_val <- ZP[ind_I == 2,]
      P_val <- P[ind_I == 2,]
      ffzf_val <- update_ffzf_svi_val(G_val, ZP_val, P_val, epsilon)
      FFa <- ffzf_val$FFa
      FFb <- ffzf_val$FFb
      ZaF <- ffzf_val$ZaF
      ZbF <- ffzf_val$ZbF
      F_val = FFa_val / (FFa_val + FFb_val)
      pre_L <- now_L
      now_L <- rcpp_psd_loss(G_val, P_val, F_val)
      L_list <- append(L_list, now_L)
    }
    if (! (abs(pre_L - now_L) > epsilon * 1e-5 && iter < maxiter) ) {break}
  }
  for (i in 1:I)
  {
    P[i,] <- PP[i,] / sum(PP[i,])
  }
  return(list(P=P, Loss = L_list, Iterations = iter))
}

# Update FF and ZF in validation set.
update_ffzf_svi_val <- function(G, ZP, P, epsilon)
{
  K <- ncol(ZP)
  J <- ncol(G)
  BETAa <- matrix(1, K, J)
  BETAb <- matrix(1, K, J)
  FFa <- matrix(1, K, J) + 0.1 * rand(K, J)
  FFb <- 10 * matrix(1, K, J) + 0.1 * rand(K, J)
  ZF <- rcpp_update_zf(FFa, FFb)
  ZaF <- ZF$ZaF
  ZbF <- ZF$ZbF
  F = FFa / (FFa + FFb)
  pre_L <- 0
  now_L <- rcpp_psd_loss(G, P, F)
  iter <- 0
  repeat
  {
    iter <- iter + 1
    FF <- rcpp_update_ff(G, ZP, ZaF, ZbF, BETAa, BETAb)
    FFa <- FF$FFa
    FFb <- FF$FFb
    ZF <- rcpp_update_zf(FFa, FFb)
    ZaF <- ZF$ZaF
    ZbF <- ZF$ZbF
    if (iter %% 10 == 0)
    {
      F = FFa / (FFa + FFb)
      pre_L <- now_L
      now_L <- rcpp_psd_loss(G, P, F)
    }
    if (! (abs(pre_L - now_L) > epsilon)) {break}
  }
  return(list(FFa=FFa, FFb=FFb, ZaF=ZaF, ZbF=ZbF))
}

# Update FF and ZF.
update_ffzf_svi <- function(G, ZP, subiter)
{
  K <- ncol(ZP)
  J <- ncol(G)
  BETAa <- matrix(1, K, J)
  BETAb <- matrix(1, K, J)
  FFa <- matrix(1, K, J) + 0.1 * rand(K, J)
  FFb <- 10 * matrix(1, K, J) + 0.1 * rand(K, J)
  ZF <- rcpp_update_zf(FFa, FFb)
  ZaF <- ZF$ZaF
  ZbF <- ZF$ZbF
  iter <- 0
  repeat
  {
    iter <- iter + 1
    FF <- rcpp_update_ff(G, ZP, ZaF, ZbF, BETAa, BETAb)
    FFa <- FF$FFa
    FFb <- FF$FFb
    ZF <- rcpp_update_zf(FFa, FFb)
    ZaF <- ZF$ZaF
    ZbF <- ZF$ZbF
    if (! (iter < subiter)) {break}
  }
  return(list(FFa=FFa, FFb=FFb, ZaF=ZaF, ZbF=ZbF))
}
