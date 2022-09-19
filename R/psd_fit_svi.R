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
#' @param val_iter The number of iterations between each validation set sampling.
#' @param maxdrop The maximum number of consecutive decreases in the loss function. Beyond this value the loop will stop.
#' @param maxiter.sample The maximum number of iterations in the sampling section.
#' @param maxiter.val The maximum number of iterations in the validation section.
#' @param val_J Sample proportion of SNPs in validation set.
#' @param val_I Sample proportion of individuals in validation set.
#' @param tau A parameter of the descending direction of SVI algorithm.
#' @param kappa A parameter of the descending direction of SVI algorithm.
#'
#' @return A List with the following parameters:
#' \describe{
#' \item{\code{P}}{The population scale matrix of the individuals.}
#' \item{\code{Loss}}{A vector represents the value of the loss function which records once for 10 iterations.}
#' \item{\code{MaxLoss}}{Maximum loss function value. Unlike other algorithms, we observe the loss function on the validation set. Therefore, monotonicity is not guaranteed, that is, the maximum value does not necessarily occur at the end, so the maximum value needs to be recorded.}
#' \item{\code{Iterations}}{An integer represents the number of iterations.}}
#'
#' @export
#'
#' @examples
#' # Refer to Articles in AwesomePackage.
psd_fit_svi <- function (G, K,
                         epsilon = 1e-5, maxiter = 5e+5, val_iter = 1e+4, maxdrop = 3,
                         maxiter.sample = 100, maxiter.val = 2000,
                         val_J = 5e-2, val_I = 1e-1,
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
  pre_L <- 0
  now_L <- -1
  L_list <- now_L
  max_L <- now_L
  drop.index <- 0
  # Loop.
  iter <- 0
  repeat
  {
    pb$tick()
    iter <- iter + 1
    # Sample.
    G_sample <- as.matrix(G_train[, sample(J, 1)])
    ffzf_sample <- update_ffzf_svi(G_sample, ZP, maxiter.sample)
    ZaF_sample <- ffzf_sample$ZaF
    ZbF_sample <- ffzf_sample$ZbF
    rho <- (tau + iter)^(-kappa)
    PP <- rcpp_update_pp_svi(G_sample, PP, ZP, ZaF_sample, ZbF_sample, ALPHA, J, rho)
    ZP <- rcpp_update_zp(PP)
    # Validation.
    if (iter %% val_iter == 0)
    {
      P <- PP / rowSums(PP)
      ZP_val <- ZP[ind_I == 2,]
      P_val <- P[ind_I == 2,]
      ffzf_val <- update_ffzf_svi(G_val, ZP_val, maxiter.val)
      FFa_val <- ffzf_val$FFa
      FFb_val <- ffzf_val$FFb
      F_val <- FFa_val / (FFa_val + FFb_val)
      pre_L <- now_L
      now_L <- rcpp_psd_loss(G_val, P_val, F_val)
      L_list <- append(L_list, now_L)
      if (now_L > pre_L && abs(pre_L - now_L) < epsilon) {break}
      else if (now_L < pre_L) {drop.index <- drop.index + 1}
      else if (now_L > pre_L) {drop.index <- 0}
      if (now_L > max_L) {max_L <- now_L}
      if (drop.index > maxdrop) {break}
    }
    if (! (iter < maxiter) ) {break}
  }
  P <- PP / rowSums(PP)
  return(list(P = P, Loss = L_list[-1], MaxLoss = max_L, Iterations = iter))
}

# Update FF and ZF in sample set or validation set.
update_ffzf_svi <- function(G, ZP, maxiter)
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
    if (! (iter < maxiter)) {break}
  }
  return(list(FFa=FFa, FFb=FFb, ZaF=ZaF, ZbF=ZbF))
}
