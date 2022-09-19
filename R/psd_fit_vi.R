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
  FFa <- matrix(1, K, J) + 0.1 * rand(K, J)
  FFb <- 10 * matrix(1, K, J) + 0.1 * rand(K, J)
  ZP <- rcpp_update_zp(PP)
  ZF <- rcpp_update_zf(FFa, FFb)
  ZaF <- ZF$ZaF
  ZbF <- ZF$ZbF
  pre_L <- 0
  now_L <- rcpp_marginal_likelihood(G, ZP, ZaF, ZbF, PP, FFa, FFb, ALPHA, BETAa, BETAb)
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
    ZP <- rcpp_update_zp(PP)
    ZF <- rcpp_update_zf(FFa, FFb)
    ZaF <- ZF$ZaF
    ZbF <- ZF$ZbF
    if (iter %% 10 == 0)
    {
      pre_L <- now_L
      now_L <- rcpp_marginal_likelihood(G, ZP, ZaF, ZbF, PP, FFa, FFb, ALPHA, BETAa, BETAb)
      L_list <- append(L_list, now_L)
      if (! (abs(pre_L - now_L) > epsilon)) {break}
    }
    if (! (iter < maxiter)) {break}
  }
  # Posterior.
  P <- PP / rowSums(PP)
  F <- FFa / (FFa + FFb)
  return(list(P = P, F = F, Loss = L_list[-1], Iterations = iter))
}
