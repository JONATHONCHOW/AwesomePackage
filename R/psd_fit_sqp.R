#' @title Use SQP algorithm to fit PSD model
#'
#' @importFrom Rcpp evalCpp
#' @importFrom progress progress_bar
#'
#' @description Fit PSD model with SQP algorithm, and use the loss
#'     function as a stopping criterion.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param K An integer 2 or greater giving the matrix rank.
#' @param epsilon Convergence criterion.
#' @param maxiter The maximum number of iterations.
#' @param initem A number of iterations when using EM algorithm for initialization.
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
#' psd_fit_sqp(G, 2, 1e-5, 10, 1)
psd_fit_sqp <- function (G, K, epsilon = 1e-5, maxiter = 50, initem = 100)
{
  I <- nrow(G)
  J <- ncol(G)
  # Add progress bar.
  pb <- progress_bar$new(
    format = '[:bar] :current/:total (:elapsed)',
    total = maxiter, clear = FALSE, width = 80
  )
  # Init.
  em <- psd_fit_em(G, K, epsilon, initem)
  P <- em$P
  F <- em$F
  pre_L <- 0
  now_L <- rcpp_psd_loss(G, P, F)
  L_list <- now_L
  # Loop.
  iter <- 0
  repeat
  {
    pb$tick()
    iter <- iter + 1
    P <- rcpp_update_p_sqp(G, P, F, 1e-9)
    F <- rcpp_update_f_sqp(G, P, F, 1e-9)
    if (iter %% 10 == 0)
    {
      pre_L <- now_L
      now_L <- rcpp_psd_loss(G, P, F)
      L_list <- append(L_list, now_L)
    }
    if (! (abs(pre_L - now_L) > epsilon && iter < maxiter) ) {break}
  }
  return(list(P=P, F=F, Loss = L_list, Iterations = iter))
}
