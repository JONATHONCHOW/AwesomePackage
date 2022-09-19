#' @title Use EM algorithm to fit PSD model
#'
#' @importFrom Rcpp evalCpp
#' @importFrom progress progress_bar
#'
#' @description Fit PSD model with EM algorithm, and use the loss
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
#' psd_fit_em(G, 2, 1e-5, 10)
psd_fit_em <- function (G, K, epsilon = 1e-5, maxiter = 500)
{
  I <- nrow(G)
  J <- ncol(G)
  # Add progress bar.
  pb <- progress_bar$new(
    format = '[:bar] :current/:total (:elapsed)',
    total = maxiter, clear = FALSE, width = 80
  )
  # Init.
  P <- rand(I, K)
  F <- rand(K, J)
  pre_L <- 0
  now_L <- rcpp_psd_loss(G, P, F)
  L_list <- now_L
  # Loop.
  iter <- 0
  repeat
  {
    pb$tick()
    iter <- iter + 1
    P <- rcpp_update_p_em(G, P, F)
    F <- rcpp_update_f_em(G, P, F)
    if (iter %% 10 == 0)
    {
      pre_L <- now_L
      now_L <- rcpp_psd_loss(G, P, F)
      L_list <- append(L_list, now_L)
      if (! (abs(pre_L - now_L) > epsilon)) {break}
    }
    if (! (iter < maxiter)) {break}
  }
  return(list(P = P, F = F, Loss = L_list[-1], Iterations = iter))
}

# Generate a random matrix with elements between 0 and 1.
#
#' @importFrom stats runif
rand <- function (m, n, min = 0, max = 1)
{
  matrix(runif(m*n, min, max), m, n)
}
