#' @title Use EM algorithm to fit PSD model
#'
#' @importFrom Rcpp evalCpp
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
#' \item{\code{Loss}}{A vector of length \code{Iterations + 1} represents the value of the loss function at each iteration.}
#' \item{\code{Iterations}}{An integer represents the number of iterations.}}
#'
#' @export
#'
#' @examples
#' G <- matrix(c(0,0,1, 0,2,1, 1,0,1, 0,1,0, 1,0,0), 3, 5)
#' psd_fit_em(G, 2, 1e-1, 10)
psd_fit_em <- function (G, K, epsilon = 1e-1, maxiter = 500)
{
  init <- init_psd_fit_em(G, K)
  rcpp_psd_fit_em(init$P, init$F, G, epsilon, maxiter)
}

# Init P and F.
init_psd_fit_em <- function (G, K)
{
  P <- rand(nrow(G), K)
  F <- rand(K, ncol(G))
  return (list(P = P, F = F))
}

# Generate a random matrix with elements between 0 and 1.
#
#' @importFrom stats runif
rand <- function (m, n, min = 0, max = 1)
{
  matrix(runif(m*n, min, max), m, n)
}
