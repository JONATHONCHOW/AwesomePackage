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
#' @param maxiter The maximum number of iterations.
#'
#' @return A matrix represents the proportion of population distribution.
#'
#' @export
#'
#' @examples
#' G <- matrix(c(0,0,1,0,2, 1,0,0,1,0, 0,1,0,0,1), nrow = 3)
#' psd_fit_em(G, 2, 10)
psd_fit_em <- function (G, K, maxiter)
{
  init <- init_psd_fit(G, K)
  rcpp_psd_fit_em(init$P, init$F, G, maxiter)
}

# Init P and F.
init_psd_fit <- function (G, K)
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
