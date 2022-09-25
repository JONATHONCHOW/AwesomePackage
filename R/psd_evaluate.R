#' @title Compute the prediction error
#'
#' @importFrom Rcpp evalCpp
#'
#' @description Compute the deviance residuals under the binomial
#'     model averaged over all entries as the prediction error.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param result Output of \code{psd_fit_em}, \code{psd_fit_sqp}, or \code{psd_fit_vi}.
#'
#' @return A value indicates the accuracy of the prediction.
#'
#' @export
#'
#' @examples
#' G <- matrix(c(0,0,1, 0,2,1, 1,0,1, 0,1,0, 1,0,0), 3, 5)
#' result <- psd_fit_em(G, 2, 1e-5, 10)
#' psd_error(G, result)
psd_error <- function (G, result)
{
  rcpp_psd_error(G, result$P, result$F)
}

#' @title Compute the loglikelihood
#'
#' @importFrom Rcpp evalCpp
#'
#' @description Compute the maximum loglikelihood function of the PSD model.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param result Output of \code{psd_fit_em}, \code{psd_fit_sqp}, or \code{psd_fit_vi}.
#'
#' @return A value indicates the accuracy of the prediction.
#'
#' @export
#'
#' @examples
#' G <- matrix(c(0,0,1, 0,2,1, 1,0,1, 0,1,0, 1,0,0), 3, 5)
#' result <- psd_fit_em(G, 2, 1e-5, 10)
#' psd_loglikelihood(G, result)
psd_loglikelihood <- function (G, result)
{
  rcpp_psd_loss(G, result$P, result$F)
}

#' @title Compute the loglikelihood of svi method in validation set
#'
#' @importFrom Rcpp evalCpp
#'
#' @description Compute the maximum loglikelihood function of the PSD model in validation set using svi method.
#'
#' @param G The I x J matrix of counts; all entries of G should
#'     be taken from \{0,1,2\}.
#' @param K An integer 2 or greater giving the matrix rank.
#' @param P The \code{P} in the output list of \code{psd_fit_svi}.
#' @param maxiter The maximum number of iterations in validation set.
#'
#' @return A value indicates the accuracy of the prediction in validation set.
#'
#' @export
#'
#' @examples
#' # Refer to Articles in AwesomePackage.
psd_loglikelihood_svi <- function (G, K, P, maxiter = 2000)
{
  I <- nrow(G)
  J <- ncol(G)
  F <- rand(K, J)
  iter <- 0
  repeat
  {
    iter <- iter + 1
    F <- rcpp_update_f_val(G, P, F, 1e-9)
    if (! (iter < maxiter)) {break}
  }
  loss <- rcpp_psd_loss(G, P, F)
  return(loss)
}
