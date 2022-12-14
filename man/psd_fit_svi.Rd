% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psd_fit_svi.R
\name{psd_fit_svi}
\alias{psd_fit_svi}
\title{Use SVI algorithm to fit PSD model}
\usage{
psd_fit_svi(
  G,
  K,
  epsilon = 1e-05,
  maxiter = 5e+05,
  val_iter = 10000,
  maxdrop = 3,
  maxiter.sample = 100,
  maxiter.val = 2000,
  val_J = 0.05,
  val_I = 0.1,
  tau = 1,
  kappa = 0.5
)
}
\arguments{
\item{G}{The I x J matrix of counts; all entries of G should
be taken from \{0,1,2\}.}

\item{K}{An integer 2 or greater giving the matrix rank.}

\item{epsilon}{Convergence criterion.}

\item{maxiter}{The maximum number of iterations.}

\item{val_iter}{The number of iterations between each validation set sampling.}

\item{maxdrop}{The maximum number of consecutive decreases in the loss function. Beyond this value the loop will stop.}

\item{maxiter.sample}{The maximum number of iterations in the sampling section.}

\item{maxiter.val}{The maximum number of iterations in the validation section.}

\item{val_J}{Sample proportion of SNPs in validation set.}

\item{val_I}{Sample proportion of individuals in validation set.}

\item{tau}{A parameter of the descending direction of SVI algorithm.}

\item{kappa}{A parameter of the descending direction of SVI algorithm.}
}
\value{
A List with the following parameters:
\describe{
\item{\code{P}}{The population scale matrix of the individuals.}
\item{\code{Loss}}{A vector represents the value of the loss function which records once for 10 iterations.}
\item{\code{MaxLoss}}{Maximum loss function value. Unlike other algorithms, we observe the loss function on the validation set. Therefore, monotonicity is not guaranteed, that is, the maximum value does not necessarily occur at the end, so the maximum value needs to be recorded.}
\item{\code{Iterations}}{An integer represents the number of iterations.}}
}
\description{
Fit PSD model with SVI algorithm, and use the loss
function as a stopping criterion.
}
\examples{
# Refer to Articles in AwesomePackage.
}
