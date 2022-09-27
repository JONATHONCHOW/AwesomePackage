#' @title Simulate data of PSD model
#'
#' @description Simulate gene data of PSD model including P, F, G.
#'
#' @param I The number of individuals to simulate.
#' @param J The number of SNPs to simulate.
#' @param K The number of populations to simulate.
#' @param type.id Choose type. Should be one of \code{"A"}, \code{"B"}.
#' @param parm_alpha A parameter of the normal distribution in the simulation of the first type of P.
#' @param parm_sd A parameter of the normal distribution in the simulation of the second type of P.
#' @param parm_F A parameter of the beta distribution in the simulation of F. It has the same length as \code{K}.
#' @param data The data needed to be provided in order to collect the suballele frequencies of real SNPs during the simulation of F.
#'
#' @return A List with the following parameters:
#' \describe{
#' \item{\code{G}}{The I x J simulated matrix of counts.}
#' \item{\code{P}}{The I x K simulated population scale matrix of the individuals.}
#' \item{\code{F}}{The K x J simulated gene scale matrix of the populations.}}
#'
#' @export
#'
#' @examples
#' data_simuA <- psd_simulation(60, 250, 3, type.id = "A", parm_F = c(0.1, 0.05, 0.01))
#' plot_structure(data_simuA$P, pops = rep(3:1))
#' data_simuB <- psd_simulation(100, 500, 5, type.id = "B")
#' plot_structure(data_simuB$P, pops = rep(5:1))
psd_simulation <- function (I, J, K, type.id = c("A", "B"), parm_alpha = 0.1, parm_sd = 2, parm_F = NULL, data = data_HGDP)
{
  if (is.null(parm_F))
  {
    parm_F <- rep(0.1, K)
  }
  if (type.id == "A")
  {
    P <- psd_simulation1_P(I, K, parm_alpha)
  }
  if (type.id == "B")
  {
    P <- psd_simulation2_P(I, K, parm_sd)
  }
  F <- psd_simulation_F(J, K, parm_F, data)
  G <- psd_simulation_G(P, F)
  # Delete cols with the same col value.
  same_col <- vector()
  for(j in 1:ncol(G))
  {
    if(length(unique(G[, j])) == 1)
    {
      same_col <- c(same_col, j)
    }
  }
  G <- G[, -same_col]
  F <- F[, -same_col]
  return(list(G = G, P = P, F = F))
}

# Simulate the first type of P.
#
#' @importFrom MCMCpack rdirichlet
psd_simulation1_P <- function (I, K, parm_alpha)
{
  P <- matrix(0, 1, K)
  temp <- rdirichlet(I, rep(parm_alpha, K))
  for (k in 1:K)
  {
    assign(paste("max", k, sep=""), temp[which(max.col(temp) == k), ])
    P <- rbind(P, get(paste("max", k, sep=""))[order(get(paste("max", k, sep=""))[, k], decreasing = T), ])
  }
  P <- P[-1, ]
  return(P)
}

# Simulate the second type of P.
#
#' @importFrom stats dnorm
psd_simulation2_P <- function (I, K, parm_sd)
{
  P <- matrix(0, I, K)
  indiv <- rep(0:(I - 1)) * (K + 1) / (I - 1)
  pop <- rep(1:K)
  for (i in 1:I)
  {
    temp <- vector()
    for (k in 1:K)
    {
      temp <- c(temp, dnorm(indiv[i], pop[k], parm_sd))
    }
    P[i, ] <- temp / sum(temp)
  }
  return(P)
}

# Simulate F using beta distribution.
#
#' @importFrom stats rbeta
psd_simulation_F <- function (J, K, parm_F, data)
{
  F <- matrix(0, K, J)
  data_sample <- data[, sample(c(1:ncol(data)), J)]
  mean <- colSums(data_sample) / nrow(data_sample) / 2
  for (k in 1:K)
  {
    beta1 <- (1 - parm_F[k]) / parm_F[k] * mean
    beta2 <- (1 - parm_F[k]) / parm_F[k] * (1 - mean)
    for (j in 1:J)
    {
      F[k, j] <- rbeta(1, beta1[j], beta2[j])
    }
  }
  return(F)
}

# Simulate G from P and F using binomial distribution.
#
#' @importFrom stats rbinom
psd_simulation_G <- function (P, F)
{
  I <- nrow(P)
  J <- ncol(F)
  G <- matrix(0, I, J)
  PF <- P %*% F
  for (i in 1:I)
  {
    for (j in 1:J)
    {
      G[i, j] <- rbinom(1, 2, PF[i, j])
    }
  }
  return(G)
}
