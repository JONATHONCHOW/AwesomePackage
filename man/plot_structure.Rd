% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_structure.R
\name{plot_structure}
\alias{plot_structure}
\title{Plot the population proportion}
\usage{
plot_structure(
  P,
  pops = NULL,
  label = NULL,
  map.indiv = NULL,
  map.pop = NULL,
  gap = NULL,
  colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#EE2C2C", "#CC79A7",
    "#8968CD", "#FF83FA", "#EECFA1", "#A52A2A", "#4169E1", "#FFFF00", "#BFEFFF",
    "#FF1493"),
  font.size = 9,
  title = NULL,
  subtitle = NULL
)
}
\arguments{
\item{P}{The proportion matrix.}

\item{pops}{Population order options.}

\item{label}{The original order of individuals. This option is only for data that needs to be grouped.}

\item{map.indiv}{The new order of individuals. This option is only for data that needs to be grouped.}

\item{map.pop}{The order of populations. This option is only for data that needs to be grouped.}

\item{gap}{Gaps between groups. This option is only for data that needs to be grouped.}

\item{colors}{Theme color options.}

\item{font.size}{Font size used in plot.}

\item{title}{Title of the plot, such as "EM", "SQP", "VI", "SVI".}

\item{subtitle}{Subtitle of the plot, such as "K = 2", "K = 3".}
}
\value{
A \code{ggplot} object.
}
\description{
Plot the population proportion of individuals using package ggplot2.
}
\examples{
P <- matrix(c(0.5,0.3,0.8, 0.5,0.7,0.2), 3, 2)
plot_structure(P, title = "FUN")
}
