
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AwesomePackage

<!-- badges: start -->
<!-- badges: end -->

The goal of AwesomePackage is to infer ancestry with some models and fit
those models with some algorithms.

We use the classical PSD model for ancestor inference, which has been
widely used, such as ADMIXTURE, FRAPPE. We illustrate the close
relationship between the PSD model, the Poisson NMF model and the
multinomial topic model, which can optimize the algorithm. We use
Expectation-Maximization algorithm (EM), co-ordinate descent algorithm
(CD), sequential quadratic programming algorithm (SQP) and stochastic
variational inference algorithm (SVI) to fit the model, and illustrate
the differences between these algorithms through simulation experiments
and real data experiments.

## Installation

You can install the development version of AwesomePackage from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JONATHONCHOW/AwesomePackage")
```

## Example

You can see some examples at Vignettes.
