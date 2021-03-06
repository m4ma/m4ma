# Minds for Mobile Agents

<!-- badges: start -->
[![R-CMD-check](https://github.com/m4ma/m4ma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/m4ma/m4ma/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An R package for Rcpp compatible functions used in the Minds for Mobile Agents project.

## Installation

You can install m4ma from GitHub using `devtools`:

```r
install.packages('devtools')

devtools::install_github('m4ma/m4ma')

```

## Example

Currently, m4ma only supports likelihood estimation. You can estimate the log likelihood of an example trace using the `m4ma::msumlogLike_rcpp()` function:

```r
# Set path to test file and load file
filepath = file.path('tests', 'testthat', 'data', 'trace_i.rda')

trace_name = load(filepath)

# Define nests and alpha lists
nests = list(
  Central = c(0, 6, 17, 28),
  NonCentral = c(0:33)[-c(6, 17, 28)],
  acc = c(1:11),
  const = c(12:22),
  dec = c(0, 23:33)
)

alpha = list(
  Central = rep(1/3, 4),
  NonCentral = c(1/3, rep(0.5, 4), 1/3, rep(0.5, 9), 1/3, 
                 rep(0.5, 9), 1/3, rep(0.5, 5)),
  acc = c(rep(0.5, 4), 1, 1/3, rep(0.5, 5)),
  const = c(rep(0.5, 4), 1, 1/3, rep(0.5, 5)),
  dec = c(1/3, rep(0.5, 4), 1, 1/3, rep(0.5, 5))
)

# Get subject parameter matrix
p = attr(get(trace_name), 'pMat')

# Compute log likelihood of trace given subject parameters
m4ma::msumlogLike_rcpp(trace_i, p, minLike = 1e-10, mult = -1)

# 176.7388

```
