# Minds for Mobile Agents

<!-- badges: start -->
[![RSD](https://img.shields.io/badge/rsd-m4ma-00a3e3.svg)](https://research-software-directory.org/software/m4ma)
[![R-CMD-check](https://github.com/m4ma/m4ma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/m4ma/m4ma/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/m4ma/m4ma/branch/development/graph/badge.svg?token=PWCVRIDAH7)](https://codecov.io/gh/m4ma/m4ma)
<!-- badges: end -->

An R package containing C++ implementations to speed up the simulation and parameter estimation of the Predictive Pedestrian model.

## How to Use m4ma

This package is currently not self-contained but should be used in combination with code from the [predped](https://github.com/CharlotteTanis/predped) repository. The m4ma package includes C++ implementations that can be used instead of R code from predped. The functions in m4ma have in most cases the same names as functions in predped, so that they can be easily substituted. Exceptions are functions to estimate parameters and likelihoods. Benchmarks that show the speed improvement of m4ma implementations compared to predped can be found [here](https://github.com/m4ma/m4ma-performance/tree/main/bench).

## Installation

You can install m4ma from GitHub using `devtools`:

```r
install.packages('devtools')

devtools::install_github('m4ma/m4ma')

```

## Getting Started

For access to the predped repository, please contact c.c.tanis@uva.nl.

### Simulating the Predictive Pedestrian Model
The easiest way to substitute predped R functions to simulate the Predictive Pedestrian model with m4ma C++ implementations is by first loading the m4ma package:

```r
library(m4ma)
```

Then, predped functions that should be replaced are removed from the global environment (if they exist):

```r
rm(list = c('example_function_name'))
```

Instead of the removed predped function, the remaining predped simulation code will call the function with the same name from m4ma. It is also possible to switch back to the original R implementation (e.g., for comparison) by creating a new environment and setting a `use` variable:

```
predped_env <- new.env()
predped_env$use <- 'r' # or use = 'cpp'
```

By default, C++ implementations are used.

**Note**: This only works for reimplemented functions that are included both in m4ma and predped.

### Estimating the Likelihood of Simulation Results

The result of the predped simulation is called *trace*. You can estimate the log likelihood of an example trace using the `m4ma::msumlogLike_rcpp()` function:

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

# Get nest indices for cells
cell_nest = m4ma::get_cell_nest()

# Transform trace into format for C++ processing
trace_rcpp = m4ma::create_rcpp_trace(get(trace_name))

# Compute log likelihood of trace given subject parameters
m4ma::msumlogLike(p, trace_rcpp, nests, alpha, cell_nest)

# 176.7388

```

Note that the estimation in m4ma requires a transformation of the trace via `m4ma::create_rcpp_trace()`.

## Documentation

The documentation of m4ma is build with [roxygen2](https://roxygen2.r-lib.org/articles/roxygen2.html) and currently only locally available. See `?m4ma` after installing the package.

## Testing

The code in m4ma is automatically tested on Windows, Mac, and Linux (Ubuntu) using GitHub actions and [testthat](https://testthat.r-lib.org/). The test coverage is calculated via [codecov](https://about.codecov.io/) and [covr](https://covr.r-lib.org/). For the entire package and new code, the coverage is required to be 80% or above.

## Maintenance

The package is maintained by Charlotte Tanis (c.c.tanis@uva.nl) and Andrew Heathcote.

## License

The code is licensed under the Apache 2.0 License. This means that m4ma can be used, modified and redistributed for free, even for commercial purposes.

## Credits

The package was developed by the Netherlands eScience Center in collaboration with the Department of Psychological Methods at the University of Amsterdam. The reimplemented code is majorly based on the predped code written by Andrew Heathcote, Charlotte Tanis, and others.
