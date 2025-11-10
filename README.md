# riemstats: Riemannian ANOVA Statistics

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/yourusername/riemstats/workflows/R-CMD-check/badge.svg)](https://github.com/nicoesve/riemstats/actions) -->
<!-- badges: end -->

## Overview

`riemstats` provides statistical tools for performing ANOVA on Riemannian manifolds. The package implements permutation-based testing procedures for comparing populations of symmetric positive definite (SPD) matrices and other manifold-valued data.

## Installation

You can install the development version of riemstats with:

```r
# install.packages("devtools")
devtools::install_github("nicoesve/riemstats")
```

## Usage

Here's a basic example using the permutation test for Riemannian ANOVA:

```r
library(riemstats)
library(riemtan)

# Create sample data on SPD manifold
data <- CSample$new(metric = "AIRM")
# ... add your SPD matrices and group labels ...

# Perform Riemannian ANOVA with permutation test
result <- riem_anova(data, nperm = 1000)
print(result$pvalue)
```

## Features

- **Permutation-based ANOVA**: Exact inference for comparing groups on Riemannian manifolds
- **Multiple test statistics**: Log-Wilks lambda, Pillai's trace
- **Harmonization methods**: ComBat and rigid alignment for batch effect correction
- **Normalization utilities**: Tools for preprocessing manifold-valued data

## Documentation
For more detailed information, check out:

* [Package website](https://nicoesve.github.io/riemstats/)
* [Function reference](https://nicoesve.github.io/riemstats/reference/)
* [Vignettes](https://nicoesve.github.io/riemstats/articles/)

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This package is licensed under the MIT License.
