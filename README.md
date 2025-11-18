<!-- badges: start -->
[![R-CMD-check](https://github.com/mingzehuang/latentcor/workflows/R-CMD-check/badge.svg)](https://github.com/mingzehuang/latentcor/actions)
[![codecov](https://codecov.io/gh/mingzehuang/latentcor/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mingzehuang/latentcor)
[![CRAN status](https://www.r-pkg.org/badges/version-last-release/latentcor)](https://CRAN.R-project.org/package=latentcor)
[![Launch binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mingzehuang/latentcor/master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03634/status.svg)](https://doi.org/10.21105/joss.03634)
[![DOI](https://zenodo.org/badge/336304814.svg)](https://zenodo.org/badge/latestdoi/336304814)
<!-- badges: end -->


# latentcor: Latent Correlation for Mixed Types of Data

`latentcor` is an `R` package for estimation of latent correlations with mixed data types (continuous, binary, truncated, and ternary) under the latent Gaussian copula model. For references on the estimation framework, see

  * [Fan, J., Liu, H., Ning, Y., and Zou, H. (2017), “High Dimensional Semiparametric Latent Graphical Model for Mixed Data.” *JRSS B*](https://arxiv.org/abs/1404.7236). **Continuous/binary** types.

  * [Quan X., Booth J.G. and Wells M.T."Rank-based approach for estimating correlations in mixed ordinal data." *arXiv*](https://arxiv.org/abs/1809.06255) **Ternary** type.

  * [Yoon G., Carroll R.J. and Gaynanova I. (2020). “Sparse semiparametric canonical correlation analysis for data of mixed types”. *Biometrika*](https://doi.org/10.1093/biomet/asaa007). **Truncated** type for zero-inflated data.

  * [Yoon G., Müller C.L. and Gaynanova I. (2021). “Fast computation of latent correlations” *JCGS*](https://doi.org/10.1080/10618600.2021.1882468). **Approximation method of computation**, see [vignette](https://mingzehuang.github.io/latentcor/articles/latentcor.html) for details.

## Statement of Need

No R software package is currently available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. The R package [`latentcor`](https://CRAN.R-project.org/package=latentcor), introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

**Multi-linear interpolation:** Earlier versions of `latentcor` used multi-linear interpolation based on functionality of R package [`chebpol`](https://cran.r-project.org/package=chebpol) written by Simen Gaure. This functionality is needed for faster computations of latent correlations with approximation method. However, `chebpol` was removed from CRAN on 2022-02-07. The current version of `latentcor` reuses the multi-linear interpolation part of the  `chebpol` (provided under Artistic-2 license) integrated directly within `latentcor`. To cite multi-linear interpolation only, please use original [`chebpol`](https://cran.r-project.org/package=chebpol).

**Accuracy:** The approximation method for ternary/ternary, truncated(zero-inflated)/ternary, and ternary/binary cases are less accurate close to boundary (zero proportions) due to size limitations of CRAN packages on the pre-stored grid.  If higher accuracy is desired and original method is computationally prohibitive, latencor is also available as [Python package](https://pypi.org/project/latentcor/) with [Github development python version](https://github.com/mingzehuang/latentcor_py)

## Installation

To use `latentcor`, you need to install [`R`](https://cran.r-project.org/). To enhance your user experience, you may use some IDE for it (e.g. [`RStudio`](https://posit.co/)).

The development version of [`latentcor`](https://github.com/mingzehuang/latentcor) is available on GitHub. You can download it with the help of the `devtools` package in `R` as follow:

```r
install.packages("devtools")
devtools::install_github("https://github.com/mingzehuang/latentcor", build_vignettes = TRUE)
```
The stable release version [`latentcor`](https://CRAN.R-project.org/package=latentcor) is available on CRAN. You can download it in `R` as follow:

```r
install.packages("latentcor")
```

## Example

A simple example estimating latent correlation is shown below.

```r
library(latentcor)

# Generate two variables of sample size 100
# The first variable is ternary (pi0 = 0.3, pi1 = 0.5, pi2 = 1-0.3-0.5 = 0.2) 
# The second variable is continuous. 
# No copula transformation is applied.
X = gen_data(n = 1000, types = c("ter", "con"), XP = list(c(0.3, .5), NA))$X

# Estimate latent correlation matrix with the original method
latentcor(X = X, types = c("ter", "con"), method = "original")$R

# Estimate latent correlation matrix with the approximation method
latentcor(X = X, types = c("ter", "con"))$R

# Speed improvement by approximation method compared with original method
library(microbenchmark)
microbenchmark(latentcor(X, types = c("ter", "con"), method = "original"),
               latentcor(X, types = c("ter", "con")))
# Unit: milliseconds
# min     lq     mean    median     uq     max     neval
# 5.3444 5.8301 7.033555 6.06740 6.74975 20.8878   100
# 1.5049 1.6245 2.009371 1.73805 1.99820  5.0027   100
# This is run on Windows 10 with Intel(R) Core(TM) i5-4570 CPU @ 3.20GHz   3.20 GHz

# Heatmap for latent correlation matrix.
latentcor(X = X, types = c("ter", "con"), showplot = TRUE)$plotR
```
Another example with the `mtcars` dataset.

```r
library(latentcor)
# Use build-in dataset mtcars
X = mtcars
# Check variable types for manual determination
apply(mtcars, 2, table)
# Or use built-in get_types function to get types suggestions
get_types(mtcars)

# Estimate latent correlation matrix with original method
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), method = "original")$R
# Estimate latent correlation matrix with approximation method
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"))$R

# Speed improvement by approximation method compared with original method
library(microbenchmark)
microbenchmark(latentcor(mtcars, types = types, method = "original"),
               latentcor(mtcars, types = types, method = "approx"))
# Unit: milliseconds
#  min       lq        mean      median        uq      max    neval
#  201.9872 215.6438   225.30385 221.5364 226.58330 411.4940   100
#   71.8457  75.1681   82.42531  80.1688  84.77845 238.3793    100
# This is run on Windows 10 with Intel(R) Core(TM) i5-4570 CPU @ 3.20GHz   3.20 GHz

# Heatmap for latent correlation matrix with approximation method.
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), showplot = TRUE)$plotR
```

Interactive heatmap see: [interactive heatmap of latent correlations (approx) for mtcars](https://rpubs.com/mingzehuang/797937)

Community Guidelines
--------------------

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution guidelines](https://github.com/mingzehuang/latentcor/blob/master/CONTRIBUTING.md) prior
    to submitting a pull request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/mingzehuang/latentcor/issues).
3.  Contributors must adhere to the [Code of Conduct](https://github.com/mingzehuang/latentcor/blob/master/CODE_OF_CONDUCT.md).

Acknowledgments
--------------

We thank Dr. Grace Yoon for providing implementation details of the [`mixedCCA`](https://github.com/irinagain/mixedCCA) R package.
