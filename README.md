<!-- badges: start -->
[![R-CMD-check](https://github.com/mingzehuang/latentcor/workflows/R-CMD-check/badge.svg)](https://github.com/mingzehuang/latentcor/actions)
[![codecov](https://codecov.io/gh/mingzehuang/latentcor/branch/master/graph/badge.svg)](https://codecov.io/gh/mingzehuang/latentcor)
<!-- badges: end -->


# latentcor: Latent Correlation for Mixed Types of Data

*latentcor* is an *R* package for estimation of latent correlations with mixed data types (continuous, binary, truncated, and ternary) under the latent Gaussian copula model. For references on the estimation framework, see

  * [Fan, J., Liu, H., Ning, Y., and Zou, H. (2017), “High Dimensional Semiparametric Latent Graphical Model for Mixed Data.” *JRSS B*](https://doi.org/10.1111/rssb.12168). **Continuous/binary** types.

  * [Quan X., Booth J.G. and Wells M.T."Rank-based approach for estimating correlations in mixed ordinal data." *arXiv*](https://arxiv.org/abs/1809.06255) **Ternary** type.

  * [Yoon G., Carroll R.J. and Gaynanova I. (2020). “Sparse semiparametric canonical correlation analysis for data of mixed types”. *Biometrika*](https://doi.org/10.1093/biomet/asaa007). **Truncated** type for zero-inflated data.

  * [Yoon G., Müller C.L. and Gaynanova I. (2021). “Fast computation of latent correlations” *JCGS*](https://doi.org/10.1080/10618600.2021.1882468). **Approximation method of computation**, see [vignette](https://mingzehuang.github.io/latentcor/articles/latentcor.html) for details.

## Installation

To use *latentcor*, you need to install [*R*](https://cran.r-project.org/). To enhance your user experience, you may use some IDE for it (e.g. [*RStudio*](https://www.rstudio.com/)).

The development version of *latentcor* is available on GitHub. You can download it with the help of the *devtools* package in *R* as follow:


```r
install.packages("devtools")
devtools::install_github("https://github.com/mingzehuang/latentcor", build_vignettes = TRUE)
```
## Example

A simple example estimating latent correlation is shown below.

```r
library(latentcor)

# Generate two variables of sample size 100
# The first variable is ternary (pi0 = 0.3, pi1 = 0.5, pi2 = 1-0.3-0.5 = 0.2) 
# The second variable is continuous. 
# No copula transformation is applied.
X = GenData(types = c("ter", "con"), XP = list(c(0.3, .5), NA))$X

# Estimate latent correlation matrix with original method
estR(X = X, types = c("ter", "con"), method = "original")$R

# Estimate latent correlation matrix with aprroximation method
estR(X = X, types = c("ter", "con"))$R

# Heatmap for latent correlation matrix.
estR(X = X, types = c("ter", "con"), showplot = TRUE)$plotR
```
Another example with real data set.

```r
library(latentcor)
# Use build-in dataset mtcars
X = mtcars
# Check variable types
apply(mtcars, 2, table)
# Estimate latent correlation matrix with original method
estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), method = "original")$R
# Estimate latent correlation matrix with approximation method
estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"))$R
# Heatmap for latent correlation matrix with approximation method.
estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), showplot = TRUE)$plotR
```

![Heatmap of of latent correlations (approx) for mtcars](https://github.com/mingzehuang/latentcor/blob/master/paper/heatmap_approx_mtcars.png)

Interactive heatmap see: [interactive heatmap of latent correlations (approx) for mtcars](https://rpubs.com/mingzehuang/797937)
