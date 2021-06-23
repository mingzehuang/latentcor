<!-- badges: start -->
[![R-CMD-check](https://github.com/mingzehuang/latentcor/workflows/R-CMD-check/badge.svg)](https://github.com/mingzehuang/latentcor/actions)
[![codecov](https://codecov.io/gh/mingzehuang/latentcor/branch/master/graph/badge.svg)](https://codecov.io/gh/mingzehuang/latentcor)
<!-- badges: end -->


# latentcor: Latent Correlation for Mixed Types of Data


*latentcor* is an *R* package provides estimation for latent correlation with mixed data types (continuous, binary, truncated and ternary).

## Installation

To use *latentcor*, you need to install [*R*](https://cran.r-project.org/). To enhance your user experience, you may use some IDE for it (e.g. [*RStudio*](https://rstudio.com/)).

The development version of *latentcor* is available on GitHub. You can download it with the help of the *devtools* package in *R* as follow:


```r
install.packages("devtools")
devtools::install_github("https://github.com/mingzehuang/latentcor", build_vignettes = TRUE)
```
## Example

```r
# Data generation
X = GenData(types = c("ter", "con"), XP = list(c(0.3, .5), NA))
# Estimate latent correlation matrix with original method
R_nc_org = estR(X = X, types = c("ter", "con"), method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_org = estR(X = X, types = c("ter", "con"), method = "approx")$R

# Heatmap for latent correlation matrix.
Heatmap_R_nc_approx = estR(X = X, types = c("ter", "con"), method = "approx", corplot = TRUE)$plotR
```
