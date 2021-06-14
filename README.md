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
## Multivariate case
n = 1000; p1 = 8; p2 = 5 # sample size and dimensions for two datasets.
rho = .9 # Autocorrelated coefficient.
# Data generation
simdata = GenData(n=n, type1 = "ternary", type2 = "continuous", p1 = p1, p2 = p2, 
                  rho = rho, copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 =  NULL)
X1 = simdata$X1; X2 = simdata$X2
# Estimate latent correlation matrix with original method
R_nc_org = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "approx")$R

# Heatmap for two latent correlation matrix for two groups of variables.
LatentPlot(R_nc_approx)
```
![\label{Heatmap for Latent Correlation Matrix}](https://rpubs.com/mingzehuang/781462)
