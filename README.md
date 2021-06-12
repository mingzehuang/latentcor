
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
library(latentcor)
### Data setting
n = 1000 # sample size
Sigma = autocor(p1 + p2, 0.7)

# Data generation
simdata = GenData(n=n, type1 = "binary", type2 = "continuous",
copula1 = "exp", copula2 = "cube",  muZ = mu, Sigma = Sigma,
c1 = matrix(rep(1, p1), ncol = p1), c2 =  NULL)
```

```r
X1 = simdata$X1; X2 = simdata$X2
# Estimate latent correlation matrix with original method
R_nc_org = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous",
                              method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous",
                              method = "approx")$R
```
