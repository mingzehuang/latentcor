---
title: 'latentcor: An R Package for estimating latent correlations from mixed data types'
tags:
- R
- Statistics
- Latent Correlation
date: "15 July 2021"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
authors:
- name: Mingze Huang
  orcid: 0000-0003-3919-1564
  affiliation: 1, 2
- name: Christian L. Müller
  orcid: 0000-0002-3821-7083
  affiliation: 3, 4, 5
- name: Irina Gaynanova
  orcid: 0000-0000-0000-0000
  affiliation: 1
affiliations:
- name: Department of Statistics, Texas A& M University, College Station, TX
  index: 1
- name: Department of Economics, Texas A& M University, College Station, TX
  index: 2
- name: Ludwig-Maximilians-Universität München, Germany
  index: 3
- name: Helmholtz Zentrum München, Germany
  index: 4
- name: Flatiron Institute, New York
- index: 5
---


# Summary

We present `latentcor`, an R package for correlation estimation from data with mixed variable types. Mixed variables types, including continuous, binary, ordinal, zero-inflated, or truncated data are routinely collected in many areas of science. Accurate estimation of correlations among such variables is often the first critical step in statistical analysis workflows. However, Pearson correlation as the default choice is not well suited for mixed data types as the underlying normality assumption is violated. The R package `latentcor` utilizes the powerful of semi-parametric latent Gaussian copula models to estimate latent correlations between mixed data types at unprecedented speed and accuracy. In its present form, the package allows to estimate correlations between any of continuous/binary/ternary/zero-inflated (truncated) variable types. The underlying implementation takes advantage of a fast multi-linear interpolation scheme with an efficient choice of interpolation grid points which gives the package a small memory footprint without compromising estimation accuracy. This makes the application of latent correlations readily available in modern data analysis. 

# Statement of need

The popular `cor` function within R package `stats` [@team2013r] allows to compute Pearson's correlation, as well as Kendall's $\tau$ and Spearman's $\rho$.  A faster algorithm for calculation of Kendall's $\tau$ is implemented in the R package `pcaPP` [@croux2013robust]. Pearson's correlation is not appropriate for skewed or ordinal data, and its use leads to invalid inference in those cases. While both Kendall's $\tau$ and Spearman's $\rho$ are more robust measures of association as they are based on ranks, the resulting values do not have correlation interpretation, and can not be used as direct substitutes in statistical methods that require correlation as input (e.g. graphical models estimation). The R package `polycor` [@fox2019poly] is designed for ordinal data and allows to computes polychoric (ordinal/ordinal) and polyserial (ordinal/continuous) correlations based on latent Gaussian model. However, the package does not have functionality for zero-inflated data, nor can it handle skewed continuous measurements as it does not allow for copula transformation. The R package `mixedCCA` [@yoon2020sparse] is based on the latent Gaussian copula model, and has functionality to computer latent correlations between continuous/binary/zero-inflated variable types. However, this functionality is an intermediate step as `mixedCCA` is specifically designed for canonical correlation analysis on mixed data rather than the latent correlation estimation by itself. Furthermore, `mixedCCA` does not allow for ordinal data types. Thus, there is a need for stand-alone R package for computation of latent correlation based on latent Gaussian copula framework that takes into account all variable types (continuous/binary/ordinal/zero-inflated), is computationally efficient and has small memory footprint. The R package `latentcor` is designed to meet this need.

# Background on latent correlations

The estimation of latent correlations consists of three steps: (i) computing Kendall's $\tau$ between each pair of variables; (ii) choosing the bridge function $F(\cdot)$ based on the types of variable pairs, the bridge function connects the Kendall's $\tau$ computed from the data, $\widehat \tau$, to the true underlying correlation $\rho$ via moment equation $\mathbb{E}(\widehat \tau) = F(\rho)$; (iii) calculating estimate of latent correlation by $F^{-1}(\widehat \tau)$. We summarize the references for the explicit form of $F(\cdot)$ for each variable combination as implemented in `latentcor` below.

 
|Type | continuous | binary | zero-inflated (truncated) | ternary |
|-----|----------|----------|----------|----------|
|continuous | @liu2009nonparanormal |- | -| - |
|binary | @fan2017high | @fan2017high | - | - |
|zero-inflated (truncated) | @yoon2020sparse | @yoon2020sparse | @yoon2020sparse | - |
|ternary | @quan2018rank | @quan2018rank | This work[^1] | @quan2018rank |
[^1]: See the accompanying `latentcor` vignette for derivation details.
 
The inversion of the bridge function $F(\cdot)$ can be done in two ways. The original approach (`method = "original"`) relies on numerical inversion for each pair of variables based on uni-root optimization. Since optimization is done separately for each pair, the original approach is computationally expensive when the number of variables is large. Figure \ref{fig:R_nc_org} displays the estimated latent correlations using the original approach versus the true values of underlying latent correlation for ternary/continuous case, the alignment of points around $y=x$ line confirms that the estimation is empirically unbiased. The second approach to invert  $F(\cdot)$ is to use approximation via multi-linear interpolation on pre-calculated fixed grid of points (`method = "approx"`). This idea has been first proposed by [@yoon2021fast], and implemented for continuous/binary/truncated pairs in `mixedCCA`. However, the implementation lacks ternary case, and the specific grid choice creates a large memory footprint. In `latentcor`, we add ternary case and optimize the choice of grid by redefining the bridge functions on a rescaled version of Kendall's $\tau$, where the scaling adopts to the smoothness of underlying $F$. type of variables by simultaneously controlling the approximation error at the same or lower level. As a result, `latentcor` has significantly smaller memory footprint and smaller approximation error compared to `mixedCCA`, and has additional functionality.

Memory footprints (in KB):

 | case | mixedCCA | latentcor |
 |-----|----------|----------|
| binary/continuous | 10.23 | 3.13 |
| binary/binary | 303.06 | 20.48 |
| truncated/continuous | 21.62 | 3.16 |
| truncated/binary | 902.32 | 28.61 | 
| truncated/truncated | 689.78 | 16.16 |
| ternary/continuous | - | 18.52 |
| ternary/binary | - | 110.93 |
| ternary/truncated | - | 191.78 |
| ternary/ternary | - | 1023.13 |

Figure \ref{fig:R_nc_approx} displays the estimated latent correlations using the approximation approach (`method = "approx"`) versus true values of underlying latent correlation for ternary/continuous case. The results are almost indistinguishable from Figure~\ref{fig:R_nc_org} at a fraction of computational cost. For reference, Figure \ref{fig:R_nc_pearson} displays the values obtained by using standard Pearson correlation, which leads to significant estimation bias.

# Usage

A simple example estimating latent correlation is shown below.

```r
library(latentcor)

# Generate two variables of sample size 100
# The first variable is ternary (pi0 = 0.3, pi1 = 0.5, pi2 = 1-0.3-0.5 = 0.2) 
# The second variable is continuous. 
# No copula transformation is applied.
X = GenData(types = c("ter", "con"), XP = list(c(0.3, .5), NA))$X

# Estimate latent correlation matrix with original method
R_nc_org = estR(X = X, types = c("ter", "con"), method = "original")$R

# Estimate latent correlation matrix with aprroximation method
R_nc_approx = estR(X = X, types = c("ter", "con"), method = "approx")$R

# Heatmap for latent correlation matrix.
Heatmap_R_nc_approx = estR(X = X, types = c("ter", "con"), method = "approx", showplot = TRUE)$plotR
```


# Rendered R Figures
Script see: [latentcor_evaluation](https://github.com/mingzehuang/latentcor_evaluation/blob/master/unbias_check.R)

![Estimated latent correlations by `latentcor` with `method = "original"` versus true population latent correlation. \label{fig:R_nc_org}](nc_org.pdf)


![Estimated latent correlations by `latentcor` with `method = "approx"` versus true population latent correlation.\label{fig:R_nc_approx}](nc_approx.pdf)


![Estimated correlations using `cor` function in `stats` package (Pearson correlation) versus true population latent correlation.\label{fig:R_nc_pearson}](nc_pearson.pdf)



# References
