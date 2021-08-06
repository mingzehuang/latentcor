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
  orcid: 0000-0002-4116-0268
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

We present `latentcor`, an R package for correlation estimation from data with mixed variable types. Mixed variables types, including continuous, binary, ordinal, zero-inflated, or truncated data are routinely collected in many areas of science. Accurate estimation of correlations among such variables is often the first critical step in statistical analysis workflows. Pearson correlation as the default choice is not well suited for mixed data types as the underlying normality assumption is violated. The concept of semi-parametric latent Gaussian copula models, on the other hand, provides a unifying way to estimate  correlations between mixed data types. The R package `latentcor` comprises a comprehensive list of these models, enabling the estimation of correlations between any of continuous/binary/ternary/zero-inflated (truncated) variable types. The underlying implementation takes advantage of a fast multi-linear interpolation scheme with an efficient choice of interpolation grid points, thus giving the package a small memory footprint without compromising estimation accuracy. This makes latent correlation estimation readily available for modern high-throughput data analysis.

# Statement of need
Currently, there is no software package available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. 
For instance, the popular `cor` function within R package `stats` [@team2013r] only allows to compute Pearson's correlation, Kendall's $\tau$ and Spearman's
$\rho$. A faster algorithm for calculating Kendall's $\tau$ is implemented in the R package `pcaPP` [@croux2013robust]. Pearson's correlation is not appropriate
for skewed or ordinal data, and its use leads to invalid inference in those cases. While the rank-based Kendall's $\tau$ and Spearman's $\rho$ are more robust
measures of association, the resulting values do not have correlation interpretation and can not be used as direct substitutes in statistical methods that
require correlation as input (e.g., graphical models estimation). The R package `polycor` [@fox2019poly] is designed for ordinal data and allows to computes
polychoric (ordinal/ordinal) and polyserial (ordinal/continuous) correlations based on latent Gaussian model. However, the package does not have functionality
for zero-inflated data, nor can it handle skewed continuous measurements as it does not allow for copula transformation. The R package `correlation`
[@makowski2020methods] in the `easystats` collection provides 16 different correlation measures, including polychoric and polyserial correlations. However, 
functionality for correlation estimation from zero-inflated data is lacking. The R package `mixedCCA` [@yoon2020sparse] is based on the latent Gaussian copula 
model and can compute latent correlations between continuous/binary/zero-inflated variable types as an intermediate step for canonical correlation analysis. 
However, `mixedCCA` does not allow for ordinal data types. The R package `latentcor`, introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

# Background on latent correlations

The estimation of latent correlations consists of three steps: 
- computing Kendall's $\tau$ between each pair of variables
- choosing the bridge function $F(\cdot)$ based on the types of variable pairs; the bridge function connects the Kendall's $\tau$ computed from the data, $\widehat \tau$, to the true underlying correlation $\rho$ via moment equation $\mathbb{E}(\widehat \tau) = F(\rho)$;
- computing estimates of latent correlation by $F^{-1}(\widehat \tau)$. 

We summarize the references for the explicit form of $F(\cdot)$ for each variable combination as implemented in `latentcor` below.

+----------------+-----------------------+-----------------+--------------------------+-----------------+
| Type           | continuous            | binary          | ternary                  | zero-inflated\  |
|                |                       |                 |                          | (truncated)     |
+================+=======================+=================+==========================+=================+
| continuous     | @liu2009nonparanormal | \-              | \-                       | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| binary         | @fan2017high          | @fan2017high    | \-                       | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| ternary        | @quan2018rank         | @quan2018rank   | @quan2018rank            | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| zero-inflated\ | @yoon2020sparse       | @yoon2020sparse | See `latentcor`\         | @yoon2020sparse |
| (truncated)    |                       |                 | vignette for derivation\ |                 |
+----------------+-----------------------+-----------------+--------------------------+-----------------+

                 
The inversion of the bridge function $F(\cdot)$ can be done in two ways. The original approach (`method = "original"`) relies on numerical inversion for each pair of variables based on uni-root optimization. Since optimization is done separately for each pair, the original approach is computationally expensive when the number of variables is large. Figure \ref{fig:R_nc_org} displays the estimated latent correlations using the original approach versus the true values of underlying latent correlation for ternary/continuous case, the alignment of points around $y=x$ line confirms that the estimation is empirically unbiased. The second approach to invert  $F(\cdot)$ is to use approximation via multi-linear interpolation on pre-calculated fixed grid of points (`method = "approx"`). This idea has been first proposed by [@yoon2021fast], and implemented for continuous/binary/truncated pairs in `mixedCCA`. However, the implementation lacks ternary case, and the specific grid choice creates a large memory footprint. In `latentcor`, we add ternary case and optimize the choice of grid by redefining the bridge functions on a rescaled version of Kendall's $\tau$, where the scaling adopts to the smoothness of underlying $F$. type of variables by simultaneously controlling the approximation error at the same or lower level. As a result, `latentcor` has significantly smaller memory footprint and smaller approximation error compared to `mixedCCA`, and has additional functionality.

Memory footprints (in KB):

 | case | mixedCCA | latentcor |
 |-----|----------|----------|
| binary/continuous | 10.08 | 4.22 |
| binary/binary | 303.04 | 69.1 |
| truncated/continuous | 20.99 | 6.16 |
| truncated/binary | 907.95 | 92.25 | 
| truncated/truncated | 687.68 | 84.33 |
| ternary/continuous | - | 125.83 |
| ternary/binary | - | 728.3 |
| ternary/truncated | - | 860.9 |
| ternary/ternary | - | 950.61 |

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
Heatmap_R_nc_approx = estR(X = X, types = c("ter", "con"),
                           method = "approx", showplot = TRUE)$plotR
```


# Rendered R Figures
Script see: [latentcor_evaluation](https://github.com/mingzehuang/latentcor_evaluation/blob/master/unbias_check.R)

![Estimated latent correlations by `latentcor` with `method = "original"` versus true population latent correlation. \label{fig:R_nc_org}](nc_org.pdf)


![Estimated latent correlations by `latentcor` with `method = "approx"` versus true population latent correlation.\label{fig:R_nc_approx}](nc_approx.pdf)


![Estimated correlations using `cor` function in `stats` package (Pearson correlation) versus true population latent correlation.\label{fig:R_nc_pearson}](nc_pearson.pdf)



# References
