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
No R software package is currently available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. 
The popular `cor` function within R package `stats` [@team2013r], for instance, allows to compute Pearson's correlation, Kendall's $\tau$ and Spearman's
$\rho$, and a faster algorithm for calculating Kendall's $\tau$ is implemented in the R package `pcaPP` [@croux2013robust]. Pearson's correlation is not
appropriate for skewed or ordinal data, and its use leads to invalid inference in those cases. While the rank-based Kendall's $\tau$ and Spearman's $\rho$ are
more robust measures of association, the resulting values do not have correlation interpretation and can not be used as direct substitutes in statistical methods
that require correlation as input (e.g., graphical model estimation [@yoon2019microbial]). The R package `polycor` [@fox2019poly] is designed for ordinal data and allows to computes
polychoric (ordinal/ordinal) and polyserial (ordinal/continuous) correlations based on latent Gaussian model. However, the package does not have functionality
for zero-inflated data, nor can it handle skewed continuous measurements as it does not allow for copula transformation. The R package `correlation`
[@makowski2020methods] in the `easystats` collection provides 16 different correlation measures, including polychoric and polyserial correlations. However, 
functionality for correlation estimation from zero-inflated data is lacking. The R package `mixedCCA` [@yoon2020sparse] is based on the latent Gaussian copula 
model and can compute latent correlations between continuous/binary/zero-inflated variable types as an intermediate step for canonical correlation analysis. 
However, `mixedCCA` does not allow for ordinal data types. The R package `latentcor`, introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

# Estimation of latent correlations

## The general estimation workflow

The estimation of latent correlations consists of three steps: 

* computing Kendall's $\tau$ between each pair of variables,

* choosing the bridge function $F(\cdot)$ based on the types of variable pairs; the bridge function connects the Kendall's $\tau$ computed from the data, $\widehat \tau$, to the true underlying correlation $\rho$ via moment equation $\mathbb{E}(\widehat \tau) = F(\rho)$;

* estimating latent correlation by calculating $F^{-1}(\widehat \tau)$. 

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


## Efficient inversion of the bridge function

In `latentcor`, the inversion of the bridge function $F(\cdot)$ can be computed in two ways. The original approach (`method = "original"`) relies on numerical
inversion for each pair of variables based on uni-root optimization [@yoon2020sparse]. Since each pair of variables requires a separate optimization run, the
original approach is computationally expensive when the number of variables is large. The second approach to invert $F(\cdot)$ is through fast multi-linear
interpolation of pre-calculated $F^{-1}$ values at specific sets of interpolation grid points (`method = "approx"`). This construction has been proposed in
[@yoon2021fast] and is available for continuous/binary/truncated pairs in the current version of `mixedCCA`. However, that implementation lacks the ternary
variable case and relies on an interpolation grid with a large memory footprint. `latentcor` includes the ternary case and provides an optimized interpolation 
grid by redefining the bridge functions on a rescaled version of Kendall's $\tau$. Here, the scaling adapts to the smoothness of the underlying type of variables 
by simultaneously controlling the approximation error at the same or lower level. As a result, `latentcor` has significantly smaller memory footprint (see
Table below) and smaller approximation error compared to `mixedCCA`.

\newpage

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

## Illustrative example 

To illustrate the excellent performance of latent correlation estimation on mixed data, we consider the simple example of estimating correlations between continuous and ternary variables. In this synthetic scenario, we have access to the true underlying correlation between the variables. Figure \ref{fig:R_all}A displays the values obtained by using standard Pearson correlation, revealing a significant estimation bias with respect to the true correlations. Figure \ref{fig:R_all}B displays the estimated latent correlations using the original approach versus the true values of underlying ternary/continuous correlations. 
The alignment of points around $y=x$ line confirms that the estimation is empirically unbiased. Figure \ref{fig:R_all}C displays the estimated latent correlations using the approximation approach (`method = "approx"`) versus true values of underlying latent correlation. The results are almost indistinguishable from Figure \ref{fig:R_all}B at a fraction of the computational cost.

![Scatter plots of estimated Pearson correlation (panel A) and latent correlations (original in panel B, approximate in panel C) vs. ground truth correlations \label{fig:R_all}](./CombinedCorrelations.pdf)

The script to reproduce the displayed results is available at [latentcor_evaluation](https://github.com/mingzehuang/latentcor_evaluation/blob/master/unbias_check.R).

# Basic Usage

We provide two basic code examples of how to use `latentcor` in R. 

The first example illustrates how to estimate latent correlation from pairs of ternary/continuous variables.

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
The second example considers the `mtcars` dataset, available in standard R. 
The `mtcars` dataset comprises continuous, binary, and ternary data types.

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
# Heatmap for latent correlation matrix.
estR(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), showplot = TRUE)$plotR
```

![Heatmap of pearson correlation, latent correlation (original) and latent correlations (approx) for mtcars \label{fig:R_cars}](./all_heatmap.pdf).

The script to reproduce Figure \ref{fig:R_cars} is available at [latentcor_cars](https://github.com/mingzehuang/latentcor_evaluation/blob/master/all_heatmap.R).

We also provide interactive heatmaps for [Pearson correlation for mtcars](https://rpubs.com/mingzehuang/797945), [latent correlation (original) for mtcars](https://rpubs.com/mingzehuang/797939), and [latent correlation (approx) for mtcars](https://rpubs.com/mingzehuang/797937).

# Availability

The R package 'latentcor' is available on [Github](https://github.com/mingzehuang/latentcor/). A comprehensive vignette with additional mathematical and computational details is available [here](https://mingzehuang.github.io/latentcor/articles/latentcor.html).

# Acknowledgments

We thank Dr. Grace Yoon for providing implementation details of the `mixedCCA` R package.

# References
