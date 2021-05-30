---
title: 'Latentcor: An R Package for Latent Correlation Estimation'
tags:
- R
- Statistics
- Latent Correlation
date: "18 May 2021"
output:
  rticles::joss_article: default
  rticles::joss_article(): default
authors:
- name: Mingze Huang
  orcid: 0000-0003-3919-1564
  affiliation: 1, 2
- name: Irina Gaynanova
  orcid: 0000-0000-0000-0000
  affiliation: 2
- name: Christian L. Muller
year: 2021
bibliography: paper.bib
affiliations:
- name: Department of Statistics, Texas A& M University
  index: 1
- name: Department of Economics, Texas A& M University
  index: 2
- name: Department of Statistics, University of Munich
  index: 3
csl: apa.csl
journal: JOSS
---

# Summary

The R package *latentcor* provides estimation for latent correlation with mixed data types (continuous, binary, truncated and ternary). Comparing to *MixedCCA*, which estimates latent correlation for canonical correlation analysis, our new package provides a standalone version for latent correlation estimation. Also we add new functionality for latent correlation between ternary/continous, ternary/binary, ternary/truncated and ternary/ternary cases.
Compare to MixedCCA, standalone, new functionality, memory footprint.

# Statement of need
No package deal with latent correlation across mixed data type.

# Usage
 Create a table
 type of variable (reference to papers)
 Definition for Kendall tau.
 Describe the theorem states expected value of kendall tau = F(r). refer to table for reference of formula.
 A sentence to describe two methods: slow original vs interpolation.
 Table to show memory improvement compare to mixedCCA.
 Generate data for pairs.
 One example of original call, one for interpolation call.
# Citations


# Rendered R Figures

# Acknowledgements

# References
