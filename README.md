# Prediction Error Variance Comparison for Single-Step GBLUP

## Overview
This repository contains code to reproduce the numerical examples presented in "Derivation of prediction error variance for non-genotyped individuals in genomic selection" submitted to *Genetics*.

The analysis demonstrates that:
  1. Direct MME inversion and Schur complement methods produce numerically identical PEV estimates
2. Schur complement is computationally more efficient
3. Parent average method provides an approximation but systematically differs from MME-based methods

## Requirements
- R version ≥ 4.0.0
- Required packages:
  - MASS
- ggplot2
- gridExtra

Install packages:
  ```r