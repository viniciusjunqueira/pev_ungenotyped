# Prediction Error Variance Comparison for Single-Step GBLUP

## Overview
This repository contains code to reproduce the numerical examples presented in *"Derivation of prediction error variance for non-genotyped individuals in genomic selection"*.

The analysis demonstrates that:
1. Direct MME inversion and Schur complement methods produce numerically identical PEV estimates.
2. The Schur complement approach is computationally more efficient.

## Requirements
- R version ≥ 4.0.0  
- Required packages:
  - MASS  
  - ggplot2  
  - gridExtra  

Install packages:
```r
install.packages(c("MASS", "ggplot2", "gridExtra"))
```

## How to Run

1. **Clone the repository**
   ```bash
   git clone https://github.com/viniciusjunqueira/pev_ungenotyped.git
   cd pev_ungenotyped
   ```

2. **Run the analysis**
   From a terminal:
   ```bash
   Rscript run_analysis.R
   ```

   Alternatively, from within R:
   ```r
   source("run_analysis.R")
   ```

3. **Output**
   - All figures and output files are saved in the `results/` folder.  
   - Console output includes timing and numerical comparisons between methods.  
   - Intermediate objects are created dynamically during execution and stored in memory (not written to disk unless specified in `config.R`).

## Folder Structure
```
pev_ungenotyped/
│
├── config.R               # Configuration parameters (paths, seeds, options)
├── generate_figures.R     # Script for generating PEV comparison plots
├── generate_results.R     # Script for computing and exporting numerical results
├── main_analysis.R        # Core implementation of PEV derivations and methods
├── simulation_functions.R # Supporting simulation and matrix computation functions
├── run_analysis.R         # Entry point script that runs the full workflow
├── results/               # Output figures and results
└── README.md
```

## Reproducibility
To ensure reproducibility, you can set a random seed and capture session info at the end of `run_analysis.R`:
```r
set.seed(123)
sessionInfo()
```

## Citation
If you use this repository, please cite the corresponding paper once published.
