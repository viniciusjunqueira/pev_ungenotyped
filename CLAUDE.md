# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R project comparing Prediction Error Variance (PEV) computation methods for non-genotyped individuals in Single-Step GBLUP genomic selection. Implements and benchmarks four methods: Direct MME inversion, Schur complement, Case 1 (no fixed effects), and Case 2 (no phenotypes).

## Commands

Run full analysis:
```bash
Rscript run_analysis.R
```

Run in test mode (faster, smaller problem sizes):
```r
# Edit config.R and set: CONFIG$test_mode = TRUE
# Then run: Rscript run_analysis.R
```

Install dependencies:
```r
install.packages(c("MASS", "ggplot2", "gridExtra", "cowplot"))
```

## Architecture

**Entry Point**: `run_analysis.R` - Orchestrates the full workflow: loads config, sources scripts, runs analysis, generates results and figures.

**Core Components**:
- `config.R` - Configuration parameters (variance components, problem sizes, seeds, test mode toggle)
- `simulation_functions.R` - Matrix building functions:
  - `build_A_matrix()` - Pedigree relationship matrix
  - `build_G_matrix()` - Genomic relationship matrix (VanRaden method)
  - `build_H_inverse()` - Single-step H-inverse matrix
  - `build_pedigree()` - Simulates pedigree with diverse family structures
  - `simulate_genotypes()` - Mendelian inheritance simulation
- `main_analysis.R` - Core `run_pev_analysis()` function implementing all four PEV methods
- `generate_results.R` - CSV output generation (timing, PEV comparisons)
- `generate_figures.R` - Publication-quality ggplot2 figures

**Data Flow**:
1. Build pedigree → 2. Compute A matrix → 3. Simulate genotypes → 4. Build G matrix → 5. Blend G with A22 → 6. Build H-inverse → 7. Compute PEV via all methods → 8. Compare results

**Output**: All results saved to `results/` folder (CSV files and PDF/PNG figures).

## Key Parameters in config.R

- `seed` - Random seed for reproducibility
- `n_replicates` - Replicates per problem size (default: 5)
- `sigma2_e`, `sigma2_u` - Variance components
- `alpha`, `beta` - Genomic blending weights (G_blended = alpha*G + beta*A22)
- `test_mode` - Set TRUE for quick testing with smaller matrices
