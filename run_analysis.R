# ============================================
# Master Script: Run Complete Analysis
# ============================================

# Load configuration
source("config.R")

# Check and install required packages
required_packages <- c("MASS", "ggplot2", "gridExtra", "cowplot")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  message("Installing required packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

# Load scripts
source("main_analysis.R")
source("generate_results.R")
source("generate_figures.R")

# Run analysis
message(strrep("=", 50))
message("PEV Comparison Analysis for Single-Step GBLUP")
message(strrep("=", 50))

problem_sizes <- get_problem_sizes()

start_time <- Sys.time()

results <- run_pev_analysis(
  problem_sizes = problem_sizes,
  n_replicates = CONFIG$n_replicates,
  seed = CONFIG$seed,
  sigma2_e = CONFIG$sigma2_e,
  sigma2_u = CONFIG$sigma2_u,
  alpha = CONFIG$alpha,
  beta = CONFIG$beta
)

tables <- generate_result_tables(results, CONFIG$output_dir)

# Create a copy of timing table with all methods for the new figures
tables$timing_all <- tables$timing
tables$timing_all$size <- factor(x = tables$timing_all$size, levels = c('Small'
                                                                        ,'Medium'
                                                                        ,'Large'
                                                                        ,'Very Large'))

# Keep original timing table with Direct/Schur comparison only
tables$timing <- tables$timing[tables$timing$method %in% c('Direct', 'Schur'), ]
tables$timing$size <- factor(x = tables$timing$size, levels = c('Small'
                                                                ,'Medium'
                                                                ,'Large'
                                                                ,'Very Large'))
tables$timing$method <- ifelse(test = tables$timing$method == 'Direct', 'Direct Inversion', 'Schur Complement')

generate_figures(results, tables, CONFIG$output_dir)

 end_time <- Sys.time()

message("\n", strrep("=", 50))
message("Analysis complete!")
message(sprintf("Total time: %.1f minutes", 
                as.numeric(difftime(end_time, start_time, units = "mins"))))
message(sprintf("Results saved to: %s/", CONFIG$output_dir))
message(strrep("=", 50))