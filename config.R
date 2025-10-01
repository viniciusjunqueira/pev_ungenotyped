# ============================================
# Configuration for PEV Analysis
# ============================================

CONFIG <- list(
  # Random seed for reproducibility
  seed = 789,
  
  # Number of replicates per problem size
  n_replicates = 5,
  
  # Output directory
  output_dir = "results",
  
  # Test mode (smaller problem sizes for quick verification)
  test_mode = FALSE,
  
  # Variance components
  sigma2_e = 1.0,
  sigma2_u = 0.5,
  
  # Genomic blending parameters
  alpha = 0.95,  # Weight on genomic relationship matrix
  beta = 0.05    # Weight on pedigree relationship matrix
)

# Problem sizes for full analysis
PROBLEM_SIZES_FULL <- data.frame(
  size_label = c("Small", "Medium", "Large", "Very Large"),
  n_founders = c(30, 50, 100, 150),
  n_sires = c(10, 15, 25, 40),
  n_dams = c(20, 35, 75, 110),
  n_progeny = c(150, 500, 1000, 2000),
  n_genotyped = c(100, 300, 600, 1200),
  n_young = c(30, 100, 200, 400),
  n_markers = c(5000, 8000, 10000, 12000)
)

# Problem sizes for test mode (faster execution)
PROBLEM_SIZES_TEST <- data.frame(
  size_label = c("Small", "Medium"),
  n_founders = c(20, 30),
  n_sires = c(8, 10),
  n_dams = c(12, 20),
  n_progeny = c(80, 150),
  n_genotyped = c(50, 100),
  n_young = c(20, 30),
  n_markers = c(2000, 3000)
)

# Select problem sizes based on mode
get_problem_sizes <- function() {
  if (CONFIG$test_mode) {
    message("Running in TEST MODE with smaller problem sizes")
    return(PROBLEM_SIZES_TEST)
  } else {
    return(PROBLEM_SIZES_FULL)
  }
}