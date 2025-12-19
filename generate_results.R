# ============================================
# Generate Results Tables
# ============================================

#' Calculate comparison statistics
calc_comparison_stats <- function(method1, method2, method1_name, method2_name) {
  diff <- method1 - method2
  data.frame(
    comparison = paste(method1_name, "vs", method2_name),
    mean_method1 = mean(method1),
    mean_method2 = mean(method2),
    mean_diff = mean(diff),
    min_diff = min(diff),
    max_diff = max(diff),
    sd_diff = sd(diff),
    rmse = sqrt(mean(diff^2)),
    correlation = cor(method1, method2)
  )
}

#' Generate all result tables
generate_result_tables <- function(results, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Timing summary
  timing_summary <- do.call(rbind, lapply(names(results$timing), function(size) {
    data <- results$timing[[size]]
    data$size <- size
    data
  }))
  
  timing_final <- aggregate(time_seconds ~ size + method, timing_summary,
                            function(x) c(mean = mean(x), sd = sd(x)))
  timing_final <- data.frame(
    size = timing_final$size,
    method = timing_final$method,
    mean_time = timing_final$time_seconds[, 1],
    sd_time = timing_final$time_seconds[, 2]
  )
  
  write.csv(timing_summary, file.path(output_dir, "timing_results_all_replicates.csv"), 
            row.names = FALSE)
  write.csv(timing_final, file.path(output_dir, "timing_summary.csv"), 
            row.names = FALSE)
  
  # PEV comparison statistics
  comparison_table <- do.call(rbind, lapply(names(results$individual), function(size) {
    pev_data <- results$individual[[size]]
    rbind(
      calc_comparison_stats(pev_data$PEV_schur, pev_data$PEV_case1,
                            "Schur", "Case1_NoFixedEffects"),
      calc_comparison_stats(pev_data$PEV_schur, pev_data$PEV_case2,
                            "Schur", "Case2_NoPhenotypes"),
      calc_comparison_stats(pev_data$PEV_case1, pev_data$PEV_case2,
                            "Case1_NoFixedEffects", "Case2_NoPhenotypes")
    ) |> transform(size = size)
  }))
  
  write.csv(comparison_table, file.path(output_dir, "pev_comparison_statistics.csv"), 
            row.names = FALSE)
  
  # Individual PEV values
  all_individual <- do.call(rbind, results$individual)
  write.csv(all_individual, file.path(output_dir, "individual_pev_values.csv"), 
            row.names = FALSE)
  
  message(sprintf("Results saved to: %s/", output_dir))
  
  return(list(
    timing = timing_final,
    comparison = comparison_table,
    individual = all_individual
  ))
}