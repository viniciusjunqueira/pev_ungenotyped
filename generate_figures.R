# ============================================
# Generate Figures
# ============================================

library(ggplot2)
library(gridExtra)

#' Generate all figures
generate_figures <- function(results, tables, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Timing comparison
  p_timing <- ggplot(tables$timing, aes(x = size, y = mean_time, 
                                        color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), 
                  width = 0.2) +
    labs(#title = "Computing Time Comparison",
         x = "Scenario", y = "Time (seconds)", color = "Method") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")
          , legend.position = 'top'
          , legend.title = element_blank()
          , text = element_text(size = 20))
  
  # PEV scatterplots
  p_direct_schur <- ggplot(tables$individual, 
                           aes(x = PEV_direct, y = PEV_schur, color = size)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = "Direct MME vs Schur Complement",
         x = "PEV (Direct)", y = "PEV (Schur)") +
    theme_minimal(base_size = 12)
  
  p_schur_pa <- ggplot(tables$individual, 
                       aes(x = PEV_schur, y = PEV_PA, color = size)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = "Schur Complement vs Parent Average",
         x = "PEV (Schur)", y = "PEV (PA)") +
    theme_minimal(base_size = 12)
  
  # Save figures
  ggsave(file.path(output_dir, "timing_comparison.pdf"), p_timing, 
         width = 8, height = 6)
  
  pev_plots <- grid.arrange(p_direct_schur, p_schur_pa, ncol = 2)
  ggsave(file.path(output_dir, "pev_comparison_plots.pdf"), pev_plots, 
         width = 12, height = 5)
  
  message(sprintf("Figures saved to: %s/", output_dir))
}
