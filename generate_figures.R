# ============================================
# Generate Figures
# ============================================

library(ggplot2)
library(gridExtra)

#' Generate all figures
generate_figures <- function(results, tables, output_dir) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Define publication-ready theme
  theme_publication <- theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 11)
    )

  # Color palette for scenarios
  scenario_colors <- c("Small" = "#1b9e77", "Medium" = "#d95f02",
                       "Large" = "#7570b3", "Very Large" = "#e7298a")

  # ============================================
  # Figure 1: PEV Comparison - Schur vs Simplified Cases
  # ============================================

  # Calculate common axis limits for both panels
  all_pev_values <- c(tables$individual$PEV_schur,
                      tables$individual$PEV_case1,
                      tables$individual$PEV_case2)
  axis_min <- floor(min(all_pev_values) * 20) / 20  # Round down to nearest 0.05

  axis_max <- ceiling(max(all_pev_values) * 20) / 20  # Round up to nearest 0.05

  # Panel A: Schur vs Case 1 (No Fixed Effects)
  p_schur_case1 <- ggplot(tables$individual,
                          aes(x = PEV_schur, y = PEV_case1, color = size)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "black", linewidth = 0.8) +
    scale_color_manual(values = scenario_colors) +
    scale_x_continuous(limits = c(axis_min, axis_max)) +
    scale_y_continuous(limits = c(axis_min, axis_max)) +
    labs(x = expression(bold("Schur Complement")),
         y = expression(bold("Case 1: No Fixed Effects")),
         title = "(A) Schur Complement vs Case 1") +
    coord_fixed(ratio = 1) +
    theme_publication

  # Panel B: Schur vs Case 2 (No Phenotypic Data)
  p_schur_case2 <- ggplot(tables$individual,
                          aes(x = PEV_schur, y = PEV_case2, color = size)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "black", linewidth = 0.8) +
    scale_color_manual(values = scenario_colors) +
    scale_x_continuous(limits = c(axis_min, axis_max)) +
    scale_y_continuous(limits = c(axis_min, axis_max)) +
    labs(x = expression(bold("Schur Complement")),
         y = expression(bold("Case 2: No Phenotypes")),
         title = "(B) Schur Complement vs Case 2") +
    coord_fixed(ratio = 1) +
    theme_publication

  # Combine panels horizontally (side by side)
  combined_plot <- grid.arrange(
    p_schur_case1 + theme(legend.position = "none"),
    p_schur_case2 + theme(legend.position = "none"),
    ncol = 2
  )

  # Extract legend
  legend <- cowplot::get_legend(
    p_schur_case1 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom",
            legend.box.margin = margin(t = 10))
  )

  # Final combined figure with legend at bottom
  final_plot <- grid.arrange(
    combined_plot,
    legend,
    nrow = 2,
    heights = c(10, 1)
  )

  # Save publication figure (wide format for side-by-side panels)
  ggsave(file.path(output_dir, "figure_pev_comparison.pdf"), final_plot, width = 12, height = 6, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "figure_pev_comparison.tiff"), final_plot, width = 12, height = 6, dpi = 600, compression = "lzw", bg = "white")
  
  # ============================================
  # Figure 2: Timing Comparison (Direct vs Schur only)
  # ============================================

  p_timing <- ggplot(tables$timing, aes(x = size, y = mean_time,
                                        color = method, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = pmax(0, mean_time - sd_time),
                      ymax = mean_time + sd_time),
                  width = 0.15, linewidth = 0.6) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Scenario",
         y = "Computation Time (seconds)") +
    theme_publication +
    theme(legend.position = "top",
          legend.title = element_blank())

  ggsave(file.path(output_dir, "timing_comparison.pdf"), p_timing,
         width = 8, height = 5, dpi = 600)
  ggsave(file.path(output_dir, "timing_comparison.tiff"), p_timing,
         width = 8, height = 5, dpi = 600, compression = "lzw")

  # ============================================
  # Figure 2b: Timing Comparison (All Methods)
  # ============================================

  if (!is.null(tables$timing_all)) {
    # Update method labels for publication
    timing_all_plot <- tables$timing_all
    timing_all_plot$method <- factor(timing_all_plot$method,
      levels = c("Direct", "Schur", "Case1_NoFixedEffects", "Case2_NoPhenotypes"),
      labels = c("Direct Inversion", "Schur Complement",
                 "Case 1: No Fixed Effects", "Case 2: No Phenotypes"))

    p_timing_all <- ggplot(timing_all_plot, aes(x = size, y = mean_time,
                                          color = method, group = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, mean_time - sd_time),
                        ymax = mean_time + sd_time),
                    width = 0.15, linewidth = 0.6) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Scenario",
           y = "Computation Time (seconds)") +
      theme_publication +
      theme(legend.position = "top",
            legend.title = element_blank())

    ggsave(file.path(output_dir, "figure_timing_all_methods.pdf"), p_timing_all,
           width = 8, height = 5, dpi = 600)
    ggsave(file.path(output_dir, "figure_timing_all_methods.tiff"), p_timing_all,
           width = 8, height = 5, dpi = 600, compression = "lzw")
  }

  # ============================================
  # Figure 3: Detailed scatter with facets (cases as rows, scenarios as cols)
  # ============================================

  # Reshape data for faceted plot
  individual_long <- reshape(tables$individual,
    direction = "long",
    varying = list(c("PEV_case1", "PEV_case2")),
    v.names = "PEV_simplified",
    timevar = "method",
    times = c("Case 1: No Fixed Effects", "Case 2: No Phenotypes"))

  individual_long$method <- factor(individual_long$method,
    levels = c("Case 1: No Fixed Effects", "Case 2: No Phenotypes"))

  individual_long$size <- factor(individual_long$size,
    levels = c("Small", "Medium", "Large", "Very Large"))

  p_faceted <- ggplot(individual_long,
                      aes(x = PEV_schur, y = PEV_simplified)) +
    geom_point(aes(color = size), alpha = 0.6, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "black", linewidth = 0.6) +
    facet_grid(method ~ size) +
    scale_color_manual(values = scenario_colors) +
    labs(x = "PEV (Schur Complement)",
         y = "PEV (Simplified Method)") +
    coord_fixed(ratio = 1) +
    theme_publication +
    theme(legend.position = "none",
          panel.spacing = unit(0.8, "lines"))

  # Squared figure: 2 rows x 4 cols
  ggsave(file.path(output_dir, "figure_pev_faceted.pdf"), p_faceted,
         width = 10, height = 6, dpi = 600)
  ggsave(file.path(output_dir, "figure_pev_faceted.tiff"), p_faceted,
         width = 10, height = 6, dpi = 600, compression = "lzw")

  message(sprintf("Figures saved to: %s/", output_dir))
}
