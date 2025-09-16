# create_final_publication_figure.R
# Create final publication-ready figure combining key visualizations

library(ggplot2)
library(dplyr)
library(patchwork)

# Source the main functions
source('create_phylopath_bias_robustness_figure.R')

# Create individual plots
cat("Creating enhanced DAG...\n")
dag <- create_enhanced_dag(
  output_file = "final_enhanced_dag.png",
  title = "Full Dataset Phylogenetic Path Analysis"
)

cat("Creating model consistency heatmap...\n")
heatmap <- create_model_consistency_heatmap(
  output_file = "final_model_consistency_heatmap.png",
  title = "Model Consistency Across Bias Corrections"
)

cat("Creating forest plot with species counts...\n")
forest <- create_forest_plot_with_counts(
  output_file = "final_forest_plot.png",
  reference_value = 0.556  # CB→FS from full dataset
)

cat("Creating all paths horizontal plot...\n")
# Extract results for a representative bias correction
results <- extract_phylopath_results()
if ("Geographic (Holarctic)" %in% names(results)) {
  all_paths <- create_all_paths_horizontal_plot(
    detailed_models = results[["Geographic (Holarctic)"]]$detailed_models,
    title = "Path Coefficient Distributions - Geographic (Holarctic) Bias Correction",
    output_file = "final_all_paths_plot.png"
  )
}

# Create a comprehensive multi-panel figure
cat("Creating comprehensive multi-panel figure...\n")

# Combine plots with patchwork
comprehensive_figure <- (dag | heatmap) / (forest | all_paths) +
  plot_annotation(
    title = "Robustness of the cooperative breeding-female song association",
    subtitle = "Analysis across multiple bias corrections demonstrates consistent positive association",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  ) +
  plot_layout(heights = c(1, 1))

# Save the comprehensive figure
ggsave("final_comprehensive_robustness_figure.png", 
       comprehensive_figure, 
       width = 20, height = 16, dpi = 300)

ggsave("final_comprehensive_robustness_figure.pdf", 
       comprehensive_figure, 
       width = 20, height = 16)

cat("\nAll figures created successfully!\n")
cat("Individual plots:\n")
cat("- final_enhanced_dag.png\n")
cat("- final_model_consistency_heatmap.png\n")
cat("- final_forest_plot.png\n")
cat("- final_all_paths_plot.png\n")
cat("\nCombined figure:\n")
cat("- final_comprehensive_robustness_figure.png/pdf\n")