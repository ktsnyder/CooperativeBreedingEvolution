# Script to run stepwise model expansion on best phyloglm models

# Load required functions
source("phyloglm_framework/stepwise_model_expansion.R")
source("phyloglm_framework/stepwise_visualization.R")

# Run stepwise expansion using the combined results with 1000 bootstrap iterations
cat("Running stepwise model expansion...\n")
cat("=================================\n\n")

# Run the batch expansion
expansion_results <- run_stepwise_batch(
  results_file = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  target_analyses = c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass"),
  output_dir = "Outputs/PhyloglmResults/stepwise_expansion",
  aic_threshold = 2,  # Require AIC improvement > 2
  n_bootstrap = 100,  # Bootstrap final models
  verbose = TRUE
)

# Create visualizations
cat("\n\nCreating visualizations...\n")

# Load the saved results (in case script was interrupted)
expansion_results <- readRDS("Outputs/PhyloglmResults/stepwise_expansion/all_expansion_results.rds")

# Create summary figure
create_stepwise_summary_figure(
  expansion_results,
  output_file = "Outputs/PhyloglmResults/stepwise_expansion/stepwise_summary"
)

# Create comparison table
comparison_table <- create_model_comparison_table(
  expansion_results,
  output_file = "Outputs/PhyloglmResults/stepwise_expansion/model_comparison.csv"
)

cat("\nModel comparison:\n")
print(comparison_table)

# Create predictor importance plot
importance_plot <- plot_predictor_importance(expansion_results)
ggsave("Outputs/PhyloglmResults/stepwise_expansion/predictor_importance.png",
       importance_plot, width = 10, height = 8, dpi = 300)

# Print summary for each analysis
cat("\n\n=================================\n")
cat("SUMMARY OF RESULTS\n")
cat("=================================\n")

for (name in names(expansion_results)) {
  cat("\n\n", name, "\n")
  cat(rep("-", nchar(name) + 2), "\n", sep = "")
  print(expansion_results[[name]])
}

cat("\n\nAll results saved to: Outputs/PhyloglmResults/stepwise_expansion/\n")