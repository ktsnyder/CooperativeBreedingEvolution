# Test updated forest plot function
source('create_phylopath_bias_robustness_figure.R')
source('create_forest_plot_with_counts.R')

# Create the forest plot for COOP->FS with updated parameters
forest_plot <- create_forest_plot_with_counts(
  rate_to_plot = 'COOP->FS',
  reference_value = 0.556,
  full_dataset_n = 875,
  output_file = 'test_forest_plot_COOP_FS_updated.png'
)

# Also check that we're using conditional averages
results <- extract_phylopath_results()
cat("\nChecking that we're using conditional averages from models with deltaCICc < 2:\n\n")

# Check one example
bias_name <- "Geographic (Holarctic)"
result <- results[[bias_name]]
if (!is.null(result$detailed_models)) {
  # Check how many models have deltaCICc < 2 per iteration
  model_counts <- result$detailed_models %>%
    filter(delta_CICc < 2) %>%
    group_by(seed) %>%
    summarise(n_models = n(), .groups = "drop")
  
  cat(paste("For", bias_name, ":\n"))
  cat(paste("  Average models with deltaCICc < 2 per iteration:", 
            round(mean(model_counts$n_models), 2), "\n"))
  cat(paste("  Range:", min(model_counts$n_models), "-", 
            max(model_counts$n_models), "\n"))
}