# Test forest plot function
source('create_phylopath_bias_robustness_figure.R')
source('create_forest_plot_with_counts.R')

# Create the forest plot for COOP->FS
forest_plot <- create_forest_plot_with_counts(
  rate_to_plot = 'COOP->FS',
  reference_value = 0.556,
  output_file = 'test_forest_plot_COOP_FS.png'
)

# Check the summary data
results <- extract_phylopath_results()

cat("\nChecking for coefficient columns in each bias correction:\n\n")

for (bias_name in names(results)) {
  result <- results[[bias_name]]
  if (!is.null(result$detailed_models)) {
    # Find COOP->FS coefficient column
    pattern_str <- 'Coop.*to.*FemaleSong.*est'
    coef_columns <- grep(pattern_str, names(result$detailed_models), value = TRUE, ignore.case = TRUE)
    
    if (length(coef_columns) > 0) {
      coef_col <- coef_columns[1]
      cat(paste('Bias:', bias_name, '- Column:', coef_col, '\n'))
      
      # Check for NA values
      coef_values <- result$detailed_models[[coef_col]]
      cat(paste('  Total rows:', nrow(result$detailed_models), '\n'))
      cat(paste('  Non-NA values:', sum(!is.na(coef_values)), '\n'))
      cat(paste('  Unique seeds:', length(unique(result$detailed_models$seed)), '\n'))
      
      # Check if all values are NA
      if (sum(!is.na(coef_values)) == 0) {
        cat('  WARNING: All values are NA!\n')
      }
    } else {
      cat(paste('Bias:', bias_name, '- NO COEFFICIENT COLUMN FOUND\n'))
    }
    cat('\n')
  }
}