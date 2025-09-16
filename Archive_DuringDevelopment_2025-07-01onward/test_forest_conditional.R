# Test updated forest plot with conditional averaging
source('create_phylopath_bias_robustness_figure.R')
source('create_forest_plot_with_counts.R')

# Create the forest plot with proper conditional averaging
forest_plot <- create_forest_plot_with_counts(
  rate_to_plot = 'COOP->FS',
  reference_value = 0.556,
  full_dataset_n = 875,
  output_file = 'test_forest_plot_conditional_avg.png'
)

# Verify the calculation for one bias correction
library(dplyr)
results <- extract_phylopath_results()
result <- results[['Geographic (Holarctic)']]

# Manually calculate conditional averages for a few seeds to verify
cat("\nVerifying conditional average calculation:\n\n")

for (seed_val in c(1001, 1002, 1003)) {
  seed_data <- result$detailed_models %>%
    filter(seed == seed_val & delta_CICc < 2)
  
  if (nrow(seed_data) > 0) {
    coef_col <- 'HighConfidence_Coop_to_FemaleSong_Agg01_est'
    valid_rows <- !is.na(seed_data[[coef_col]])
    
    if (sum(valid_rows) > 0) {
      # Calculate weights
      weights <- exp(-0.5 * seed_data$delta_CICc[valid_rows])
      weights <- weights / sum(weights)
      
      # Conditional average
      cond_avg <- sum(seed_data[[coef_col]][valid_rows] * weights)
      
      cat(paste("Seed", seed_val, ": Conditional avg =", round(cond_avg, 4), 
                "(", sum(valid_rows), "models with coefficient)\n"))
    }
  }
}

cat("\nThe forest plot now uses proper conditional averaging!\n")