source('create_phylopath_bias_robustness_figure.R')
library(dplyr)

results <- extract_phylopath_results()
result <- results[['Geographic (Holarctic)']]

# Check with seed 1001
seed_data <- result$detailed_models %>%
  filter(seed == 1001 & delta_CICc < 2)

cat(paste('For seed 1001, models with delta_CICc < 2:', nrow(seed_data), '\n'))

if (nrow(seed_data) > 0) {
  # Calculate CICc weights
  seed_data$CICc_weight <- exp(-0.5 * seed_data$delta_CICc)
  seed_data$CICc_weight <- seed_data$CICc_weight / sum(seed_data$CICc_weight)
  
  coef_col <- 'HighConfidence_Coop_to_FemaleSong_Agg01_est'
  
  cat('\nModels with their weights:\n')
  print(seed_data[, c('model', 'delta_CICc', 'CICc_weight', coef_col)])
  
  # Calculate conditional average
  # Only use rows where the coefficient is not NA
  valid_rows <- !is.na(seed_data[[coef_col]])
  if (sum(valid_rows) > 0) {
    # Re-normalize weights for valid rows only
    valid_weights <- seed_data$CICc_weight[valid_rows]
    valid_weights <- valid_weights / sum(valid_weights)
    
    cond_avg <- sum(seed_data[[coef_col]][valid_rows] * valid_weights)
    cat(paste('\nConditional average:', round(cond_avg, 4), '\n'))
    
    simple_mean <- mean(seed_data[[coef_col]][valid_rows])
    cat(paste('Simple mean:', round(simple_mean, 4), '\n'))
    cat(paste('Difference:', round(cond_avg - simple_mean, 4), '\n'))
  }
}