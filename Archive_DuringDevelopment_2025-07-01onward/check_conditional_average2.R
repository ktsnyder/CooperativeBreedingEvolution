# Check how to properly do conditional averaging
source('create_phylopath_bias_robustness_figure.R')

# Extract results
results <- extract_phylopath_results()

# Look at one example
bias_name <- "Geographic (Holarctic)"
result <- results[[bias_name]]

# Check a few rows to understand the data
cat("\nFirst few rows of detailed_models:\n")
print(head(result$detailed_models[, c("seed", "model", "delta_CICc", "CICc")], 10))

# Calculate CICc weights for conditional averaging
# For models with delta_CICc < 2
library(dplyr)

# Check one iteration to see how conditional averaging should work
seed_1_data <- result$detailed_models %>%
  filter(seed == 1 & delta_CICc < 2)

cat(paste("\nFor seed 1, models with delta_CICc < 2:", nrow(seed_1_data), "\n"))

if (nrow(seed_1_data) > 0) {
  # Calculate CICc weights using the formula from phylopath
  # weight = exp(-0.5 * delta_CICc) / sum(exp(-0.5 * delta_CICc))
  seed_1_data$CICc_weight <- exp(-0.5 * seed_1_data$delta_CICc)
  seed_1_data$CICc_weight <- seed_1_data$CICc_weight / sum(seed_1_data$CICc_weight)
  
  coef_col <- grep("Coop.*to.*FemaleSong.*est", names(seed_1_data), value = TRUE)[1]
  
  cat("\nModels with their weights:\n")
  print(seed_1_data[, c("model", "delta_CICc", "CICc_weight", coef_col)])
  
  # Calculate conditional average
  cond_avg <- sum(seed_1_data[[coef_col]] * seed_1_data$CICc_weight, na.rm = TRUE)
  cat(paste("\nConditional average for seed 1:", round(cond_avg, 4), "\n"))
  
  # Compare to simple mean
  simple_mean <- mean(seed_1_data[[coef_col]], na.rm = TRUE)
  cat(paste("Simple mean for seed 1:", round(simple_mean, 4), "\n"))
  cat(paste("Difference:", round(cond_avg - simple_mean, 4), "\n"))
}