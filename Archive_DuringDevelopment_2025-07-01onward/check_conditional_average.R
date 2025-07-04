# Check how to properly do conditional averaging
source('create_phylopath_bias_robustness_figure.R')

# Extract results
results <- extract_phylopath_results()

# Look at one example
bias_name <- "Geographic (Holarctic)"
result <- results[[bias_name]]

# Check the structure of detailed_models
cat("Column names in detailed_models:\n")
print(names(result$detailed_models))

# Check a few rows to understand the data
cat("\nFirst few rows of detailed_models:\n")
print(head(result$detailed_models[, c("seed", "model", "delta_CICc", "CICc", "CICc_weight")], 10))

# Check if CICc_weight is available
if ("CICc_weight" %in% names(result$detailed_models)) {
  cat("\nCICc_weight column is available for proper conditional averaging\n")
} else {
  cat("\nNeed to calculate CICc weights from delta_CICc\n")
}

# Look for the coefficient column
coef_col <- grep("Coop.*to.*FemaleSong.*est", names(result$detailed_models), value = TRUE)[1]
cat(paste("\nCoefficient column:", coef_col, "\n"))

# Check one iteration to see how conditional averaging should work
seed_1_data <- result$detailed_models[result$detailed_models$seed == 1 & result$detailed_models$delta_CICc < 2, ]
cat(paste("\nFor seed 1, models with delta_CICc < 2:", nrow(seed_1_data), "\n"))
if (nrow(seed_1_data) > 0) {
  print(seed_1_data[, c("model", "delta_CICc", "CICc_weight", coef_col)])
}