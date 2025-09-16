# Debug Territoriality Level 3 extraction
source('create_phylopath_bias_robustness_figure.R')

# Extract results
results <- extract_phylopath_results()

# Check Year-round Territoriality
if ("Year-round Territoriality" %in% names(results)) {
  terr3 <- results[["Year-round Territoriality"]]
  cat("Year-round Territoriality found in results\n")
  cat(paste("Number of rows:", nrow(terr3$detailed_models), "\n"))
  
  # Check for COOP->FS coefficient
  coef_cols <- grep("Coop.*to.*FemaleSong", names(terr3$detailed_models), value = TRUE)
  if (length(coef_cols) > 0) {
    coef_col <- coef_cols[1]
    coef_values <- terr3$detailed_models[[coef_col]]
    valid_values <- coef_values[!is.na(coef_values)]
    cat(paste("\nCoefficient column:", coef_col, "\n"))
    cat(paste("Sample of coefficient values:", paste(head(valid_values), collapse = ", "), "\n"))
    cat(paste("Mean coefficient:", round(mean(valid_values), 3), "\n"))
  }
} else {
  cat("Year-round Territoriality NOT found in results\n")
  cat("Available keys:", paste(names(results), collapse = ", "), "\n")
}