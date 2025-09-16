# Test that we're loading the correct Territoriality (Level 3) file
source('create_phylopath_bias_robustness_figure.R')
source('create_forest_plot_with_counts.R')
library(phylopath)

# Extract results and check which file is being used
results <- extract_phylopath_results()

# Check Territoriality (Level 3)
terr3_result <- results[["Territoriality (Level 3)"]]

if (!is.null(terr3_result)) {
  cat("Territoriality (Level 3) data loaded successfully\n")
  cat(paste("Number of rows:", nrow(terr3_result$detailed_models), "\n"))
  cat(paste("Number of iterations:", terr3_result$n_iterations, "\n"))
  
  # Check which territoriality variable is in the columns
  terr_cols <- grep("Territor", names(terr3_result$detailed_models), value = TRUE)
  cat("\nTerritoriality-related columns found:\n")
  print(unique(terr_cols))
  
  # Specifically check for the correct variable
  if (any(grepl("TerritorialityWeakVsStrong", terr_cols))) {
    cat("\n✓ Correct! Using TerritorialityWeakVsStrong variable\n")
  } else {
    cat("\n✗ Warning: Not using TerritorialityWeakVsStrong variable\n")
  }
}

# Now create the forest plot with the updated data
cat("\nCreating forest plot with updated Territoriality (Level 3) data...\n")
forest_plot <- create_forest_plot_with_counts(
  rate_to_plot = "COOP->FS",
  output_file = "forest_plot_updated_terr3.png"
)

cat("Forest plot created successfully!\n")