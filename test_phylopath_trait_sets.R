# Test script for updated phylopath bias robustness functions
# This script demonstrates how to use the updated functions with different trait combinations

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Source the updated functions
source("phylopath_helper_functions.R")
source("create_phylopath_bias_robustness_figure_updated.R")
source("create_bias_model_consistency_heatmap_updated.R")
source("create_forest_plot_with_counts_updated.R")

# Define the two trait set combinations
trait_set1 <- "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
trait_set2 <- "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET"

# Test 1: Extract results for both trait sets
cat("\n=== Testing extract_phylopath_results() ===\n")

# Extract results for trait set 1 (TerritorialityWeakVsStrong)
cat("\nExtracting results for trait set 1 (TerritorialityWeakVsStrong):\n")
results1 <- extract_phylopath_results(trait_set = trait_set1)
cat("Found", length(results1), "bias corrections\n")
cat("Bias corrections found:", paste(names(results1), collapse = ", "), "\n")

# Extract results for trait set 2 (Territory_12vs3)
cat("\nExtracting results for trait set 2 (Territory_12vs3):\n")
results2 <- extract_phylopath_results(trait_set = trait_set2)
cat("Found", length(results2), "bias corrections\n")
cat("Bias corrections found:", paste(names(results2), collapse = ", "), "\n")

# Test 2: Create bias model consistency heatmaps
cat("\n\n=== Testing create_bias_model_consistency_heatmap() ===\n")

# Create heatmap for trait set 1
cat("\nCreating heatmap for trait set 1 (TerritorialityWeakVsStrong):\n")
if (length(results1) > 0) {
  heatmap1 <- create_bias_model_consistency_heatmap(
    bias_results_list = results1,
    trait_set = trait_set1,
    title = "Model Consistency - TerritorialityWeakVsStrong"
  )
  cat("Heatmap created successfully\n")
  cat("Output file will be saved to:", 
      file.path("Outputs", "PhylopathFigures", get_trait_set_label(trait_set1), 
                "model_consistency_heatmap_byrate.pdf"), "\n")
} else {
  cat("No results found for trait set 1\n")
}

# Create heatmap for trait set 2
cat("\nCreating heatmap for trait set 2 (Territory_12vs3):\n")
if (length(results2) > 0) {
  heatmap2 <- create_bias_model_consistency_heatmap(
    bias_results_list = results2,
    trait_set = trait_set2,
    title = "Model Consistency - Territory_12vs3"
  )
  cat("Heatmap created successfully\n")
  cat("Output file will be saved to:", 
      file.path("Outputs", "PhylopathFigures", get_trait_set_label(trait_set2), 
                "model_consistency_heatmap_byrate.pdf"), "\n")
} else {
  cat("No results found for trait set 2\n")
}

# Test 3: Create forest plots
cat("\n\n=== Testing create_forest_plot_with_counts() ===\n")

# Define rates to plot
rates_to_plot <- c("COOP->FS", "TERR->FS", "TERR->COOP")

# Create forest plots for trait set 1
cat("\nCreating forest plots for trait set 1 (TerritorialityWeakVsStrong):\n")
if (length(results1) > 0) {
  for (rate in rates_to_plot) {
    cat("  Creating forest plot for", rate, "\n")
    forest_plot1 <- create_forest_plot_with_counts(
      bias_results_list = results1,
      rate_to_plot = rate,
      trait_set = trait_set1
    )
    cat("  Output file will be saved to:", 
        file.path("Outputs", "PhylopathFigures", get_trait_set_label(trait_set1), 
                  paste0("forest_plot_", gsub("->", "_", rate), "_with_counts.pdf")), "\n")
  }
} else {
  cat("No results found for trait set 1\n")
}

# Create forest plots for trait set 2
cat("\nCreating forest plots for trait set 2 (Territory_12vs3):\n")
if (length(results2) > 0) {
  for (rate in rates_to_plot) {
    cat("  Creating forest plot for", rate, "\n")
    forest_plot2 <- create_forest_plot_with_counts(
      bias_results_list = results2,
      rate_to_plot = rate,
      trait_set = trait_set2
    )
    cat("  Output file will be saved to:", 
        file.path("Outputs", "PhylopathFigures", get_trait_set_label(trait_set2), 
                  paste0("forest_plot_", gsub("->", "_", rate), "_with_counts.pdf")), "\n")
  }
} else {
  cat("No results found for trait set 2\n")
}

# Test 4: Check for full dataset results
cat("\n\n=== Testing find_or_create_full_dataset_result() ===\n")

cat("\nChecking for full dataset result for trait set 1:\n")
full_result1 <- find_or_create_full_dataset_result(trait_set1, run_if_missing = TRUE)
if (!is.null(full_result1)) {
  cat("Full dataset result found for trait set 1\n")
  cat("Number of species:", full_result1$n_species, "\n")
} else {
  cat("No full dataset result found for trait set 1\n")
}

cat("\nChecking for full dataset result for trait set 2:\n")
full_result2 <- find_or_create_full_dataset_result(trait_set2, run_if_missing = TRUE)
if (!is.null(full_result2)) {
  cat("Full dataset result found for trait set 2\n")
  cat("Number of species:", full_result2$n_species, "\n")
} else {
  cat("No full dataset result found for trait set 2\n")
}

cat("\n\n=== Test Summary ===\n")
cat("All functions have been updated to support multiple trait combinations.\n")
cat("The functions now:\n")
cat("1. Accept a trait_set parameter to specify which combination to use\n")
cat("2. Automatically generate output filenames that include the full trait set\n")
cat("3. Look for results in the appropriate subdirectories\n")
cat("4. Handle both regular downsampling and JackknifeSpecies results\n")
cat("5. Save outputs to trait-specific subdirectories in Outputs/PhylopathFigures/\n")
cat("\nTo use these functions with your data:\n")
cat("1. Make sure your phylopath results are in the expected directory structure\n")
cat("2. Run the full dataset phylopath for each trait combination if needed\n")
cat("3. Call the functions with the appropriate trait_set parameter\n")