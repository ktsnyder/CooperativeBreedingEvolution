# Example script for running phylopath bias robustness analysis
# for both trait combinations with specific iteration counts

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

# Define your two trait set combinations
trait_set_weak_vs_strong <- "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
trait_set_12vs3 <- "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET"

# ========================================================================
# EXAMPLE 1: Run analysis for TerritorialityWeakVsStrong with 500 iterations
# ========================================================================

cat("\n==== Running analysis for TerritorialityWeakVsStrong (500 iterations) ====\n")

# Extract results for 500 iterations only
results_weak_500 <- extract_phylopath_results(
  results_dir = "Outputs",
  trait_set = trait_set_weak_vs_strong,
  n_iterations = 500  # This will only find files with "_500_" in the name
)

cat("Found", length(results_weak_500), "bias corrections with 500 iterations\n")

# Create bias model consistency heatmap
if (length(results_weak_500) > 0) {
  heatmap_weak_500 <- create_bias_model_consistency_heatmap(
    bias_results_list = results_weak_500,
    trait_set = trait_set_weak_vs_strong,
    title = "Model Consistency - TerritorialityWeakVsStrong (500 iterations)"
  )
  cat("Created heatmap for TerritorialityWeakVsStrong\n")
}

# Create forest plots for key rates
rates_to_analyze <- c("COOP->FS", "TERR->FS", "TERR->COOP")

for (rate in rates_to_analyze) {
  forest_weak_500 <- create_forest_plot_with_counts(
    bias_results_list = results_weak_500,
    rate_to_plot = rate,
    trait_set = trait_set_weak_vs_strong
  )
  cat("Created forest plot for", rate, "\n")
}

# ========================================================================
# EXAMPLE 2: Run analysis for Territory_12vs3 with 10 iterations (practice run)
# ========================================================================

cat("\n==== Running analysis for Territory_12vs3 (10 iterations - practice) ====\n")

# Extract results for 10 iterations only
results_12vs3_10 <- extract_phylopath_results(
  results_dir = "Outputs",
  trait_set = trait_set_12vs3,
  n_iterations = 10  # This will only find files with "_10_" in the name
)

cat("Found", length(results_12vs3_10), "bias corrections with 10 iterations\n")

# Create visualizations for the practice run
if (length(results_12vs3_10) > 0) {
  # Note: You might want to add a note to the title about this being a practice run
  heatmap_12vs3_10 <- create_bias_model_consistency_heatmap(
    bias_results_list = results_12vs3_10,
    trait_set = trait_set_12vs3,
    title = "Model Consistency - Territory_12vs3 (10 iterations - PRACTICE RUN)"
  )
  cat("Created heatmap for Territory_12vs3 (10 iterations)\n")
}

# ========================================================================
# EXAMPLE 3: Run analysis for Territory_12vs3 with 500 iterations (if available)
# ========================================================================

cat("\n==== Checking for Territory_12vs3 (500 iterations) ====\n")

# Try to extract results for 500 iterations
results_12vs3_500 <- extract_phylopath_results(
  results_dir = "Outputs",
  trait_set = trait_set_12vs3,
  n_iterations = 500
)

if (length(results_12vs3_500) > 0) {
  cat("Found", length(results_12vs3_500), "bias corrections with 500 iterations\n")
  
  # Create full analysis
  heatmap_12vs3_500 <- create_bias_model_consistency_heatmap(
    bias_results_list = results_12vs3_500,
    trait_set = trait_set_12vs3,
    title = "Model Consistency - Territory_12vs3 (500 iterations)"
  )
  
  for (rate in rates_to_analyze) {
    forest_12vs3_500 <- create_forest_plot_with_counts(
      bias_results_list = results_12vs3_500,
      rate_to_plot = rate,
      trait_set = trait_set_12vs3
    )
  }
  cat("Created all visualizations for Territory_12vs3 (500 iterations)\n")
} else {
  cat("No 500-iteration results found for Territory_12vs3\n")
  cat("You may need to run the full 500-iteration analysis\n")
}

# ========================================================================
# EXAMPLE 4: Extract ALL results regardless of iteration count
# ========================================================================

cat("\n==== Extracting ALL results (any iteration count) ====\n")

# Don't specify n_iterations to get all available results
results_weak_all <- extract_phylopath_results(
  results_dir = "Outputs",
  trait_set = trait_set_weak_vs_strong
  # n_iterations not specified - will find all files
)

cat("\nFor TerritorialityWeakVsStrong, found results for:\n")
for (bias_name in names(results_weak_all)) {
  n_iter <- length(unique(results_weak_all[[bias_name]]$detailed_models$seed))
  cat("  -", bias_name, ":", n_iter, "iterations\n")
}

# ========================================================================
# NOTES ON FILE SELECTION
# ========================================================================

cat("\n==== How file selection works ====\n")
cat("1. When n_iterations is specified, only files with '_[n]_' in the filename are selected\n")
cat("2. Files are expected to be named like: detailed_models_Remove83HolarcticNoncoop_500_2025-06-11.csv\n")
cat("3. If multiple files match (e.g., different dates), the most recent is used\n")
cat("4. If n_iterations is not specified, all matching files are considered\n")
cat("\n")
cat("Output files will be saved to:\n")
cat("  Outputs/PhylopathFigures/[full_trait_set_name]/\n")
cat("This keeps results organized by trait combination\n")

# ========================================================================
# CHECK FOR FULL DATASET RESULTS
# ========================================================================

cat("\n==== Checking for full dataset results ====\n")

# Check for full dataset result for each trait set
full_weak <- find_or_create_full_dataset_result(trait_set_weak_vs_strong, run_if_missing = TRUE)
full_12vs3 <- find_or_create_full_dataset_result(trait_set_12vs3, run_if_missing = TRUE)

if (is.null(full_weak)) {
  cat("\nNote: No full dataset result found for TerritorialityWeakVsStrong\n")
  cat("The forest plots will not show the reference line without this\n")
}

if (is.null(full_12vs3)) {
  cat("\nNote: No full dataset result found for Territory_12vs3\n")
  cat("The heatmap will not show full dataset model rankings without this\n")
}