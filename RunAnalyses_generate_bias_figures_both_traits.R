# Generate phylopath bias robustness figures for both trait combinations
# Uses already-generated "detailed_models" results files
# This script creates publication-ready figures for the phylopath analysis with downsampling

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

# Define your trait sets
trait_sets <- list(
  weak_vs_strong = "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET",
  terr_12vs3 = "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET"
)

# Specify which iteration count to use (change this as needed)
n_iterations_to_use <- 500  # Use 500 for final figures, 10 for testing

# Process each trait set
for (trait_name in names(trait_sets)) {
  trait_set <- trait_sets[[trait_name]]
  
  cat("\n========================================\n")
  cat("Processing:", trait_name, "\n")
  cat("Trait set:", trait_set, "\n")
  cat("Using", n_iterations_to_use, "iterations\n")
  cat("========================================\n")
  
  # Extract results for specified iteration count
  results <- extract_phylopath_results(
    results_dir = "Outputs",
    trait_set = trait_set,
    n_iterations = n_iterations_to_use
  )
  
  if (length(results) == 0) {
    cat("WARNING: No results found for", n_iterations_to_use, "iterations\n")
    cat("Trying without iteration filter...\n")
    
    # Try without iteration filter
    results <- extract_phylopath_results(
      results_dir = "Outputs",
      trait_set = trait_set
    )
  }
  
  cat("Found", length(results), "bias corrections\n")
  
  if (length(results) > 0) {
    # 1. Create model consistency heatmap
    cat("\nCreating model consistency heatmap...\n")
    heatmap <- create_bias_model_consistency_heatmap(
      bias_results_list = results,
      top_n = 20,
      trait_set = trait_set,
      title = paste0("Model Consistency Across Bias Corrections\n", 
                     "(", trait_name, ", n=", n_iterations_to_use, " iterations)")
    )
    cat("✓ Heatmap saved\n")
    
    # 2. Create forest plots for important rates
    rates_to_plot <- c("COOP->FS", "FS->COOP", "TERR->FS", "TERR->COOP", "MASS->COOP", "MASS->FS", "MASS->TERR")
    
    cat("\nCreating forest plots...\n")
    for (rate in rates_to_plot) {
      cat("  -", rate, "...")
      forest <- create_forest_plot_with_counts(
        bias_results_list = results,
        rate_to_plot = rate,
        trait_set = trait_set
      )
      cat(" ✓\n")
    }
    
    # 3. Report output location
    output_dir <- file.path("Outputs", "PhylopathFigures", get_trait_set_label(trait_set))
    cat("\nAll figures saved to:\n", output_dir, "\n")
    
  } else {
    cat("\nERROR: No results found for this trait set\n")
  }
}

# Summary
cat("\n\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
cat("Processed", length(trait_sets), "trait combinations\n")
cat("Used", n_iterations_to_use, "iterations for each\n")
cat("\nTo use different iteration counts:\n")
cat("  - Change 'n_iterations_to_use' at the top of this script\n")
cat("  - Use 10 for quick testing\n")
cat("  - Use 500 for publication figures\n")
cat("\nOutput structure:\n")
cat("  Outputs/PhylopathFigures/\n")
cat("    ├── FemaleSong_Agg01_HighConfidence_Coop_TerritorialityWeakVsStrong_logMass_AVONET/\n")
cat("    │   ├── model_consistency_heatmap_byrate.pdf\n")
cat("    │   ├── forest_plot_COOP_FS_with_counts.pdf\n")
cat("    │   ├── forest_plot_TERR_FS_with_counts.pdf\n")
cat("    │   └── forest_plot_TERR_COOP_with_counts.pdf\n")
cat("    └── FemaleSong_Agg01_HighConfidence_Coop_Territory_12vs3_logMass_AVONET/\n")
cat("        └── (same files)\n")