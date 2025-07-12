# Test the clean plotting functions

library(ggplot2)
library(dplyr)
library(patchwork)

# Source the clean plotting functions
source("phyloglm_framework/phyloglm_summary_plots_clean.R")

if (exists("batch_results")) {
  
  # Test 1: Forest plot only
  cat("Creating clean forest plot...\n")
  forest <- create_clean_forest_plot(batch_results)
  ggsave("clean_forest_plot.png", forest, width = 10, height = 6, dpi = 300)
  cat("Saved to clean_forest_plot.png\n")
  
  # Test 2: Complete summary figure
  cat("\nCreating complete summary figure...\n")
  summary_fig <- create_complete_summary_figure(
    batch_results,
    output_file = "clean_summary_figure.png"
  )
  cat("Saved to clean_summary_figure.png\n")
  
  # Test 3: Individual interaction plots
  cat("\nTesting individual interaction plots...\n")
  
  # Forward direction
  if ("FS_CB_TerrWS_Mass" %in% names(batch_results$results)) {
    result <- batch_results$results[["FS_CB_TerrWS_Mass"]]
    if (result$success) {
      int_plot <- create_clean_interaction_plot(
        result,
        predictor_var = "HighConfidence_Coop",
        terr_var = "TerritorialityWeakVsStrong", 
        response_var = "FemaleSong_Agg01",
        predictor_label = "Cooperative Breeding",
        response_label = "Female Song"
      )
      ggsave("clean_interaction_cb_fs.png", int_plot, width = 6, height = 5, dpi = 300)
      cat("Saved CB->FS interaction to clean_interaction_cb_fs.png\n")
    }
  }
  
  # Reverse direction
  if ("CB_FS_TerrWS_Mass" %in% names(batch_results$results)) {
    result <- batch_results$results[["CB_FS_TerrWS_Mass"]]
    if (result$success) {
      int_plot <- create_clean_interaction_plot(
        result,
        predictor_var = "FemaleSong_Agg01",
        terr_var = "TerritorialityWeakVsStrong",
        response_var = "HighConfidence_Coop",
        predictor_label = "Female Song",
        response_label = "Cooperative Breeding"
      )
      ggsave("clean_interaction_fs_cb.png", int_plot, width = 6, height = 5, dpi = 300)
      cat("Saved FS->CB interaction to clean_interaction_fs_cb.png\n")
    }
  }
  
} else {
  cat("No batch_results found. Please load your results first.\n")
}