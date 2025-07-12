# Test the simplified plotting functions

library(ggplot2)
library(dplyr)
library(patchwork)

# Load batch results
batch_results <- readRDS('./Outputs/PhyloglmResults/PhyloGLM_Batch_20250710_235142_boot100/all_results.rds')

# Source the simplified plotting functions
source("phyloglm_framework/phyloglm_summary_plots_simple.R")

# Test 1: Forest plot only
cat("Creating simplified forest plot...\n")
tryCatch({
  forest <- create_simple_forest_plot(batch_results)
  ggsave("simple_forest_plot.png", forest, width = 10, height = 6, dpi = 300)
  cat("Saved to simple_forest_plot.png\n")
}, error = function(e) {
  cat("Error in forest plot:", e$message, "\n")
})

# Test 2: Complete summary figure
cat("\nCreating complete summary figure...\n")
tryCatch({
  summary_fig <- create_simple_summary_figure(
    batch_results,
    output_file = "simple_summary_figure.png"
  )
  cat("Saved to simple_summary_figure.png\n")
}, error = function(e) {
  cat("Error in summary figure:", e$message, "\n")
})

# Test 3: Individual interaction plots
cat("\nTesting individual interaction plots...\n")

# Forward direction
if ("FS_CB_TerrWS_Mass" %in% names(batch_results)) {
  result <- batch_results[["FS_CB_TerrWS_Mass"]]
  if (result$success) {
    tryCatch({
      int_plot <- create_simple_interaction_plot(
        result,
        predictor_var = "HighConfidence_Coop",
        response_var = "FemaleSong_Agg01",
        predictor_label = "Cooperative Breeding",
        response_label = "Female Song"
      )
      ggsave("simple_interaction_cb_fs.png", int_plot, width = 6, height = 5, dpi = 300)
      cat("Saved CB->FS interaction to simple_interaction_cb_fs.png\n")
    }, error = function(e) {
      cat("Error in CB->FS interaction:", e$message, "\n")
    })
  }
}

# Reverse direction
if ("CB_FS_TerrWS_Mass" %in% names(batch_results)) {
  result <- batch_results[["CB_FS_TerrWS_Mass"]]
  if (result$success) {
    tryCatch({
      int_plot <- create_simple_interaction_plot(
        result,
        predictor_var = "FemaleSong_Agg01",
        response_var = "HighConfidence_Coop",
        predictor_label = "Female Song",
        response_label = "Cooperative Breeding"
      )
      ggsave("simple_interaction_fs_cb.png", int_plot, width = 6, height = 5, dpi = 300)
      cat("Saved FS->CB interaction to simple_interaction_fs_cb.png\n")
    }, error = function(e) {
      cat("Error in FS->CB interaction:", e$message, "\n")
    })
  }
}

cat("\nDone!\n")