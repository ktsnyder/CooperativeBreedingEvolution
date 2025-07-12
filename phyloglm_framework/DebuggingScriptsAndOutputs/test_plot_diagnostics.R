# Diagnostic script to understand plot issues

# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Source plotting functions
source("phyloglm_framework/phyloglm_summary_plots.R")

# Function to diagnose batch results
diagnose_batch_results <- function(batch_results) {
  cat("\n=== Batch Results Diagnostic ===\n")
  
  # Check structure
  cat("\nBatch results structure:\n")
  cat("- Number of analyses:", length(batch_results$results), "\n")
  cat("- Analysis names:\n")
  
  for (name in names(batch_results$results)) {
    result <- batch_results$results[[name]]
    if (result$success) {
      cat("  ", name, ":\n")
      cat("    Response:", result$config$response, "\n")
      cat("    Predictors:", paste(result$config$predictors, collapse = ", "), "\n")
      cat("    Controls:", paste(result$config$controls, collapse = ", "), "\n")
      cat("    N models:", length(result$models$models), "\n")
      cat("    Best model:", result$comparison$comparison$Model[1], "\n")
    }
  }
  
  # Look for bidirectional analyses
  cat("\n=== Checking for Bidirectional Analyses ===\n")
  
  # Forward: CB -> FS
  forward_found <- FALSE
  for (name in names(batch_results$results)) {
    result <- batch_results$results[[name]]
    if (result$success && 
        result$config$response == "FemaleSong_Agg01" &&
        "HighConfidence_Coop" %in% result$config$predictors &&
        "TerritorialityWeakVsStrong" %in% result$config$predictors) {
      cat("\nForward analysis found:", name, "\n")
      forward_found <- TRUE
    }
  }
  
  # Reverse: FS -> CB
  reverse_found <- FALSE
  for (name in names(batch_results$results)) {
    result <- batch_results$results[[name]]
    if (result$success && 
        result$config$response == "HighConfidence_Coop" &&
        "FemaleSong_Agg01" %in% result$config$predictors &&
        "TerritorialityWeakVsStrong" %in% result$config$predictors) {
      cat("\nReverse analysis found:", name, "\n")
      reverse_found <- TRUE
    }
  }
  
  if (!forward_found) cat("\nWARNING: No forward analysis (FS ~ CB + Terr) found!\n")
  if (!reverse_found) cat("\nWARNING: No reverse analysis (CB ~ FS + Terr) found!\n")
  
  # Check for problematic analyses
  cat("\n=== Checking for Analyses to Exclude ===\n")
  for (name in names(batch_results$results)) {
    if (grepl("Terr3", name)) {
      cat("- Found Territory_12vs3 analysis (should be excluded):", name, "\n")
    }
    result <- batch_results$results[[name]]
    if (result$success && 
        !any(c("TerritorialityWeakVsStrong", "Territory") %in% result$config$predictors)) {
      cat("- Found analysis without territoriality:", name, "\n")
    }
  }
}

# Test individual plot functions
test_individual_plots <- function(batch_results) {
  cat("\n=== Testing Individual Plot Functions ===\n")
  
  # Test forest plot
  cat("\nTesting forest plot...\n")
  tryCatch({
    forest <- create_bidirectional_forest_plot(batch_results)
    ggsave("test_forest_plot_only.png", forest, width = 10, height = 6)
    cat("Forest plot saved to test_forest_plot_only.png\n")
  }, error = function(e) {
    cat("ERROR in forest plot:", e$message, "\n")
  })
  
  # Test finding specific analyses
  cat("\nLooking for specific analyses for interaction plots...\n")
  
  # Find a CB->FS analysis
  cb_fs_analysis <- NULL
  for (name in names(batch_results$results)) {
    result <- batch_results$results[[name]]
    if (result$success && 
        result$config$response == "FemaleSong_Agg01" &&
        "HighConfidence_Coop" %in% result$config$predictors &&
        "TerritorialityWeakVsStrong" %in% result$config$predictors) {
      cb_fs_analysis <- name
      break
    }
  }
  
  if (!is.null(cb_fs_analysis)) {
    cat("Found CB->FS analysis:", cb_fs_analysis, "\n")
    tryCatch({
      int_plot <- create_interaction_plot(
        batch_results$results[[cb_fs_analysis]],
        predictor_var = "HighConfidence_Coop",
        terr_var = "TerritorialityWeakVsStrong",
        response_var = "FemaleSong_Agg01",
        predictor_label = "Cooperative Breeding",
        response_label = "Female Song"
      )
      ggsave("test_interaction_cb_fs.png", int_plot, width = 8, height = 6)
      cat("Interaction plot saved to test_interaction_cb_fs.png\n")
    }, error = function(e) {
      cat("ERROR in interaction plot:", e$message, "\n")
    })
  }
}

# Main diagnostic function
if (exists("batch_results")) {
  diagnose_batch_results(batch_results)
  test_individual_plots(batch_results)
} else {
  cat("No batch_results found. Load your results first:\n")
  cat('batch_results <- readRDS("path/to/all_results.rds")\n')
}