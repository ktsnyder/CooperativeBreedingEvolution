# Test the new summary plots

# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Source the new plotting functions
source("phyloglm_framework/phyloglm_summary_plots.R")
source("phyloglm_framework/visualization_framework.R")

# Load previous results (replace with your actual results path)
# For this example, I'll show how to use the functions

# Example 1: If you have batch results already loaded
if (exists("batch_results")) {
  
  # Create the combined summary figure with bidirectional interaction plots
  # This will create Panel A (forest plot) + Panel B (CB->FS interaction) + Panel C (FS->CB interaction)
  summary_fig <- create_phyloglm_summary_figure(
    batch_results,
    predictor_var = "HighConfidence_Coop",
    response_var = "FemaleSong_Agg01",
    predictor_label = "Cooperative Breeding",
    response_label = "Female Song",
    output_file = "phyloglm_summary_figure_3panels.png"
  )
  
  print(summary_fig)
  
  # Create just the forest plot
  forest_plot <- create_bidirectional_forest_plot(
    batch_results,
    predictor_var = "HighConfidence_Coop",
    response_var = "FemaleSong_Agg01",
    predictor_label = "Cooperative Breeding",
    response_label = "Female Song"
  )
  
  ggsave("forest_plot_only.png", forest_plot, width = 8, height = 6, dpi = 300)
  
  # Create model selection heatmap with improved text size
  heatmap <- create_model_selection_heatmap(batch_results, top_n = 5)
  ggsave("model_selection_heatmap_fixed.png", heatmap, width = 12, height = 8, dpi = 300)
  
} else {
  cat("No batch_results found. Please run an analysis first or load saved results.\n")
  cat("\nTo load saved results:\n")
  cat('batch_results <- readRDS("path/to/your/all_results.rds")\n')
}

# Example 2: Loading from a specific batch directory
load_and_plot <- function(batch_dir) {
  # Load results
  all_results_file <- file.path(batch_dir, "all_results.rds")
  
  if (!file.exists(all_results_file)) {
    stop("Results file not found: ", all_results_file)
  }
  
  all_results <- readRDS(all_results_file)
  
  # Create batch_results object in expected format
  batch_results <- list(
    results = all_results,
    output_dir = batch_dir
  )
  
  # Check if bootstrap_n is saved
  if (file.exists(file.path(batch_dir, "analysis_metadata.txt"))) {
    metadata <- readLines(file.path(batch_dir, "analysis_metadata.txt"))
    boot_line <- grep("Bootstrap iterations:", metadata, value = TRUE)
    if (length(boot_line) > 0) {
      batch_results$bootstrap_n <- as.numeric(gsub(".*Bootstrap iterations: ", "", boot_line))
    }
  }
  
  # Create plots
  summary_fig <- create_phyloglm_summary_figure(
    batch_results,
    output_file = file.path(batch_dir, "summary_figure_panels_AB.png")
  )
  
  cat("Summary figure saved to:", file.path(batch_dir, "summary_figure_panels_AB.png"), "\n")
  
  return(batch_results)
}

# Usage example:
# batch_results <- load_and_plot("Outputs/PhyloglmResults/PhyloGLM_Batch_20250110_143022_boot100")