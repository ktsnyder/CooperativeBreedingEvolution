# Example: How to use the stepwise model expansion framework
# This example shows different ways to expand your phyloglm models

# Method 1: Simple interface (recommended for most users)
# ======================================================

source("phyloglm_framework/run_stepwise_simple.R")

# Basic usage - expands best models for FS->CB and CB->FS
results <- run_stepwise_simple()

# Include Familial Living (note: this reduces sample size)
results_with_fam <- run_stepwise_simple(include_familial = TRUE)

# Skip bootstrap for faster testing
results_quick <- run_stepwise_simple(bootstrap_n = 0)


# Method 2: Direct function calls (more control)
# ==============================================

source("phyloglm_framework/stepwise_model_expansion.R")
source("phyloglm_framework/stepwise_visualization.R")

# Load your existing results
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Expand a specific analysis
fs_cb_result <- all_results[["FS_vs_CB_TerrWS_Mass"]]

expanded <- stepwise_expand_model(
  base_result = fs_cb_result,
  data = data,
  tree = tree,
  aic_threshold = 2,      # Require AIC improvement > 2
  n_bootstrap = 100,      # Bootstrap iterations
  verbose = TRUE          # Print progress
)

# Print summary
print(expanded)

# Create visualizations
forest_plot <- plot_expanded_model_forest(expanded)
improvement_plot <- plot_stepwise_improvement(expanded)

# Save plots
ggsave("expanded_model_forest.png", forest_plot, width = 10, height = 8)
ggsave("stepwise_improvement.png", improvement_plot, width = 10, height = 10)


# Method 3: Custom candidate predictors
# =====================================

# Define your own set of predictors to test
my_candidates <- list(
  continuous = list(
    Territory_numeric = list(
      name = "Territory",
      type = "continuous",
      transform = function(x) as.numeric(x)
    )
  ),
  binary = list(
    GeographicRegion = list(
      name = "GeographicRegion_Jetz",
      type = "binary",
      transform = NULL
    )
  )
)

# Run with custom candidates
expanded_custom <- stepwise_expand_model(
  base_result = fs_cb_result,
  data = data,
  tree = tree,
  candidates = my_candidates,
  aic_threshold = 2
)


# Method 4: Batch processing multiple analyses
# ============================================

# Expand all analyses that include "TerrWS"
target_analyses <- grep("TerrWS", names(all_results), value = TRUE)

expansion_results <- run_stepwise_batch(
  results_file = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  target_analyses = target_analyses,
  output_dir = "Outputs/PhyloglmResults/stepwise_all_TerrWS"
)


# Understanding the output
# ========================

# The expansion result contains:
# - base_model: Original best model
# - final_model: Expanded model with added predictors
# - expansion_history: Step-by-step record of additions
# - improvement_tests: All predictors tested with their AIC improvements
# - bootstrap: Bootstrap results for final model (if requested)
# - effects: Effect sizes with confidence intervals

# Access specific components:
expanded$total_improvement  # Total AIC improvement
expanded$n_predictors_added # Number of predictors added
expanded$improvement_tests  # Table of all tests

# Get formula of final model
formula(expanded$final_model)