# Test script for PhyloGLM framework
# This runs the analysis step by step to help debug any issues

# Clear workspace
rm(list = ls())

# Load required libraries
library(ape)
library(phylolm)
library(ggplot2)
library(dplyr)

# Source framework components
cat("Loading framework components...\n")
source("phyloglm_framework/variable_classification.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/config_builder.R")

# Set paths - update these to match your files
data_file <- "Data_R_2025-06-09.csv"
tree_file <- "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
output_base_dir <- "Outputs/PhyloglmResults"

# Step 1: Load and check data
cat("\n=== Step 1: Loading Data ===\n")
data <- read.csv(data_file)
cat("Data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")

# Step 2: Load and check tree
cat("\n=== Step 2: Loading Tree ===\n")
tree <- read.nexus(tree_file)
cat("Tree has", length(tree$tip.label), "tips\n")

# Step 3: Create a simple test configuration
cat("\n=== Step 3: Creating Configuration ===\n")
config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop"),
  controls = c("logMass_AVONET"),
  complexity_levels = c("null", "main", "additive"),
  name = "Simple_Test"
)
print(config)

# Step 4: Check variable availability
cat("\n=== Step 4: Checking Variables ===\n")
all_vars <- c(config$response, config$predictors, config$controls)
for (var in all_vars) {
  if (var %in% colnames(data)) {
    cat(var, "- Found (", sum(!is.na(data[[var]])), "non-missing values)\n")
  } else {
    cat(var, "- NOT FOUND!\n")
  }
}

# Step 5: Try data preparation
cat("\n=== Step 5: Preparing Data ===\n")
prepared <- tryCatch({
  prepare_analysis_data(data, tree, config)
}, error = function(e) {
  cat("ERROR in data preparation:", e$message, "\n")
  NULL
})

if (!is.null(prepared)) {
  cat("Data preparation successful!\n")
  cat("Final dataset has", prepared$n_species, "species\n")
  
  # Step 6: Check data quality
  cat("\n=== Step 6: Checking Data Quality ===\n")
  quality <- check_data_quality(prepared)
  cat("Data quality acceptable:", quality$is_acceptable, "\n")
  if (!quality$is_acceptable) {
    cat("Issues:", paste(names(quality$issues), collapse = ", "), "\n")
  }
  
  # Step 7: Try running single analysis
  cat("\n=== Step 7: Running Analysis ===\n")
  cat("Creating output directory...\n")
  test_output_dir <- file.path(output_base_dir, "test_run")
  dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- tryCatch({
    run_single_analysis(
      config = config,
      data = data,
      tree = tree,
      output_dir = test_output_dir,
      bootstrap_n = 10,  # Small number for testing
      save_outputs = TRUE
    )
  }, error = function(e) {
    cat("ERROR in analysis:", e$message, "\n")
    NULL
  })
  
  if (!is.null(result) && result$success) {
    cat("\nAnalysis completed successfully!\n")
    cat("Results saved to:", test_output_dir, "\n")
    cat("\nBest model:", result$comparison$comparison$Model[1], "\n")
    cat("Best AIC:", round(result$comparison$comparison$AIC[1], 2), "\n")
  }
}

# Step 8: Test with more complex configuration
cat("\n=== Step 8: Testing Complex Configuration ===\n")
complex_config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  name = "Complex_Test"
)

# Validate configuration
validation <- validate_config(complex_config, data)
cat("Configuration valid:", validation$valid, "\n")
if (!validation$valid) {
  cat("Issues:", paste(validation$issues, collapse = "; "), "\n")
}
if (length(validation$warnings) > 0) {
  cat("Warnings:", paste(validation$warnings, collapse = "; "), "\n")
}

cat("\n=== Testing Complete ===\n")
cat("If all steps completed without errors, the framework is working correctly.\n")
cat("You can now run the full examples using example_usage.R\n")