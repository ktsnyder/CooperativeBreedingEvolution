# Test script to verify bootstrap tracking is working

# Load required libraries
library(ape)
library(phylolm)
library(ggplot2)
library(dplyr)

# Source framework components
source("phyloglm_framework/variable_classification.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/config_builder.R")
source("phyloglm_framework/visualization_framework.R")

# Set paths
data_file <- "Data_R_2025-06-09.csv"
tree_file <- "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
output_base_dir <- "Outputs/PhyloglmResults"

# Load data and tree
cat("Loading data and tree...\n")
data <- read.csv(data_file)
tree <- read.nexus(tree_file)

# Test 1: Single analysis with custom bootstrap count
cat("\n=== Test 1: Single Analysis with 50 Bootstrap ===\n")
config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop"),
  controls = c("logMass_AVONET"),
  name = "Test_Bootstrap_50"
)

result <- run_single_analysis(
  config = config,
  data = data,
  tree = tree,
  output_dir = file.path(output_base_dir, "test_bootstrap_50"),
  bootstrap_n = 50,
  save_outputs = TRUE
)

# Check if bootstrap count is in saved files
if (result$success) {
  cat("\nChecking saved files for bootstrap info...\n")
  
  # Check model comparison CSV
  model_comp <- read.csv(file.path(result$output_dir, "model_comparison.csv"))
  if ("bootstrap_n" %in% colnames(model_comp)) {
    cat("✓ Bootstrap count found in model_comparison.csv:", unique(model_comp$bootstrap_n), "\n")
  } else {
    cat("✗ Bootstrap count NOT found in model_comparison.csv\n")
  }
  
  # Check complete results RDS
  complete <- readRDS(file.path(result$output_dir, "complete_analysis_results.rds"))
  if (!is.null(complete$bootstrap_n)) {
    cat("✓ Bootstrap count found in complete_analysis_results.rds:", complete$bootstrap_n, "\n")
  } else {
    cat("✗ Bootstrap count NOT found in complete_analysis_results.rds\n")
  }
}

# Test 2: Batch analysis with different bootstrap count
cat("\n=== Test 2: Batch Analysis with 200 Bootstrap ===\n")
configs <- list(
  test1 = create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop"),
    controls = c("logMass_AVONET"),
    name = "Batch_Test_1"
  ),
  test2 = create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01"),
    controls = c("logMass_AVONET"),
    name = "Batch_Test_2"
  )
)

batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  output_dir = output_base_dir,
  bootstrap_n = 200,
  parallel = FALSE,
  save_intermediate = TRUE
)

# Check batch directory naming
cat("\nBatch output directory:", batch_results$output_dir, "\n")
if (grepl("boot200", batch_results$output_dir)) {
  cat("✓ Bootstrap count included in directory name\n")
} else {
  cat("✗ Bootstrap count NOT in directory name\n")
}

# Check metadata file
metadata_file <- file.path(batch_results$output_dir, "analysis_metadata.txt")
if (file.exists(metadata_file)) {
  cat("\n✓ Metadata file exists. Contents:\n")
  cat(readLines(metadata_file), sep = "\n")
} else {
  cat("\n✗ Metadata file NOT found\n")
}

# Test visualization with bootstrap info
cat("\n=== Test 3: Visualization with Bootstrap Info ===\n")
if (length(batch_results$results) > 0) {
  # Create a forest plot
  p <- create_comparative_forest_plot(batch_results, "HighConfidence_Coop")
  
  # Check if subtitle contains bootstrap info
  build <- ggplot_build(p)
  labels <- p$labels
  if (!is.null(labels$subtitle) && grepl("bootstrap", labels$subtitle)) {
    cat("✓ Bootstrap info found in plot subtitle\n")
  } else {
    cat("✗ Bootstrap info NOT in plot subtitle\n")
  }
  
  # Save plot
  ggsave(file.path(batch_results$output_dir, "test_forest_plot.png"), p, width = 8, height = 6)
}

cat("\n=== Bootstrap Tracking Test Complete ===\n")
cat("Check the output directories to verify all bootstrap information is saved correctly.\n")