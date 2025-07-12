# Example Usage of PhyloGLM Framework
# This script demonstrates how to use the flexible phyloglm analysis framework

# Load required libraries and source files
library(ape)
library(phylolm)
library(ggplot2)
library(dplyr)

# Source framework components
source("phyloglm_framework/variable_classification.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/config_builder.R")

# Set paths
data_file <- "Data_R_2025-06-09.csv"
tree_file <- "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
output_base_dir <- "Outputs/PhyloglmResults"

#' Example 1: Single Analysis
#' Analyze female song vs cooperative breeding with territorial moderation
run_single_analysis_example <- function() {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Create configuration
  config <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FemaleSong_CB_Territorial"
  )
  
  # Validate configuration
  validation <- validate_config(config, data)
  if (!validation$valid) {
    stop("Configuration issues:", paste(validation$issues, collapse = "; "))
  }
  
  # Run analysis
  result <- run_single_analysis(
    config = config,
    data = data,
    tree = tree,
    output_dir = file.path(output_base_dir, "single_example"),
    bootstrap_n = 100,
    save_outputs = TRUE
  )
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  cat("N species:", result$prepared_data$n_species, "\n")
  cat("Best model:", result$comparison$comparison$Model[1], "\n")
  cat("Best AIC:", round(result$comparison$comparison$AIC[1], 2), "\n")
  
  return(result)
}

#' Example 2: Batch Analysis with Standard Configurations
#' Run multiple standard analyses at once
run_batch_standard_example <- function() {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Create standard configurations
  configs <- create_standard_configs(
    include_sets = c("basic", "territorial")
  )
  
  # Run batch analysis
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir,
    parallel = FALSE,  # Set to TRUE for faster processing
    bootstrap_n = 100,
    save_intermediate = TRUE
  )
  
  return(batch_results)
}

#' Example 3: Custom Batch Analysis
#' Analyze effects of different cooperative breeding definitions
run_custom_batch_example <- function() {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Create custom configurations for different CB definitions
  cb_vars <- c("HighConfidence_Coop", "AnyCoopEqualsCoop", 
               "MeanCoopTie2Coop", "Griesser2017FamilialLiving")
  
  configs <- list()
  
  # Create bidirectional analyses for each CB variable
  for (cb_var in cb_vars) {
    # FS -> CB direction
    configs[[paste0("fs_", cb_var)]] <- create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c(cb_var, "TerritorialityWeakVsStrong"),
      controls = c("logMass_AVONET"),
      name = paste0("FS_vs_", cb_var)
    )
    
    # CB -> FS direction  
    configs[[paste0(cb_var, "_fs")]] <- create_analysis_config(
      response = cb_var,
      predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
      controls = c("logMass_AVONET"),
      name = paste0(cb_var, "_vs_FS")
    )
  }
  
  # Run batch
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir,
    parallel = FALSE,
    bootstrap_n = 100
  )
  
  return(batch_results)
}

#' Example 4: Using YAML Configuration
#' Load analyses from a configuration file
run_yaml_config_example <- function() {
  
  # First create a template if it doesn't exist
  config_file <- "phyloglm_analyses.yaml"
  if (!file.exists(config_file)) {
    create_config_template(config_file)
    stop("Template created at ", config_file, ". Please edit it and run again.")
  }
  
  # Load configurations from YAML
  configs <- load_configs_from_yaml(config_file)
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Validate all configurations
  for (i in seq_along(configs)) {
    validation <- validate_config(configs[[i]], data)
    if (!validation$valid) {
      stop(paste("Config", i, "has issues:", 
                 paste(validation$issues, collapse = "; ")))
    }
  }
  
  # Run batch
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir
  )
  
  return(batch_results)
}

#' Example 5: Variable Type Analysis
#' Examine different variable types in the dataset
examine_variables_example <- function() {
  
  # Load data
  data <- read.csv(data_file)
  
  # Get classifications
  classifications <- get_variable_classifications()
  
  # Analyze all classified variables
  all_vars <- unique(unlist(classifications))
  var_summary <- summarize_variables(data, all_vars)
  
  # Print summary by type
  cat("\n=== Variable Summary by Type ===\n")
  for (type in unique(var_summary$type)) {
    cat("\n", toupper(type), "VARIABLES:\n")
    type_vars <- var_summary[var_summary$type == type, ]
    print(type_vars[, c("variable", "label", "n_complete", "missing_prop")])
  }
  
  # Save full summary
  write.csv(var_summary, "variable_summary.csv", row.names = FALSE)
  
  return(var_summary)
}

#' Example 6: Complex Analysis with Multiple Variable Types
#' Demonstrate handling of categorical predictors
run_complex_variable_example <- function() {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Configuration with 3-level categorical variable
  config <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Territory", "Migration_AVONET"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway"),
    max_interactions = 2,
    name = "Complex_Categorical_Analysis"
  )
  
  # Run analysis
  result <- run_single_analysis(
    config = config,
    data = data,
    tree = tree,
    output_dir = file.path(output_base_dir, "complex_example"),
    bootstrap_n = 50  # Fewer for complex models
  )
  
  return(result)
}

#' Example 7: Geographic and Dimorphism Analysis
#' Include continuous moderators
run_geographic_dimorphism_example <- function() {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Create configurations
  configs <- list(
    # Wing dimorphism as control
    wing_dim = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
      controls = c("PercentAbsLogWingDimorphism"),
      name = "FS_CB_WingDimorphism"
    ),
    
    # Absolute latitude
    latitude = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
      controls = c("logMass_AVONET", "Centroid.Latitude_AVONET"),
      transformations = list(Centroid.Latitude_AVONET = "abs"),
      name = "FS_CB_AbsLatitude"
    ),
    
    # Geographic region as predictor
    region = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "GeographicRegion_Jetz"),
      controls = c("logMass_AVONET"),
      complexity_levels = c("main", "additive", "twoway"),
      name = "FS_CB_GeographicRegion"
    )
  )
  
  # Run batch
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir
  )
  
  return(batch_results)
}

# Main function to run all examples
if (interactive()) {
  cat("PhyloGLM Framework Examples\n")
  cat("===========================\n")
  cat("Choose an example to run:\n")
  cat("1. Single analysis example\n")
  cat("2. Batch analysis with standard configs\n")
  cat("3. Custom batch analysis\n")
  cat("4. YAML configuration example\n")
  cat("5. Variable examination\n")
  cat("6. Complex categorical analysis\n")
  cat("7. Geographic and dimorphism analysis\n")
  cat("\nEnter number (1-7): ")
  
  choice <- readline()
  
  result <- switch(choice,
    "1" = run_single_analysis_example(),
    "2" = run_batch_standard_example(),
    "3" = run_custom_batch_example(),
    "4" = run_yaml_config_example(),
    "5" = examine_variables_example(),
    "6" = run_complex_variable_example(),
    "7" = run_geographic_dimorphism_example(),
    stop("Invalid choice")
  )
}