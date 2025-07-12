# Comprehensive Example of PhyloGLM Framework Usage
# This script demonstrates the complete workflow for multiple PhyloGLM analyses

# Load required libraries
library(ape)
library(phylolm)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source all framework components
source("phyloglm_framework/variable_classification.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/config_builder.R")
source("phyloglm_framework/visualization_framework.R")

# Set paths
data_file <- "Data_R_2025-06-09.csv"
tree_file <- "2024-05-26ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex"
output_base_dir <- "Outputs/PhyloGLM_Results"

#' Complete workflow example
#' Demonstrates all the requested analysis types
run_comprehensive_analysis <- function() {
  
  # Load data and tree
  cat("Loading data and tree...\n")
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Step 1: Examine available variables
  cat("\n=== Step 1: Variable Examination ===\n")
  classifications <- get_variable_classifications()
  
  # Analyze variables of interest
  vars_to_analyze <- c(
    # Response variables
    "FemaleSong_Agg01", "HighConfidence_Coop", "Griesser2017FamilialLiving",
    # Predictors
    "Territory", "TerritorialityWeakVsStrong",
    # Controls and moderators
    "logMass_AVONET", "PercentAbsLogWingDimorphism", "logMaleFemalePlumageDiffAbs",
    "Centroid.Latitude_AVONET", "Migration_AVONET", "GeographicRegion_Jetz"
  )
  
  var_summary <- summarize_variables(data, vars_to_analyze)
  print(var_summary[, c("variable", "type", "label", "n_complete", "missing_prop")])
  
  # Step 2: Create comprehensive analysis configurations
  cat("\n=== Step 2: Creating Analysis Configurations ===\n")
  
  configs <- list()
  
  # A) Analyses with three-state Territory variable
  configs$fs_cb_territory3 <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Territory"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway"),
    name = "FS_CB_Territory3State"
  )
  
  # B) Familial Living analyses (both directions)
  configs$fs_fl <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_FamilialLiving_Territorial"
  )
  
  configs$fl_fs <- create_analysis_config(
    response = "Griesser2017FamilialLiving",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    name = "FamilialLiving_FS_Territorial"
  )
  
  # C) Analyses with different continuous variables
  # Wing dimorphism
  configs$fs_cb_wingdim <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("PercentAbsLogWingDimorphism"),
    name = "FS_CB_WingDimorphism"
  )
  
  # Plumage dimorphism
  configs$fs_cb_plumdim <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("logMaleFemalePlumageDiffAbs"),
    name = "FS_CB_PlumageDimorphism"
  )
  
  # Absolute latitude
  configs$fs_cb_latitude <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET", "Centroid.Latitude_AVONET"),
    transformations = list(Centroid.Latitude_AVONET = "abs"),
    name = "FS_CB_AbsoluteLatitude"
  )
  
  # D) Analyses with additional factors
  # Geographic region
  configs$fs_cb_region <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "GeographicRegion_Jetz"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("main", "additive", "twoway"),
    name = "FS_CB_GeographicRegion"
  )
  
  # Migration
  configs$fs_cb_migration <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Migration_AVONET"),
    controls = c("logMass_AVONET"),
    name = "FS_CB_Migration"
  )
  
  # Complex model with multiple factors
  configs$fs_cb_complex <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Territory", "Migration_AVONET"),
    controls = c("logMass_AVONET", "PercentAbsLogWingDimorphism"),
    complexity_levels = c("main", "additive", "twoway"),
    max_interactions = 2,
    name = "FS_CB_Complex_MultiFactors"
  )
  
  # Save configurations to YAML for future use
  save_configs_to_yaml(configs, "phyloglm_comprehensive_analyses.yaml")
  
  # Step 3: Run batch analysis
  cat("\n=== Step 3: Running Batch Analysis ===\n")
  cat("This will run", length(configs), "separate analyses\n")
  
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir,
    parallel = FALSE,  # Set to TRUE for faster processing
    bootstrap_n = 100,
    save_intermediate = TRUE
  )
  
  # Step 4: Create comprehensive visualizations
  cat("\n=== Step 4: Creating Visualizations ===\n")
  viz_dir <- file.path(batch_results$output_dir, "comprehensive_visualizations")
  dir.create(viz_dir, showWarnings = FALSE)
  
  # Forest plot comparing cooperative breeding effects across analyses
  cb_forest <- create_comparative_forest_plot(
    batch_results,
    parameter_name = "HighConfidence_Coop",
    analysis_labels = c(
      "FS_CB_Territory3State" = "3-State Territory",
      "FS_CB_WingDimorphism" = "Wing Dimorphism Control",
      "FS_CB_PlumageDimorphism" = "Plumage Dimorphism Control",
      "FS_CB_AbsoluteLatitude" = "Latitude Control",
      "FS_CB_GeographicRegion" = "Geographic Region",
      "FS_CB_Migration" = "Migration Status",
      "FS_CB_Complex_MultiFactors" = "Complex Multi-Factor"
    )
  )
  ggsave(file.path(viz_dir, "cb_effects_comparison.png"), 
         cb_forest, width = 10, height = 8)
  
  # Model selection heatmap
  model_heatmap <- create_model_selection_heatmap(batch_results, top_n = 3)
  ggsave(file.path(viz_dir, "model_selection_heatmap.png"), 
         model_heatmap, width = 12, height = 8)
  
  # Effect size matrix
  effect_matrix <- create_effect_matrix_plot(
    batch_results,
    parameters = c("HighConfidence_Coop", "TerritorialityWeakVsStrong", 
                   "Territory1", "Territory2")
  )
  ggsave(file.path(viz_dir, "effect_size_matrix.png"), 
         effect_matrix, width = 10, height = 8)
  
  # Sample size comparison
  sample_plot <- create_sample_size_plot(batch_results)
  ggsave(file.path(viz_dir, "sample_sizes.png"), 
         sample_plot, width = 10, height = 6)
  
  # Create comprehensive summary figure
  summary_fig <- create_batch_summary_figure(
    batch_results,
    main_parameter = "HighConfidence_Coop",
    output_file = file.path(viz_dir, "comprehensive_summary.png")
  )
  
  # Step 5: Extract and compare specific results
  cat("\n=== Step 5: Extracting Key Results ===\n")
  
  # Compare effect of cooperative breeding across different control variables
  cb_effects <- data.frame()
  
  for (analysis_name in names(batch_results$results)) {
    result <- batch_results$results[[analysis_name]]
    
    if (result$success && !is.null(result$effects)) {
      # Get CB effect from best model
      best_model <- result$comparison$comparison$Model[1]
      cb_effect <- result$effects[
        result$effects$Model == best_model & 
        result$effects$Parameter == "HighConfidence_Coop",
      ]
      
      if (nrow(cb_effect) > 0) {
        cb_effect$Analysis <- analysis_name
        cb_effect$N_species <- result$prepared_data$n_species
        cb_effects <- rbind(cb_effects, cb_effect[1, ])
      }
    }
  }
  
  # Save comparison table
  write.csv(cb_effects, 
            file.path(viz_dir, "cooperative_breeding_effects_comparison.csv"),
            row.names = FALSE)
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  cat("Analyses completed:", sum(sapply(batch_results$results, function(r) r$success)), 
      "of", length(batch_results$results), "\n")
  
  # Show CB effects across analyses
  if (nrow(cb_effects) > 0) {
    cat("\nCooperative Breeding Effects (Odds Ratios):\n")
    cb_summary <- cb_effects[, c("Analysis", "OddsRatio", "p_value", "N_species")]
    cb_summary$Significant <- ifelse(cb_summary$p_value < 0.05, "*", "")
    print(cb_summary)
  }
  
  return(batch_results)
}

#' Run specific analysis scenario
#' Allows running individual analysis types
run_specific_scenario <- function(scenario = "territory3") {
  
  # Load data and tree
  data <- read.csv(data_file)
  tree <- read.nexus(tree_file)
  
  # Define scenario-specific configuration
  config <- switch(scenario,
    "territory3" = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "Territory"),
      controls = c("logMass_AVONET"),
      complexity_levels = c("null", "main", "additive", "twoway"),
      name = "Territory_3Level_Analysis"
    ),
    
    "familial" = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
      controls = c("logMass_AVONET"),
      complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
      name = "Familial_Living_Analysis"
    ),
    
    "dimorphism" = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
      controls = c("PercentAbsLogWingDimorphism", "logMaleFemalePlumageDiffAbs"),
      name = "Sexual_Dimorphism_Analysis"
    ),
    
    "geographic" = create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c("HighConfidence_Coop", "GeographicRegion_Jetz", "Migration_AVONET"),
      controls = c("logMass_AVONET"),
      name = "Geographic_Factors_Analysis"
    ),
    
    stop("Unknown scenario. Choose from: territory3, familial, dimorphism, geographic")
  )
  
  # Run single analysis
  result <- run_single_analysis(
    config = config,
    data = data,
    tree = tree,
    output_dir = file.path(output_base_dir, paste0("scenario_", scenario)),
    bootstrap_n = 100,
    save_outputs = TRUE
  )
  
  return(result)
}

# Main execution
if (interactive()) {
  cat("PhyloGLM Framework - Comprehensive Analysis\n")
  cat("==========================================\n")
  cat("Choose an option:\n")
  cat("1. Run complete comprehensive analysis (all scenarios)\n")
  cat("2. Run specific scenario\n")
  cat("3. Load and visualize previous results\n")
  cat("\nEnter choice (1-3): ")
  
  choice <- readline()
  
  if (choice == "1") {
    results <- run_comprehensive_analysis()
  } else if (choice == "2") {
    cat("\nAvailable scenarios:\n")
    cat("- territory3: Three-state territory analysis\n")
    cat("- familial: Familial living analysis\n")
    cat("- dimorphism: Sexual dimorphism controls\n")
    cat("- geographic: Geographic factors\n")
    cat("\nEnter scenario name: ")
    scenario <- readline()
    results <- run_specific_scenario(scenario)
  } else if (choice == "3") {
    cat("\nEnter path to previous results directory: ")
    results_dir <- readline()
    if (file.exists(file.path(results_dir, "all_results.rds"))) {
      results <- readRDS(file.path(results_dir, "all_results.rds"))
      cat("Results loaded successfully!\n")
    } else {
      stop("Results file not found")
    }
  } else {
    stop("Invalid choice")
  }
}