# Script to run phyloglm analyses with Territory_12vs3 as numeric variable

# Load required libraries and functions
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/model_fitting.R")
source("phyloglm_framework/model_comparison.R")
source("phyloglm_framework/effect_calculation.R")
source("phyloglm_framework/bootstrap_helpers.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/visualization_framework.R")

library(phylolm)
library(dplyr)

# Load the processed data
cat("Loading data...\n")
data <- readRDS("ProcessedData/phylogenetic_data_processed.rds")

# Load phylogenetic tree
tree <- readRDS("ProcessedData/matched_tree.rds")

# Create analysis configuration for FS_vs_CB_Terr3cont_Mass
config_fs_cb <- list(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "Territory_12vs3"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  numeric_vars = c("Territory_12vs3"),  # Treat Territory_12vs3 as numeric
  name = "FS_vs_CB_Terr3cont_Mass"
)

# Create analysis configuration for CB_vs_FS_Terr3cont_Mass
config_cb_fs <- list(
  response = "HighConfidence_Coop",
  predictors = c("FemaleSong_Agg01", "Territory_12vs3"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  numeric_vars = c("Territory_12vs3"),  # Treat Territory_12vs3 as numeric
  name = "CB_vs_FS_Terr3cont_Mass"
)

# Function to run a single analysis with numeric territory
run_numeric_territory_analysis <- function(config, data, tree, bootstrap_n = 1000) {
  
  cat("\n========================================\n")
  cat("Running analysis:", config$name, "\n")
  cat("========================================\n")
  
  # Prepare data - ensure Territory_12vs3 is numeric
  prepared <- prepare_phyloglm_data(data, config)
  
  # Convert Territory_12vs3 to numeric if it's not already
  if ("Territory_12vs3" %in% names(prepared$data)) {
    prepared$data$Territory_12vs3 <- as.numeric(as.character(prepared$data$Territory_12vs3))
    cat("Converted Territory_12vs3 to numeric. Values:", 
        sort(unique(prepared$data$Territory_12vs3)), "\n")
  }
  
  # Build formulas
  formulas <- build_formulas(config)
  cat("Built", length(formulas), "model formulas\n")
  
  # Fit models
  models <- fit_phyloglm_models(formulas, prepared$data, prepared$tree)
  cat("Successfully fit", sum(sapply(models$models, function(x) !is.null(x))), "models\n")
  
  # Compare models
  comparison <- compare_phyloglm_models(models)
  cat("\nModel comparison:\n")
  print(head(comparison$comparison, 10))
  
  # Calculate coefficients
  coefficients <- extract_all_coefficients(models, formulas)
  
  # Calculate effects for all models
  effects <- list()
  for (model_name in names(models$models)) {
    if (!is.null(models$models[[model_name]])) {
      effects[[model_name]] <- calculate_effect_sizes(
        models$models[[model_name]], 
        model_name
      )
    }
  }
  
  # Bootstrap if requested
  bootstrap_results <- NULL
  if (bootstrap_n > 0) {
    cat("\nRunning bootstrap with", bootstrap_n, "iterations...\n")
    bootstrap_results <- run_bootstrap_analysis(
      formulas = formulas,
      data = prepared$data,
      tree = prepared$tree,
      n_bootstrap = bootstrap_n,
      parallel = FALSE
    )
    
    # Update coefficients with bootstrap CIs
    coefficients <- add_bootstrap_ci(coefficients, bootstrap_results)
  }
  
  # Create output directory
  output_dir <- file.path("Outputs/PhyloglmResults/numeric_territory", config$name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save all results
  results <- list(
    success = TRUE,
    config = config,
    prepared_data = prepared,
    models = models,
    comparison = comparison,
    coefficients = coefficients,
    effects = effects,
    bootstrap = bootstrap_results,
    output_dir = output_dir
  )
  
  saveRDS(results, file.path(output_dir, "complete_analysis_results.rds"))
  
  # Save CSV summaries
  write.csv(comparison$comparison, 
            file.path(output_dir, "model_comparison.csv"), 
            row.names = FALSE)
  write.csv(coefficients, 
            file.path(output_dir, "coefficients.csv"), 
            row.names = FALSE)
  
  # Create visualizations
  cat("\nCreating visualizations...\n")
  visualize_phyloglm_results(results, output_prefix = file.path(output_dir, config$name))
  
  return(results)
}

# Run the analyses
cat("Starting analyses with Territory as numeric variable...\n\n")

# Run FS_vs_CB_Terr3cont_Mass
results_fs_cb <- run_numeric_territory_analysis(config_fs_cb, data, tree, bootstrap_n = 1000)

# Run CB_vs_FS_Terr3cont_Mass  
results_cb_fs <- run_numeric_territory_analysis(config_cb_fs, data, tree, bootstrap_n = 1000)

cat("\n\nAnalyses complete! Results saved to:\n")
cat("- Outputs/PhyloglmResults/numeric_territory/FS_vs_CB_Terr3cont_Mass/\n")
cat("- Outputs/PhyloglmResults/numeric_territory/CB_vs_FS_Terr3cont_Mass/\n")

# Create a combined summary plot
cat("\nCreating combined summary plot...\n")
source("phyloglm_framework/phyloglm_summary_plots_simple.R")

combined_results <- list(
  FS_vs_CB_Terr3cont_Mass = results_fs_cb,
  CB_vs_FS_Terr3cont_Mass = results_cb_fs
)

# Create summary figure
create_simple_summary_figure(
  combined_results, 
  output_file = "Outputs/PhyloglmResults/numeric_territory/territory_numeric_summary"
)