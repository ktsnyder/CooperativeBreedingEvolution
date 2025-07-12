# Run ALL PhyloGLM Analyses - Complete Set
# This script runs all standard analyses including dimorphism and geographic variables

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
output_dir <- "Outputs/PhyloglmResults"

# Load data and tree
cat("Loading data and tree...\n")
data <- read.csv(data_file)
tree <- read.nexus(tree_file)

# Create comprehensive configuration list
cat("\nCreating comprehensive analysis configurations...\n")

configs <- list()

# 1. BASIC ANALYSES (with logMass control)
# FS vs CB with territoriality
configs$fs_cb_terr_mass <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  name = "FS_CB_TerrWS_Mass"
)

# CB vs FS with territoriality (reverse)
configs$cb_fs_terr_mass <- create_analysis_config(
  response = "HighConfidence_Coop",
  predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  name = "CB_FS_TerrWS_Mass"
)

# REMOVED 3-state territory analysis as requested

# 2. FAMILIAL LIVING ANALYSES
# FS vs Familial Living with territoriality
configs$fs_fl_terr_mass <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  name = "FS_FL_TerrWS_Mass"
)

# FL vs FS with territoriality (reverse)
configs$fl_fs_terr_mass <- create_analysis_config(
  response = "Griesser2017FamilialLiving",
  predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  name = "FL_FS_TerrWS_Mass"
)

# 3. WING DIMORPHISM ANALYSES
# FS vs CB with wing dimorphism control
configs$fs_cb_terr_wingdim <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("PercentAbsLogWingDimorphism"),
  name = "FS_CB_TerrWS_WingDim"
)

# REMOVED 3-state territory analysis as requested

# 4. PLUMAGE DIMORPHISM ANALYSES
# FS vs CB with plumage dimorphism control
configs$fs_cb_terr_plumdim <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMaleFemalePlumageDiffAbs"),
  name = "FS_CB_TerrWS_PlumDim"
)

# REMOVED 3-state territory analysis as requested

# 5. LATITUDE ANALYSES
# Create configurations that include latitude
configs$fs_cb_terr_lat <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  name = "FS_CB_TerrWS_Lat"
)
# Note: This needs special handling to add latitude as additional control

# 6. MIGRATION ANALYSES
# FS vs CB with migration
configs$fs_cb_migr_mass <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "Migration_AVONET"),
  controls = c("logMass_AVONET"),
  name = "FS_CB_Migration_Mass"
)

# 7. GEOGRAPHIC REGION ANALYSES
# FS vs CB with geographic region
configs$fs_cb_region_mass <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "GeographicRegion_Jetz"),
  controls = c("logMass_AVONET"),
  name = "FS_CB_Region_Mass"
)

# Print summary of analyses to run
cat("\nTotal analyses configured:", length(configs), "\n")
cat("Analysis names:\n")
for (name in names(configs)) {
  cat("  -", configs[[name]]$name, "\n")
}

# Ask for confirmation
cat("\nThis will run", length(configs), "analyses, each with 15 models.\n")
cat("Total models to fit:", length(configs) * 15, "\n")
cat("Proceed? (y/n): ")
response <- readline()

if (tolower(response) == "y") {
  # Run batch analysis
  cat("\nStarting batch analysis...\n")
  
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_dir,
    bootstrap_n = 100,  # Adjust as needed
    parallel = FALSE,   # Set to TRUE for faster processing
    save_intermediate = TRUE
  )
  
  # Create comprehensive visualizations
  cat("\nCreating comprehensive visualizations...\n")
  viz_dir <- file.path(batch_results$output_dir, "comprehensive_visualizations")
  dir.create(viz_dir, showWarnings = FALSE)
  
  # Forest plot for cooperative breeding effects
  if ("HighConfidence_Coop" %in% unlist(lapply(configs, function(c) c$predictors))) {
    cb_forest <- create_comparative_forest_plot(
      batch_results,
      parameter_name = "HighConfidence_Coop"
    )
    ggsave(file.path(viz_dir, "cooperative_breeding_effects.png"), 
           cb_forest, width = 12, height = 10)
  }
  
  # Model selection heatmap
  model_heatmap <- create_model_selection_heatmap(batch_results, top_n = 5)
  ggsave(file.path(viz_dir, "model_selection_all_analyses.png"), 
         model_heatmap, width = 14, height = 10)
  
  # Sample size comparison
  sample_plot <- create_sample_size_plot(batch_results)
  ggsave(file.path(viz_dir, "sample_sizes_all_analyses.png"), 
         sample_plot, width = 12, height = 8)
  
  # Create summary figure
  summary_fig <- create_batch_summary_figure(
    batch_results,
    main_parameter = "HighConfidence_Coop",
    output_file = file.path(viz_dir, "all_analyses_summary.png")
  )
  
  cat("\nAnalysis complete! Results saved to:", batch_results$output_dir, "\n")
  
} else {
  cat("Analysis cancelled.\n")
}