# Simpler script using batch runner for numeric territory analyses

source("phyloglm_framework/batch_runner.R")

# Create configurations for numeric territory analyses
configs <- list()

# FS_vs_CB_Terr3cont_Mass
configs$FS_vs_CB_Terr3cont_Mass <- list(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "Territory_12vs3"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  numeric_vars = c("Territory_12vs3"),  # This tells the system to treat it as numeric
  name = "FS_vs_CB_Terr3cont_Mass"
)

# CB_vs_FS_Terr3cont_Mass
configs$CB_vs_FS_Terr3cont_Mass <- list(
  response = "HighConfidence_Coop", 
  predictors = c("FemaleSong_Agg01", "Territory_12vs3"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  numeric_vars = c("Territory_12vs3"),  # This tells the system to treat it as numeric
  name = "CB_vs_FS_Terr3cont_Mass"
)

# Run the batch analysis
results <- run_phyloglm_batch(
  configs = configs,
  data_path = "ProcessedData/phylogenetic_data_processed.rds",
  tree_path = "ProcessedData/matched_tree.rds",
  output_dir = "Outputs/PhyloglmResults/numeric_territory",
  n_bootstrap = 1000,
  parallel = FALSE,
  save_intermediate = TRUE
)

cat("\nAnalyses complete!\n")