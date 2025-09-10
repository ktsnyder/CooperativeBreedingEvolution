# PhyloGLM runner
# Kate Snyder
# 

data = read.csv("Data_R.csv")
tree = read.nexus(treefile)

source(file.path("PhyloGLM_functions", "batch_runner.R"))
source(file.path("PhyloGLM_functions", "config_builder.R"))

## Base analyses ----
## Results akin to Table 2 are in Outputs/PhyloGLM_outputs/PhyloGLM_Batch_[Date]_[Time]_boot100/individual_analyses/[analysis name]/best_model_effects.csv
configs <- list()
# 1. BASIC ANALYSES (with logMass_normalized control)
# FS vs CB with territoriality weak strong
configs$fs_cb_terrws_massnorm <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_normalized"),
  name = "FS_CB_TerrWS_MassNorm"
)

# CB vs FS with territoriality weak strong (reverse)
configs$cb_fs_terrws_massnorm <- create_analysis_config(
  response = "HighConfidence_Coop",
  predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
  controls = c("logMass_normalized"),
  name = "CB_FS_TerrWS_MassNorm"
)

# FS vs CB with territoriality year round
configs$fs_cb_terr3_massnorm <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "Territory_12vs3"),
  controls = c("logMass_normalized"),
  name = "FS_CB_Terr3_MassNorm"
)

# CB vs FS with territoriality year round (reverse)
configs$cb_fs_terr3_massnorm <- create_analysis_config(
  response = "HighConfidence_Coop",
  predictors = c("FemaleSong_Agg01", "Territory_12vs3"),
  controls = c("logMass_normalized"),
  name = "CB_FS_Terr3_MassNorm"
)

batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  output_dir = file.path("Outputs", "PhyloGLM_outputs"),
  bootstrap_n = 100,  
  save_intermediate = TRUE
)

## Direct-comparison of cooperative breeding versus other sociality variables ----



## Stepwise iterative model expansion ----


