# PhyloGLM runner
# Kate Snyder
# 

data = read.csv("Data_R.csv")
tree = read.nexus(treefile)

source(file.path("PhyloGLM_functions", "batch_runner.R"))
source(file.path("PhyloGLM_functions", "config_builder.R"))

if (!exists("nBoot")) {
  nBoot = 100
  print("Defaulting to nBoot = 100. To run the tests as performed in the manuscript, set nBoot = 500.")
}

## Base analyses ----
## Results akin to Table 1 are in Outputs/PhyloGLM_outputs/PhyloGLM_Batch_[Date]_[Time]_boot100/individual_analyses/[analysis name]/best_model_effects.csv
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
  bootstrap_n = nBoot,  
  save_intermediate = TRUE
)

## Base analyses using alternative cooperative breeding classifications ----
# Generates results akin to Supplemental Tables 15 & 16
if (!exists("run_alt_coop")) {
  run_alt_coop = FALSE
  print("Defaulting to not running base phyloglm analyses with each alternative cooperative breeding classification method. To run these analyses, define run_alt_coop = TRUE")
}

if (run_alt_coop) {
  source("subsettreedata.R")
  
  AltCoops = c("MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnCoop", "CockburnInferred", "JetzCoopInclCockburn", "Griesser2017Coop", "DaleCoop", "CornwallisCoop")
  
  output_dir = file.path("Outputs", "PhyloGLM_outputs", "AltCoops")
  
  for (i in 1:length(AltCoops)) {
    tempcoop = AltCoops[i]
    
    subsetout <- subsettreedata(columns = c(tempcoop, "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong"), newdata = "Data_R.csv", newtree = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
    subsetdf = subsetout$subsetdf
    rownames(subsetdf) <- subsetdf$species
    subsettree = subsetout$subsettree
    
    FS_response_formula = paste("FemaleSong_Agg01 ~", tempcoop ,"* TerritorialityWeakVsStrong + logMass_normalized")
    print(FS_response_formula)
    
    FS_vs_CBxTerrWS_Mass <- run_bootstrap_model(formula = as.formula(FS_response_formula), data = subsetdf, tree = subsettree, n_boot = nBoot, method = "logistic_MPLE", save_prefix = paste0("FS_vs_", tempcoop, "xTerrWS_Mass"), save_matrices = TRUE, matrix_dir = output_dir, save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)
    
    CB_response_formula = paste(tempcoop, "~ FemaleSong_Agg01 * TerritorialityWeakVsStrong")
    print(CB_response_formula)
    
    CB_vs_FSxTerrWS <- run_bootstrap_model(formula = as.formula(CB_response_formula), data = subsetdf, tree = subsettree, n_boot = nBoot, method = "logistic_MPLE", save_prefix = paste0(tempcoop,"_vs_FSxTerrWS"), save_matrices = TRUE, matrix_dir = output_dir, save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)
  }
  
  source(file.path("PhyloGLM_functions", "compile_phyloglm_altcoops_tables.R"))
} else {
  print("Not running base phyloglm analyses with each alternative cooperative breeding classification method. To run these analyses, define run_alt_coop = TRUE")
}

## Direct-comparison of the effect of cooperative breeding versus other sociality variables on prediction of female song  ----
# Performs analyses and generates tables with results akin to those in Table 2, Supplemental Table 18
# Results located in Outputs/PhyloGLM_outputs/Replace_Coop_With_Soc/Replace_Coop_Additive/Replace_Coop_Summary_Table_boot[nBoot].csv 
source(file.path("PhyloGLM_functions", "replace_Coop_phyloglms.R"))


## Stepwise iterative model expansion ----
# Performs forward stepwise expansion starting from the best base models
# Results located in Outputs/PhyloGLM_outputs/Stepwise/stepwise_[timestamp]/stepwise_summary_boot[nBoot].csv
if (!exists("run_stepwise")) {
  run_stepwise = FALSE
  print("Defaulting to not running stepwise expansion. To run stepwise expansion, define run_stepwise = TRUE")
}

if (run_stepwise) {
  source(file.path("PhyloGLM_functions", "stepwise_expansion.R"))
  
  # Define starting formulas based on the best base models
  starting_formulas <- list(
    # Female song as response with territoriality weak/strong
    FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized,
    # Female song as response with territoriality year-round  
    FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized,
    # Cooperative breeding as response with territoriality weak/strong
    HighConfidence_Coop ~ FemaleSong_Agg01 * TerritorialityWeakVsStrong,
    # Cooperative breeding as response with territoriality year-round
    HighConfidence_Coop ~ FemaleSong_Agg01 * Territory_12vs3
  )
  
  stepwise_results <- run_stepwise_expansion(
    formulas = starting_formulas,
    data_path = "Data_R.csv",
    tree_path = treefile,
    aic_threshold = 2,
    n_bootstrap = nBoot,
    max_iterations = 100,
    output_dir = file.path("Outputs", "PhyloGLM_outputs"),
    verbose = TRUE
  )
} else {
  print("Not running stepwise expansion. To run stepwise expansion, define run_stepwise = TRUE")
}


