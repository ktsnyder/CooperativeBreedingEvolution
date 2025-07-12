# Debug the batch runner error

# Add detailed error tracking
options(error = function() {
  calls <- sys.calls()
  cat("\n=== ERROR TRACEBACK ===\n")
  for(i in seq_along(calls)) {
    cat(i, ": ", deparse(calls[[i]])[1], "\n", sep = "")
  }
  cat("======================\n")
})

# Load everything
library(phylolm)
library(ape)
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/config_builder.R")

# Create configs and run single analysis with detailed output
configs <- create_standard_configs()
config <- configs$FS_vs_CB_TerrWS_absLat

# Manually run the steps from run_single_analysis
data <- read.csv("Data_R_2025-06-09.csv")
tree <- ape::read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

cat("Step 1: Data preparation\n")
tryCatch({
  prepared <- prepare_analysis_data(
    data = data,
    tree = tree,
    config = config,
    verbose = TRUE
  )
  cat("  Success - prepared data has", nrow(prepared$data), "species\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  stop("Failed at data preparation")
})

cat("\nStep 2: Quality checks\n")
tryCatch({
  quality <- check_phyloglm_data_quality(prepared$data, config)
  cat("  Success\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  stop("Failed at quality check")
})

cat("\nStep 3: Build formulas\n")
tryCatch({
  variable_types <- detect_variable_types(prepared$data)
  
  formulas <- build_phyloglm_formulas(
    response_var = config$response,
    predictor_vars = config$predictors,
    control_vars = config$controls,
    variable_types = variable_types,
    complexity_levels = config$complexity_levels,
    max_interactions = config$max_interactions
  )
  cat("  Success - built", length(formulas), "formulas\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  stop("Failed at formula building")
})

cat("\nStep 4: Fit models (first one only)\n")
tryCatch({
  first_formula <- formulas[[1]]
  cat("  Trying formula:", deparse(first_formula), "\n")
  
  model <- phyloglm(
    formula = first_formula,
    data = prepared$data,
    phy = prepared$tree,
    method = "logistic_MPLE",
    btol = 30,
    log.alpha.bound = 4
  )
  cat("  Success\n")
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  cat("  Error class:", class(e), "\n")
})