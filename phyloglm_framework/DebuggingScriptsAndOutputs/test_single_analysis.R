# Test single analysis to debug the connection error

# Load required packages and functions
library(phylolm)
library(ape)
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/config_builder.R")

# Load data
data <- read.csv("Data_R_2025-06-09.csv")
tree <- ape::read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Create a single config
config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  complexity_levels = c("null", "main", "additive"),
  max_interactions = 2,
  name = "Test_FS_CB_TerrWS_Mass"
)

# Try running without saving outputs
cat("Testing without save_outputs...\n")
result1 <- tryCatch({
  run_single_analysis(
    config = config,
    data = data,
    tree = tree,
    output_dir = tempdir(),
    bootstrap_n = 0,  # No bootstrap for quick test
    save_outputs = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
  list(success = FALSE, error = e$message)
})

cat("\nResult:", result1$success, "\n")
if (!result1$success) {
  cat("Error details:", result1$error, "\n")
}

# If that works, try with save_outputs
if (result1$success) {
  cat("\nTesting with save_outputs...\n")
  result2 <- tryCatch({
    run_single_analysis(
      config = config,
      data = data,
      tree = tree,
      output_dir = "test_output",
      bootstrap_n = 0,
      save_outputs = TRUE
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    list(success = FALSE, error = e$message)
  })
  
  cat("Result:", result2$success, "\n")
}