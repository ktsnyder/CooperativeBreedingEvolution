# Simple test to find connection error

library(phylolm)
source('phyloglm_framework/batch_runner.R')
source('phyloglm_framework/config_builder.R')

configs <- create_standard_configs()
data <- read.csv('Data_R_2025-06-09.csv')
tree <- ape::read.nexus('2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

# Run just the first analysis with minimal bootstrap
cat("Running single analysis with bootstrap=10...\n")
result <- run_phyloglm_batch(
    analysis_configs = list(configs[[1]]),  # Just first config
    data = data,
    tree = tree,
    output_dir = "test_output_debug",
    bootstrap_n = 10,  # Very small
    save_intermediate = FALSE
)

cat("Success:", result[[1]]$success, "\n")
if (!result[[1]]$success) {
  cat("Error:", result[[1]]$error, "\n")
}