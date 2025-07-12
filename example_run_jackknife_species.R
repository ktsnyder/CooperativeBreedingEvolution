# Example script for running jackknife species analysis
# This script demonstrates how to use run_jackknife_species_phylopath()

# Load required libraries
library(phylopath)
library(phytools)
library(dplyr)

# Source the functions
source("run_phylopath_fxns.R")
source("run_jackknife_species_phylopath.R")

# Load data
df <- read.csv("Data_R_2025-06-09.csv", stringsAsFactors = FALSE)
tree <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Run jackknife analysis
# This will take a long time as it runs 875 iterations (one for each species)
jackknife_results <- run_jackknife_species_phylopath(
  dfIn = df,
  tree = tree,
  female_song_var = "FemaleSong_Agg01",
  coop_breeding_var = "HighConfidence_Coop", 
  territoriality_var = "TerritorialityWeakVsStrong",
  mass_var = "logMass_AVONET",
  save_conditional_plots = FALSE,  # Set to TRUE if you want plots (will create large PDF)
  save_path_coefficients = TRUE,
  output_dir = file.path("Outputs", "PhylopathJackknife")
)

# The function will save several files:
# 1. results_JackknifeSpecies_n875_[date].csv - Main results with best models
# 2. detailed_models_JackknifeSpecies_n875_[date].csv - All models with delta_CICc < 2
# 3. model_frequencies_JackknifeSpecies_n875_[date].csv - How often each model appears
# 4. edge_summary_JackknifeSpecies_n875_[date].csv - Summary of path coefficients

# You can also access the results directly:
# jackknife_results$results_df - Main results data frame
# jackknife_results$model_frequencies - Model frequency summary
# jackknife_results$detailed_models - Detailed model information
# jackknife_results$edge_summary - Edge coefficient summary

# Example: Find species whose removal most affects the results
results_summary <- jackknife_results$results_df %>%
  mutate(
    # Check if best model changed from the most common model
    model_changed = best_model != names(which.max(table(jackknife_results$results_df$best_model)))
  ) %>%
  filter(model_changed) %>%
  select(RemovedSpecies, best_model, best_CICc, nSub2dCIC_Models)

print("Species whose removal changed the best model:")
print(results_summary)