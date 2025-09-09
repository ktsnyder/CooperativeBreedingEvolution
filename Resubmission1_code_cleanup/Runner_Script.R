#### Full Runner Script ----
### First revision, 9/1/2025

# Generates Supplemental Table 3 - phylANOVA
source("run_01_phylANOVA.R")

# Generates components used in: Figure 1; Extended Data Figures 1, 2; Supplemental Tables 4, 5, 6
# Performs only 20 simulations per analysis by default for the purpose of example. Change "nsim" in script to run for a different number of simulations
# Will produce many CSV files and PDF plots in subfolder "OutputFiles"
source("run_02_brownie.R")

# Generates Supplemental Table 19 - ARD vs ER rates for binary traits
source("run_03_binary_transition_rates.R")

# Generates plots and components of: Figures 3B-D; Supplemental Tables 7, 8, 10, 11; Extended Data Figures 4A&B, 5; Extended Data Table 1 ----
count_transitions = FALSE # change to TRUE if you want transition counts and arrow plot outputs; note: takes longer, especially if the number of simulations (nsims_real) is high
jackknife_families_above = 65 # Jackknife analysis: Defaults to only testing removal of families with at least 65 species present, for the sake of example. For the publication, we performed jackknife analyses iteratively removing each family with at least 3 species in our dataset; change this value to "3" to perform the full jackknifing analysis
run_all_sociality_traits = FALSE # if FALSE, will run only binary traits HighConfidence_Coop, MeanCoopTie2Noncoop, Griesser2017FamilialLiving, Griesser2023.MoreThanTwoCaretakers versus FemaleSong_Agg01; if TRUE, will run all sociality traits and alternative cooperative breeding classifications
source("run_04_simmap_overlap.R")

# Generates Supplemental Table 13, 14 and calculates number of species to downsample in next script
source("run_05_bias_analyses.R")

### Phylopath ----
# Note: you will likely get a warning when running phylopath that says something like "In check_models_data_tree(model_set, data, tree, na.rm) : Column HighConfidence_Coop appears to have binary data, but was not recognized as binary. If it should be treated as binary, convert it to a factor first." This seems to appear even when features are converted to factors, at least in phylopath package version 1.3.0. Thus, please disregard this warning. 

# Generates Figure 4A, Extended Data Figure 7A:
source("run_06_phylopath_main.R")

# Generates Figure 4B-D, Extended Data Figure 7B-D:
# Performs analyses used and calls script to make figures
# Must run "run_05_bias_analyses.R" and "run_06_phylopath_main.R" before running
n_iterations = 2 
source("run_07_phylopath_downsampling.R")

# Generates Supplemental Table 15; Extended Data Figure 6:
source("run_08_phylopath_altCoop.R")


