## Full Runner Script ----
### First revision, 9/1/2025

## To run:
# Set working directory to the main folder in this codebase - the folder should contain subdirectory "Source Data Process_CB".
newdata = "Data_R.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

## Check for and load required packages
source("check_required_packages.R")

### Compile data ----
# (OPTIONAL - all analyses use a data csv already compiled using this method - newdata="Data_R.csv")
# See PDF "Source Data Process_CB/Unaltered from publication/Source Data Processes.pdf" for detailed notes on how files containing cooperative breeding classification data were processed and standardized for use in merge_data_allcolumns.R
source(file.path("Source Data Process_CB","merge_data_allcolumns.R"))

### Re-generate consensus trees ----
## (OPTIONAL - all analyses use a consensus tree previously calculated using this method - treefile="ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
## Commented out due to it being one of the more computationally intensive and lengthy processes. 

## Load multitree. Note: This file must be in your working directory. It is too large to be indexed on GitHub, so if you obtained the repo from GitHub, you will need to download trees directly from BirdTree.org.
#passermultitree = read.nexus("PasserineMultiphy1000Hackett4_nondicho.nex")
#nTrees = 10 # set to 1000 to calculate consensus trees as in the manuscript. NOTE: 10 trees, across two consensus methods, may take 15min-1hour. 1000 trees may take multiple days, depending on your system computing capabilities.
#source("run_00_make various consensus trees.R")

### PhylANOVA ----
# Generates Supplemental Table 3 - phylANOVA
source("run_01_phylANOVA.R")

### Brownie ----
# Generates components used in: Figure 1; Extended Data Figures 1, 2; Supplemental Tables 4, 5, 6
# Performs only 20 simulations per analysis by default for the purpose of example. Change "nsim" in script to run for a different number of simulations
# Will produce many CSV files and PDF plots in subfolder "Outputs/Brownie_outputs"
source("run_02_brownie.R")

### Estimating binary state transition rates ----
# Generates Supplemental Table 19 - ARD vs ER rates for binary traits
source("run_03_binary_transition_rates.R")

### Simmap overlap ----
## Generates plots and components of: Figures 3B-D; Supplemental Tables 7, 8, 10, 11; Extended Data Figures 4A&B, 5; Extended Data Table 1 ---
count_transitions = FALSE # change to TRUE if you want transition counts and arrow plot outputs; note: takes longer, especially if the number of simulations (nsims_real) is high
jackknife_families_above = 65 # Jackknife analysis: Defaults to only testing removal of families with at least 65 species present, for the sake of example. For the publication, we performed jackknife analyses iteratively removing each family with at least 3 species in our dataset; change this value to "3" to perform the full jackknifing analysis
run_all_sociality_traits = FALSE # if FALSE, will run only binary traits HighConfidence_Coop, MeanCoopTie2Noncoop, Griesser2017FamilialLiving, Griesser2023.MoreThanTwoCaretakers versus FemaleSong_Agg01; if TRUE, will run all sociality traits and alternative cooperative breeding classifications
source("run_04_simmap_overlap.R")

### Analysis of data biases ----
# Generates Supplemental Table 13, 14 and calculates number of species to downsample in next script
source("run_05_bias_analyses.R")

### Phylopath ----
# Note: you will likely get a warning when running phylopath that says something like "In check_models_data_tree(model_set, data, tree, na.rm) : Column HighConfidence_Coop appears to have binary data, but was not recognized as binary. If it should be treated as binary, convert it to a factor first." This seems to appear even when features are converted to factors, at least in phylopath package version 1.3.0. Thus, please disregard this warning. 

# Generates Figure 4A, Extended Data Figure 7A:
source("run_06_phylopath_main.R")

# Generates Figure 4B-D, Extended Data Figure 7B-D:
# Performs analyses used and calls script to make figures
# Must run "run_05_bias_analyses.R" and "run_06_phylopath_main.R" before running
n_iterations = 2 # Change to 500 to perform analyses as in publication
territoriality_var = "TerritorialityWeakVsStrong"
source("run_07_phylopath_downsampling.R")
territoriality_var = "Territory_12vs3"
source("run_07_phylopath_downsampling.R")

# Generates Supplemental Table 15; Extended Data Figure 6:
source("run_08_phylopath_altCoop.R")


### PhyloGLM ----
## Generates Tables 1-2, Supplemental Tables 15-18, 21-24
# 
# Base phyloglm (Table 1, Supplemental Table 17) and with alt Coops (Supplemental Tables 15-16)
# # Results output to Outputs/PhyloGLM_outputs/PhyloGLM_Batch_[Date]_[Time]_boot100/individual_analyses/[analysis name]/best_model_effects.csv
# 
# Direct-comparison Coop vs Sociality factors (Table 2, Supplemental Table 18)
# # Results output to Outputs/PhyloGLM_outputs/Replace_Coop_With_Soc/Replace_Coop_Additive/Replace_Coop_Summary_Table_boot[nBoot].csv 
# 
# Stepwise iterative (Supplemental Tables 21-24)
# # Results output to Outputs/PhyloGLM_outputs/stepwise_[date]_[time]
nBoot = 100 # number of bootstrap iterations for test run. Set to 500 to replicate analyses performed for the manuscript
run_alt_coop = FALSE # set to TRUE to run base phyloglm analyses for alternative cooperative breeding classifications
run_stepwise = TRUE
source("run_09_phyloglm.R")

### Multitree tests ----
## Performs Brownie and Simmap Overlap tests across many randomly sampled individual trees taken from BirdTree.org
## Generates Extended Data Figures 1A, 3A
nTreesToSample = 10 # set to 200 to run the full analyses performed for the manuscript
nSimsPerTree = 5 # set to 20 to run the full analyses performed for the manuscript
# Load multitree. Note: This file must be in your working directory. It is too large to be indexed on GitHub, so if you obtained the repo from GitHub, you will need to download trees directly from BirdTree.org.
multitree = read.tree("BirdzillaHackett4_Stage2_1000trees.tre")
source("run_10_multitree_tests.R")

### Species counts tables ----
## Generates Supplemental Tables 1, 2, 9, 12, and 25
source("run_11_generate_species_counts_tables.R")
