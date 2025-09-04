#### Full Runner Script ----

# Generates Supplemental Table 3 - phylANOVA
source("run_01_phylANOVA.R")

# Generates components used in: Figure 1; Extended Data Figures 1, 2; Supplemental Tables 4, 5, 6
# Performs only 20 simulations per analysis by default for the purpose of example. Change "nsim" in script to run for a different number of simulations
# Will produce many CSV files and PDF plots in subfolder "OutputFiles"
source("run_02_brownie.R")

# Generates Supplemental Table 19 - ARD vs ER rates for binary traits
source("run_03_binary_transition_rates.R")

# Generates components of: Supplemental Tables 7, 8, 10, 11; Extended Data Figure 5A&B; Extended Data Table 1 ----
source("run_04_simmap_overlap.R")
