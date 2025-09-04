#### Full Runner Script ----

# Generates Supplemental Table 3 - phylANOVA
source("run_01_phylANOVA.R")

# Generates components used in: Figure 1; Extended Data Figures 1, 2; Supplemental Tables 4, 5, 6
# Performs only 20 simulations per analysis by default for the purpose of example. Change "nsim" in script to run for a different number of simulations
# Will produce many CSV files and PDF plots in subfolder "OutputFiles"
source("run_02_brownie.R")