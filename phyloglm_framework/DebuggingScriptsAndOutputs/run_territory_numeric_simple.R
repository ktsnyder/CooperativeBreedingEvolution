# Simple script to run phyloglm with Territory_12vs3 as numeric

# Load required functions
source("phyloglm_framework/batch_runner.R")

# Load and prepare data
cat("Loading data...\n")
data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Convert Territory_12vs3 to numeric
# First check current values
cat("\nCurrent Territory_12vs3 values:\n")
print(table(data$Territory_12vs3, useNA = "always"))

# Convert to numeric - assuming values are "1vs2" and "3"
data$Territory_12vs3_numeric <- as.numeric(factor(data$Territory_12vs3))
cat("\nNumeric Territory_12vs3 values:\n")
print(table(data$Territory_12vs3_numeric, useNA = "always"))

# If the values are actually "1vs2" and "3", let's map them to 1.5 and 3
if (all(unique(na.omit(data$Territory_12vs3)) %in% c("1vs2", "3"))) {
  data$Territory_12vs3_numeric <- ifelse(data$Territory_12vs3 == "1vs2", 1.5, 
                                         ifelse(data$Territory_12vs3 == "3", 3, NA))
  cat("\nRemapped Territory_12vs3_numeric values (1vs2->1.5, 3->3):\n")
  print(table(data$Territory_12vs3_numeric, useNA = "always"))
}

# Save the modified data temporarily
saveRDS(data, "temp_data_with_numeric_territory.rds")

# Create configurations using the numeric variable
configs <- list()

# FS_vs_CB_Terr3cont_Mass
configs$FS_vs_CB_Terr3cont_Mass <- list(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "Territory_12vs3_numeric"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  name = "FS_vs_CB_Terr3cont_Mass"
)

# CB_vs_FS_Terr3cont_Mass
configs$CB_vs_FS_Terr3cont_Mass <- list(
  response = "HighConfidence_Coop",
  predictors = c("FemaleSong_Agg01", "Territory_12vs3_numeric"),
  controls = c("logMass_AVONET"),
  transformations = list(),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  name = "CB_vs_FS_Terr3cont_Mass"
)

# Run the batch analysis
cat("\nRunning analyses with numeric territory variable...\n")
results <- run_phyloglm_batch(
  configs = configs,
  data_path = "temp_data_with_numeric_territory.rds",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = "Outputs/PhyloglmResults/numeric_territory",
  n_bootstrap = 1000,
  parallel = FALSE,
  save_intermediate = TRUE
)

cat("\nAnalyses complete! Results saved to Outputs/PhyloglmResults/numeric_territory/\n")

# Clean up temporary file
file.remove("temp_data_with_numeric_territory.rds")