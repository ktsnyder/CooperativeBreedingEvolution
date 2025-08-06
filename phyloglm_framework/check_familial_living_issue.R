# Check the Familial Living sample size issue

# Load the comprehensive results
results <- readRDS("Outputs/PhyloglmResults/stepwise_comprehensive/FS_vs_CB_TerrWS_Mass/expansion_result.rds")

cat("BASE MODEL:\n")
cat("N species:", results$base_n_species, "\n")
cat("Base AIC:", results$base_aic, "\n\n")

cat("PREDICTORS AND SAMPLE SIZES:\n")
print(results$improvements[, c("Predictor", "N_Species", "AIC_Improvement")])

# Calculate percent of species retained
results$improvements$Percent_Species <- 
  round(100 * results$improvements$N_Species / results$base_n_species, 1)

cat("\n\nPERCENT OF SPECIES RETAINED:\n")
print(results$improvements[, c("Predictor", "N_Species", "Percent_Species", "AIC_Improvement")])

# Load the actual data to check Familial Living
data <- read.csv("Data_R_2025-06-09.csv")
cat("\n\nFAMILIAL LIVING DATA AVAILABILITY:\n")
cat("Total species in dataset:", nrow(data), "\n")
cat("Species with Familial Living data:", sum(!is.na(data$Griesser2017FamilialLiving)), "\n")
cat("Percent with data:", round(100 * sum(!is.na(data$Griesser2017FamilialLiving)) / nrow(data), 1), "%\n")

# Show what happens with AIC
cat("\n\nWHY THIS IS A PROBLEM:\n")
cat("When comparing models with different N:\n")
cat("- Base model: 875 species, AIC = 894.11\n")
cat("- FamilialLiving model: 407 species, AIC = 446.76\n")
cat("- 'Improvement' = 894.11 - 446.76 = 447.35\n")
cat("- But these AICs are NOT comparable!\n")
cat("\nThe correct approach:\n")
cat("1. Fit base model to same 407 species that have FamilialLiving data\n")
cat("2. Then add FamilialLiving and compare AICs\n")
cat("3. This gives the true improvement from adding the predictor\n")