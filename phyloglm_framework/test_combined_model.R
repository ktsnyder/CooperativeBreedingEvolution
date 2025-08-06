# Test combined model with Territory and Latitude

library(phylolm)
library(ape)

# Load data
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
full_data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Get base model for FS_vs_CB_TerrWS_Mass
base_result <- all_results[["FS_vs_CB_TerrWS_Mass"]]
best_model_name <- base_result$comparison$comparison$Model[1]
best_model <- base_result$models$models[[best_model_name]]
base_aic <- base_result$comparison$comparison$AIC[1]
base_data <- base_result$prepared_data$data
base_tree <- base_result$prepared_data$tree

cat("Base model:", best_model_name, "\n")
cat("Base AIC:", base_aic, "\n")
cat("Base formula:", deparse(formula(best_model)), "\n")
cat("N species:", nrow(base_data), "\n\n")

# First, fit Territory model
cat("Step 1: Adding Territory\n")
territory_numeric <- as.numeric(full_data$Territory[match(base_data$species, full_data$species)])
base_data$Territory_num <- territory_numeric
test_data <- base_data[!is.na(base_data$Territory_num), ]
test_tree <- keep.tip(base_tree, test_data$species)
test_data <- test_data[match(test_tree$tip.label, test_data$species), ]

territory_formula <- update(formula(best_model), ~ . + Territory_num)
territory_model <- phyloglm(
  formula = territory_formula,
  data = test_data,
  phy = test_tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4
)

territory_aic <- -2 * territory_model$logLik + 2 * territory_model$d
cat("Territory model AIC:", territory_aic, "\n")
cat("Improvement from base:", base_aic - territory_aic, "\n\n")

# Now add Latitude to the Territory model
cat("Step 2: Adding Latitude to Territory model\n")
abs_lat <- abs(full_data$Centroid.Latitude_AVONET[match(test_data$species, full_data$species)])
test_data$abs_Latitude <- abs_lat

# Check NAs
cat("Species with Latitude data:", sum(!is.na(test_data$abs_Latitude)), "\n")

# Filter for complete cases
combined_data <- test_data[!is.na(test_data$abs_Latitude), ]
cat("N species after filtering:", nrow(combined_data), "\n")

combined_tree <- keep.tip(test_tree, combined_data$species)
combined_data <- combined_data[match(combined_tree$tip.label, combined_data$species), ]

# Fit combined model
combined_formula <- update(territory_formula, ~ . + abs_Latitude)
cat("Combined formula:", deparse(combined_formula), "\n")

combined_model <- phyloglm(
  formula = combined_formula,
  data = combined_data,
  phy = combined_tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4
)

combined_aic <- -2 * combined_model$logLik + 2 * combined_model$d
cat("\nCombined model AIC:", combined_aic, "\n")
cat("Improvement from base:", base_aic - combined_aic, "\n")
cat("Improvement from Territory model:", territory_aic - combined_aic, "\n")

# Show coefficients
cat("\nCombined model coefficients:\n")
print(summary(combined_model)$coefficients)