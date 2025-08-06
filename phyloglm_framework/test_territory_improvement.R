# Quick test to verify Territory improvement

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

# Add Territory as numeric
territory_numeric <- as.numeric(full_data$Territory[match(base_data$species, full_data$species)])
base_data$Territory_num <- territory_numeric

# Check how many species have Territory data
cat("Species with Territory data:", sum(!is.na(base_data$Territory_num)), "\n")

# Keep only species with Territory data
test_data <- base_data[!is.na(base_data$Territory_num), ]
test_tree <- keep.tip(base_tree, test_data$species)
test_data <- test_data[match(test_tree$tip.label, test_data$species), ]

cat("N species for test:", nrow(test_data), "\n")

# Update formula
new_formula <- update(formula(best_model), ~ . + Territory_num)
cat("New formula:", deparse(new_formula), "\n")

# Fit model
new_model <- phyloglm(
  formula = new_formula,
  data = test_data,
  phy = test_tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4
)

new_aic <- -2 * new_model$logLik + 2 * new_model$d
improvement <- base_aic - new_aic

cat("\nNew AIC:", new_aic, "\n")
cat("AIC improvement:", improvement, "\n")
cat("Territory coefficient:", coef(new_model)["Territory_num"], "\n")
cat("Territory p-value:", summary(new_model)$coefficients["Territory_num", "p.value"], "\n")

# Show all coefficients
cat("\nAll coefficients:\n")
print(summary(new_model)$coefficients)