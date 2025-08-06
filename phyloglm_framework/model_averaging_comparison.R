# Compare top two models and show weighted average effects
library(phylolm)
library(ape)

# Load results
results <- readRDS('Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds')
fs_cb_terr3 <- results[['FS_vs_CB_Terr3_Mass']]

# Get model comparison
comparison <- fs_cb_terr3$comparison$comparison
cat("Model comparison (top 5):\n")
print(comparison[1:5,])

# Get the two top models
model1 <- fs_cb_terr3$models$models[["Terr_Mass"]]
model2 <- fs_cb_terr3$models$models[["MainPred_Terr_Mass"]]

# Extract coefficients
coef1 <- summary(model1)$coefficients
coef2 <- summary(model2)$coefficients

cat("\n\nModel 1 (Terr_Mass) coefficients:\n")
print(coef1)

cat("\n\nModel 2 (MainPred_Terr_Mass) coefficients:\n")
print(coef2)

# Calculate weighted averages for shared coefficients
weights <- comparison$weight[1:2]
cat("\n\nModel weights:", round(weights, 3), "\n")
cat("Weight ratio:", round(weights[1]/weights[2], 2), "\n")

# For Territory and Mass effects (present in both models)
shared_params <- intersect(rownames(coef1), rownames(coef2))
cat("\n\nWeighted average coefficients for shared parameters:\n")
for (param in shared_params) {
  avg_est <- weights[1] * coef1[param, "Estimate"] + weights[2] * coef2[param, "Estimate"]
  avg_se <- sqrt(weights[1]^2 * coef1[param, "StdErr"]^2 + weights[2]^2 * coef2[param, "StdErr"]^2)
  cat(param, ":\n")
  cat("  Weighted estimate:", round(avg_est, 4), "\n")
  cat("  Weighted SE:", round(avg_se, 4), "\n")
  cat("  Model 1:", round(coef1[param, "Estimate"], 4), "±", round(coef1[param, "StdErr"], 4), "\n")
  cat("  Model 2:", round(coef2[param, "Estimate"], 4), "±", round(coef2[param, "StdErr"], 4), "\n\n")
}

# Show CB effect (only in model 2)
cat("HighConfidence_Coop effect (only in Model 2):\n")
cat("  Estimate:", round(coef2["HighConfidence_Coop", "Estimate"], 4), "\n")
cat("  SE:", round(coef2["HighConfidence_Coop", "StdErr"], 4), "\n")
cat("  p-value:", round(coef2["HighConfidence_Coop", "p.value"], 4), "\n")

# Calculate model-averaged prediction for CB effect
# This represents the "average" effect accounting for model uncertainty
cb_avg_effect <- weights[1] * 0 + weights[2] * coef2["HighConfidence_Coop", "Estimate"]
cb_avg_se <- sqrt(weights[2]^2 * coef2["HighConfidence_Coop", "StdErr"]^2 + 
                  weights[1] * weights[2] * coef2["HighConfidence_Coop", "Estimate"]^2)

cat("\nModel-averaged CB effect:\n")
cat("  Averaged estimate:", round(cb_avg_effect, 4), "\n")
cat("  Averaged SE (including model uncertainty):", round(cb_avg_se, 4), "\n")
cat("  Effective weight:", round(weights[2], 3), "\n")