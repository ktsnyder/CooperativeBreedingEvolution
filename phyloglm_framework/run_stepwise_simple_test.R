# Simplified stepwise expansion test
# This version uses the exact data from the original analysis

library(phylolm)
library(dplyr)

# Load results
cat("Loading results...\n")
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")

# Get base analysis
base_result <- all_results[["FS_vs_CB_TerrWS_Mass"]]
best_model_name <- base_result$comparison$comparison$Model[1]
best_model <- base_result$models$models[[best_model_name]]

cat("\nBase model:", best_model_name, "\n")
cat("Base AIC:", base_result$comparison$comparison$AIC[1], "\n")
cat("N species:", nrow(best_model$data), "\n")

# Use the exact data and tree from the base model
base_data <- best_model$data
base_tree <- best_model$phy

# Load full dataset to get additional predictors
full_data <- read.csv("Data_R_2025-06-09.csv")

# Function to test adding a predictor
test_add_predictor <- function(base_model, new_var_name, new_var_data) {
  # Get base formula
  base_formula <- formula(base_model)
  
  # Add new variable to data
  model_data <- base_model$data
  
  # Match species to ensure alignment
  if ("species" %in% names(new_var_data)) {
    new_var_data <- new_var_data[match(rownames(model_data), new_var_data$species), ]
  }
  
  model_data[[new_var_name]] <- new_var_data[[new_var_name]]
  
  # Remove rows with NA in new variable
  complete_rows <- complete.cases(model_data)
  if (sum(complete_rows) < nrow(model_data)) {
    cat("  Reduced from", nrow(model_data), "to", sum(complete_rows), "species due to missing data\n")
  }
  
  model_data <- model_data[complete_rows, ]
  model_tree <- keep.tip(base_model$phy, rownames(model_data))
  
  # Ensure proper ordering
  model_data <- model_data[match(model_tree$tip.label, rownames(model_data)), ]
  
  # Create new formula
  new_formula <- update(base_formula, paste("~ . +", new_var_name))
  
  # Try fitting
  tryCatch({
    new_model <- phyloglm(
      formula = new_formula,
      data = model_data,
      phy = model_tree,
      method = "logistic_MPLE",
      btol = 30,
      log.alpha.bound = 4
    )
    
    base_aic <- -2 * base_model$logLik + 2 * base_model$d
    new_aic <- -2 * new_model$logLik + 2 * new_model$d
    improvement <- base_aic - new_aic
    
    cat("  Success! AIC improvement:", round(improvement, 2), "\n")
    
    return(list(
      success = TRUE,
      model = new_model,
      aic_improvement = improvement,
      n_species = nrow(model_data)
    ))
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}

# Test adding predictors
cat("\n=== Testing Additional Predictors ===\n")

# 1. Territory as numeric
cat("\n1. Testing Territory (1-3 numeric):\n")
territory_data <- full_data[, c("species", "Territory")]
territory_data$Territory <- as.numeric(territory_data$Territory)
result_territory <- test_add_predictor(best_model, "Territory", territory_data)

# 2. Migration as numeric  
cat("\n2. Testing Migration (1-3 numeric):\n")
migration_data <- full_data[, c("species", "Migration_AVONET")]
migration_data$Migration_AVONET <- as.numeric(migration_data$Migration_AVONET)
result_migration <- test_add_predictor(best_model, "Migration_AVONET", migration_data)

# 3. Geographic Region
cat("\n3. Testing Geographic Region:\n")
region_data <- full_data[, c("species", "GeographicRegion_Jetz")]
result_region <- test_add_predictor(best_model, "GeographicRegion_Jetz", region_data)

# 4. Absolute Latitude
cat("\n4. Testing Absolute Latitude:\n")
lat_data <- full_data[, c("species", "Centroid.Latitude_AVONET")]
lat_data$abs_latitude <- abs(lat_data$Centroid.Latitude_AVONET)
result_lat <- test_add_predictor(best_model, "abs_latitude", lat_data[, c("species", "abs_latitude")])

# 5. Plumage Dimorphism (if not already in model)
if (!("logMaleFemalePlumageDiffAbs" %in% names(base_data))) {
  cat("\n5. Testing Plumage Dimorphism:\n")
  plum_data <- full_data[, c("species", "logMaleFemalePlumageDiffAbs")]
  result_plum <- test_add_predictor(best_model, "logMaleFemalePlumageDiffAbs", plum_data)
}

# 6. Wing Dimorphism (if not already in model)
if (!("PercentAbsLogWingDimorphism" %in% names(base_data))) {
  cat("\n6. Testing Wing Dimorphism:\n")
  wing_data <- full_data[, c("species", "PercentAbsLogWingDimorphism")]
  result_wing <- test_add_predictor(best_model, "PercentAbsLogWingDimorphism", wing_data)
}

# Summary
cat("\n=== Summary ===\n")
results <- list(
  Territory = result_territory,
  Migration = result_migration,
  GeographicRegion = result_region,
  AbsLatitude = result_lat
)

kept <- sapply(results, function(x) x$success && x$aic_improvement > 2)
cat("\nPredictors to keep (AIC improvement > 2):\n")
for (name in names(kept)[kept]) {
  cat("  -", name, "(AIC improved by", round(results[[name]]$aic_improvement, 2), ")\n")
}

# If Territory was successful, build final model with it
if (result_territory$success && result_territory$aic_improvement > 2) {
  cat("\n=== Building Final Model with Territory ===\n")
  
  # Prepare data
  model_data <- best_model$data
  territory_data_matched <- territory_data[match(rownames(model_data), territory_data$species), ]
  model_data$Territory <- territory_data_matched$Territory
  
  # Final formula
  final_formula <- update(formula(best_model), ~ . + Territory)
  cat("Final formula:", deparse(final_formula), "\n")
  
  # Fit final model
  final_model <- phyloglm(
    formula = final_formula,
    data = model_data,
    phy = base_tree,
    method = "logistic_MPLE",
    btol = 30,
    log.alpha.bound = 4
  )
  
  cat("\nFinal model coefficients:\n")
  print(summary(final_model)$coefficients)
}