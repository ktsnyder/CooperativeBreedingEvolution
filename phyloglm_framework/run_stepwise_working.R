# Working stepwise expansion for phyloglm models
# This version properly handles the phyloglm object structure

source("phyloglm_framework/batch_runner_helpers.R")  # For calculate_effect_sizes

library(phylolm)
library(ape)
library(dplyr)

# Load data and results
cat("Loading data and results...\n")
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
full_data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Function to expand a model
expand_model <- function(base_result, predictor_name, predictor_values = NULL) {
  # Get best model
  best_model_name <- base_result$comparison$comparison$Model[1]
  best_model <- base_result$models$models[[best_model_name]]
  base_aic <- base_result$comparison$comparison$AIC[1]
  
  # Get species from the model
  model_species <- rownames(best_model$X)
  
  # Get original data used
  orig_data <- base_result$prepared_data$data
  orig_tree <- base_result$prepared_data$tree
  
  # If predictor values not provided, get from full_data
  if (is.null(predictor_values)) {
    if (!(predictor_name %in% names(full_data))) {
      return(list(success = FALSE, error = "Predictor not found in data"))
    }
    pred_data <- full_data[, c("species", predictor_name)]
    names(pred_data)[2] <- predictor_name
  } else {
    pred_data <- data.frame(
      species = full_data$species,
      pred = predictor_values
    )
    names(pred_data)[2] <- predictor_name
  }
  
  # Merge with original data
  expanded_data <- merge(orig_data, pred_data, by = "species", all.x = TRUE)
  
  # Remove rows with NA in new predictor
  expanded_data <- expanded_data[!is.na(expanded_data[[predictor_name]]), ]
  
  if (nrow(expanded_data) < nrow(orig_data)) {
    cat("  Data reduced from", nrow(orig_data), "to", nrow(expanded_data), "species\n")
  }
  
  # Match tree
  keep_species <- intersect(expanded_data$species, orig_tree$tip.label)
  
  if (length(keep_species) == 0) {
    return(list(success = FALSE, error = "No species overlap between data and tree"))
  }
  
  expanded_data <- expanded_data[expanded_data$species %in% keep_species, ]
  expanded_tree <- keep.tip(orig_tree, keep_species)
  
  # Order data to match tree
  expanded_data <- expanded_data[match(expanded_tree$tip.label, expanded_data$species), ]
  rownames(expanded_data) <- expanded_data$species
  
  # Create new formula
  base_formula <- formula(best_model)
  new_formula <- update(base_formula, paste("~ . +", predictor_name))
  
  # Fit new model
  tryCatch({
    new_model <- phyloglm(
      formula = new_formula,
      data = expanded_data,
      phy = expanded_tree,
      method = "logistic_MPLE",
      btol = 30,
      log.alpha.bound = 4
    )
    
    new_aic <- -2 * new_model$logLik + 2 * new_model$d
    improvement <- base_aic - new_aic
    
    return(list(
      success = TRUE,
      model = new_model,
      formula = new_formula,
      aic = new_aic,
      aic_improvement = improvement,
      n_species = nrow(expanded_data),
      data = expanded_data,
      tree = expanded_tree
    ))
    
  }, error = function(e) {
    return(list(success = FALSE, error = e$message))
  })
}

# Test expansions for both directions
analyses <- c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass")

for (analysis_name in analyses) {
  cat("\n========================================\n")
  cat("Analyzing:", analysis_name, "\n")
  cat("========================================\n")
  
  base_result <- all_results[[analysis_name]]
  base_aic <- base_result$comparison$comparison$AIC[1]
  
  cat("\nBase model AIC:", base_aic, "\n")
  cat("Base model formula:", deparse(formula(base_result$models$models[[base_result$comparison$comparison$Model[1]]])), "\n")
  
  # Test predictors
  cat("\nTesting additional predictors:\n")
  
  # 1. Territory as numeric
  cat("\n1. Territory (1-3 numeric):\n")
  territory_numeric <- as.numeric(full_data$Territory)
  result_terr <- expand_model(base_result, "Territory_num", territory_numeric)
  if (result_terr$success) {
    cat("   Success! AIC improvement:", round(result_terr$aic_improvement, 2), "\n")
    print(result_terr$aic_improvement)  # Force output
  } else {
    cat("   Failed:", result_terr$error, "\n")
    print(result_terr$error)  # Force output
  }
  
  # 2. Migration as numeric
  cat("\n2. Migration (1-3 numeric):\n")
  migration_numeric <- as.numeric(full_data$Migration_AVONET)
  result_mig <- expand_model(base_result, "Migration_num", migration_numeric)
  if (result_mig$success) {
    cat("   Success! AIC improvement:", round(result_mig$aic_improvement, 2), "\n")
  } else {
    cat("   Failed:", result_mig$error, "\n")
  }
  
  # 3. Geographic Region
  cat("\n3. Geographic Region:\n")
  result_region <- expand_model(base_result, "GeographicRegion_Jetz")
  if (result_region$success) {
    cat("   Success! AIC improvement:", round(result_region$aic_improvement, 2), "\n")
  } else {
    cat("   Failed:", result_region$error, "\n")
  }
  
  # 4. Absolute Latitude (if not already in model)
  if (!grepl("absLat|Latitude", analysis_name)) {
    cat("\n4. Absolute Latitude:\n")
    abs_lat <- abs(full_data$Centroid.Latitude_AVONET)
    result_lat <- expand_model(base_result, "abs_Latitude", abs_lat)
    if (result_lat$success) {
      cat("   Success! AIC improvement:", round(result_lat$aic_improvement, 2), "\n")
    } else {
      cat("   Failed:", result_lat$error, "\n")
    }
  }
  
  # 5. Plumage Dimorphism (if not already in model)
  if (!grepl("PlumDim", analysis_name)) {
    cat("\n5. Plumage Dimorphism:\n")
    result_plum <- expand_model(base_result, "logMaleFemalePlumageDiffAbs")
    if (result_plum$success) {
      cat("   Success! AIC improvement:", round(result_plum$aic_improvement, 2), "\n")
    } else {
      cat("   Failed:", result_plum$error, "\n")
    }
  }
  
  # 6. Wing Dimorphism (if not already in model)
  if (!grepl("WingDim", analysis_name)) {
    cat("\n6. Wing Dimorphism:\n")
    result_wing <- expand_model(base_result, "PercentAbsLogWingDimorphism")
    if (result_wing$success) {
      cat("   Success! AIC improvement:", round(result_wing$aic_improvement, 2), "\n")
    } else {
      cat("   Failed:", result_wing$error, "\n")
    }
  }
  
  # Build final model with significant improvements
  results <- list(
    Territory = result_terr,
    Migration = result_mig,
    Region = result_region
  )
  if (exists("result_lat")) results$Latitude <- result_lat
  if (exists("result_plum")) results$Plumage <- result_plum
  if (exists("result_wing")) results$Wing <- result_wing
  
  # Find predictors with AIC improvement > 2
  improvements <- sapply(results, function(x) if(x$success) x$aic_improvement else 0)
  keep <- names(improvements)[improvements > 2]
  
  if (length(keep) > 0) {
    cat("\n\nBuilding expanded model with:", paste(keep, collapse = ", "), "\n")
    
    # Start with base data
    final_data <- base_result$prepared_data$data
    
    # Add each kept predictor
    for (pred in keep) {
      if (pred == "Territory") {
        final_data$Territory_num <- as.numeric(full_data$Territory[match(final_data$species, full_data$species)])
      } else if (pred == "Migration") {
        final_data$Migration_num <- as.numeric(full_data$Migration_AVONET[match(final_data$species, full_data$species)])
      } else if (pred == "Region") {
        final_data$GeographicRegion_Jetz <- full_data$GeographicRegion_Jetz[match(final_data$species, full_data$species)]
      } else if (pred == "Latitude") {
        final_data$abs_Latitude <- abs(full_data$Centroid.Latitude_AVONET[match(final_data$species, full_data$species)])
      }
    }
    
    # Remove NAs
    final_data <- final_data[complete.cases(final_data), ]
    
    # Check we still have species
    if (nrow(final_data) == 0) {
      cat("ERROR: No complete cases after adding predictors\n")
      next
    }
    
    # Match with tree
    keep_species <- intersect(final_data$species, base_result$prepared_data$tree$tip.label)
    if (length(keep_species) == 0) {
      cat("ERROR: No species overlap with tree\n")
      next
    }
    
    final_tree <- keep.tip(base_result$prepared_data$tree, keep_species)
    final_data <- final_data[final_data$species %in% keep_species, ]
    final_data <- final_data[match(final_tree$tip.label, final_data$species), ]
    rownames(final_data) <- final_data$species
    
    # Build formula
    base_formula <- formula(base_result$models$models[[base_result$comparison$comparison$Model[1]]])
    additions <- character()
    if ("Territory" %in% keep) additions <- c(additions, "Territory_num")
    if ("Migration" %in% keep) additions <- c(additions, "Migration_num") 
    if ("Region" %in% keep) additions <- c(additions, "GeographicRegion_Jetz")
    if ("Latitude" %in% keep) additions <- c(additions, "abs_Latitude")
    
    final_formula <- update(base_formula, paste("~ . +", paste(additions, collapse = " + ")))
    
    cat("Final formula:", deparse(final_formula), "\n")
    
    # Fit final model
    final_model <- phyloglm(
      formula = final_formula,
      data = final_data,
      phy = final_tree,
      method = "logistic_MPLE",
      btol = 30,
      log.alpha.bound = 4
    )
    
    final_aic <- -2 * final_model$logLik + 2 * final_model$d
    total_improvement <- base_aic - final_aic
    
    cat("\nFinal model AIC:", round(final_aic, 2), "\n")
    cat("Total AIC improvement:", round(total_improvement, 2), "\n")
    cat("N species:", nrow(final_data), "\n")
    
    # Calculate and display effects
    effects <- calculate_effect_sizes(final_model, "Expanded_Model")
    cat("\nSignificant effects (p < 0.05):\n")
    sig_effects <- effects[effects$P.Value < 0.05 & effects$Variable != "(Intercept)", ]
    for (i in 1:nrow(sig_effects)) {
      cat("  ", sig_effects$Variable[i], ": OR =", round(sig_effects$OR[i], 2), 
          "(p =", round(sig_effects$P.Value[i], 4), ")\n")
    }
  }
  
  # Clean up for next iteration
  to_remove <- intersect(c("result_lat", "result_plum", "result_wing"), ls())
  if (length(to_remove) > 0) rm(list = to_remove)
}