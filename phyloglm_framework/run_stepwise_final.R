# Final stepwise expansion script with visualization
# This version saves results and creates plots

source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/stepwise_visualization.R")

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

# Main stepwise expansion function
run_stepwise_expansion_simple <- function(output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data and results
  cat("Loading data and results...\n")
  all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
  full_data <- read.csv("Data_R_2025-06-09.csv")
  tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  
  # Function to expand a model (from run_stepwise_working.R)
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
        btol = 50,  # Increased for stability
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
  
  # Store all results
  expansion_results <- list()
  
  # Test expansions for both directions
  analyses <- c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass")
  
  for (analysis_name in analyses) {
    cat("\n========================================\n")
    cat("Analyzing:", analysis_name, "\n")
    cat("========================================\n")
    
    base_result <- all_results[[analysis_name]]
    base_aic <- base_result$comparison$comparison$AIC[1]
    
    cat("\nBase model AIC:", base_aic, "\n")
    
    # Initialize result storage
    analysis_results <- list(
      base_result = base_result,
      base_aic = base_aic,
      tests = list(),
      improvements = data.frame()
    )
    
    # Test predictors
    cat("\nTesting additional predictors:\n")
    
    # 1. Territory as numeric
    cat("\n1. Territory (1-3 numeric):\n")
    territory_numeric <- as.numeric(full_data$Territory)
    result_terr <- expand_model(base_result, "Territory_num", territory_numeric)
    analysis_results$tests$Territory <- result_terr
    
    # 2. Migration as numeric
    cat("\n2. Migration (1-3 numeric):\n")
    migration_numeric <- as.numeric(full_data$Migration_AVONET)
    result_mig <- expand_model(base_result, "Migration_num", migration_numeric)
    analysis_results$tests$Migration <- result_mig
    
    # 3. Geographic Region
    cat("\n3. Geographic Region:\n")
    result_region <- expand_model(base_result, "GeographicRegion_Jetz")
    analysis_results$tests$Region <- result_region
    
    # 4. Absolute Latitude (if not already in model)
    if (!grepl("absLat|Latitude", analysis_name)) {
      cat("\n4. Absolute Latitude:\n")
      abs_lat <- abs(full_data$Centroid.Latitude_AVONET)
      result_lat <- expand_model(base_result, "abs_Latitude", abs_lat)
      analysis_results$tests$Latitude <- result_lat
    }
    
    # 5. Plumage Dimorphism (if not already in model)
    if (!grepl("PlumDim", analysis_name)) {
      cat("\n5. Plumage Dimorphism:\n")
      result_plum <- expand_model(base_result, "logMaleFemalePlumageDiffAbs")
      analysis_results$tests$Plumage <- result_plum
    }
    
    # 6. Wing Dimorphism (if not already in model)
    if (!grepl("WingDim", analysis_name)) {
      cat("\n6. Wing Dimorphism:\n")
      result_wing <- expand_model(base_result, "PercentAbsLogWingDimorphism")
      analysis_results$tests$Wing <- result_wing
    }
    
    # Create improvement summary
    for (test_name in names(analysis_results$tests)) {
      test <- analysis_results$tests[[test_name]]
      if (test$success) {
        cat("   ", test_name, "- Success! AIC improvement:", round(test$aic_improvement, 2), "\n")
        analysis_results$improvements <- rbind(
          analysis_results$improvements,
          data.frame(
            Predictor = test_name,
            AIC_Improvement = test$aic_improvement,
            N_Species = test$n_species,
            Kept = test$aic_improvement > 2,
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("   ", test_name, "- Failed:", test$error, "\n")
        analysis_results$improvements <- rbind(
          analysis_results$improvements,
          data.frame(
            Predictor = test_name,
            AIC_Improvement = NA,
            N_Species = NA,
            Kept = FALSE,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    # Build final model with improvements > 2
    kept_predictors <- analysis_results$improvements$Predictor[
      !is.na(analysis_results$improvements$AIC_Improvement) & 
      analysis_results$improvements$AIC_Improvement > 2
    ]
    
    if (length(kept_predictors) > 0) {
      cat("\n\nBuilding expanded model with:", paste(kept_predictors, collapse = ", "), "\n")
      
      # Start with base model and iteratively add kept predictors
      current_result <- base_result
      
      for (pred in kept_predictors) {
        test_result <- analysis_results$tests[[pred]]
        if (test_result$success) {
          current_result <- list(
            models = list(models = list(best = test_result$model)),
            comparison = list(comparison = data.frame(Model = "best", AIC = test_result$aic)),
            prepared_data = list(data = test_result$data, tree = test_result$tree)
          )
        }
      }
      
      analysis_results$final_model <- current_result$models$models[[1]]
      analysis_results$final_aic <- -2 * analysis_results$final_model$logLik + 
                                   2 * analysis_results$final_model$d
      analysis_results$total_improvement <- base_aic - analysis_results$final_aic
      
      cat("Final AIC:", round(analysis_results$final_aic, 2), "\n")
      cat("Total improvement:", round(analysis_results$total_improvement, 2), "\n")
      
      # Calculate effects
      # First extract coefficients from the model
      coef_summary <- summary(analysis_results$final_model)$coefficients
      final_coefs <- data.frame(
        Parameter = rownames(coef_summary),
        Estimate = coef_summary[, "Estimate"],
        StdErr = coef_summary[, "StdErr"],
        z_value = coef_summary[, "z.value"],
        P.Value = coef_summary[, "p.value"],
        stringsAsFactors = FALSE
      )
      
      analysis_results$final_effects <- calculate_effect_sizes(final_coefs)
    }
    
    expansion_results[[analysis_name]] <- analysis_results
  }
  
  # Save results
  saveRDS(expansion_results, file.path(output_dir, "stepwise_results.rds"))
  
  # Create visualizations
  cat("\n\nCreating visualizations...\n")
  
  # 1. Improvement barplot
  all_improvements <- data.frame()
  for (name in names(expansion_results)) {
    imp <- expansion_results[[name]]$improvements
    imp$Analysis <- name
    all_improvements <- rbind(all_improvements, imp)
  }
  
  p_improvements <- ggplot(all_improvements, 
                          aes(x = reorder(Predictor, AIC_Improvement), 
                              y = AIC_Improvement, 
                              fill = Analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "Stepwise Model Improvements",
         subtitle = "AIC improvement when adding each predictor (threshold = 2)",
         x = "", y = "AIC Improvement") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")
  
  ggsave(file.path(output_dir, "stepwise_improvements.png"), 
         p_improvements, width = 10, height = 6, dpi = 300)
  
  # 2. Final model summaries
  summary_text <- "\nSTEPWISE EXPANSION SUMMARY\n"
  summary_text <- paste0(summary_text, "==========================\n\n")
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    summary_text <- paste0(summary_text, name, "\n")
    summary_text <- paste0(summary_text, paste(rep("-", nchar(name)), collapse = ""), "\n")
    summary_text <- paste0(summary_text, "Base AIC: ", round(result$base_aic, 2), "\n")
    
    if ("final_model" %in% names(result)) {
      summary_text <- paste0(summary_text, "Final AIC: ", round(result$final_aic, 2), "\n")
      summary_text <- paste0(summary_text, "Total improvement: ", round(result$total_improvement, 2), "\n")
      summary_text <- paste0(summary_text, "Predictors added: ", 
                            paste(result$improvements$Predictor[result$improvements$Kept], 
                                  collapse = ", "), "\n")
    } else {
      summary_text <- paste0(summary_text, "No predictors improved model by >2 AIC units\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "stepwise_summary.txt"))
  cat(summary_text)
  
  # 3. Forest plots if we have expanded models
  forest_plots <- list()
  for (name in names(expansion_results)) {
    if ("final_effects" %in% names(expansion_results[[name]])) {
      effects <- expansion_results[[name]]$final_effects
      effects <- effects[effects$Variable != "(Intercept)", ]
      
      # Clean variable names
      effects$Variable <- effects$Parameter  # Use Parameter column
      effects$Variable <- gsub("HighConfidence_Coop", "Cooperative Breeding", effects$Variable)
      effects$Variable <- gsub("FemaleSong_Agg01", "Female Song", effects$Variable)
      effects$Variable <- gsub("TerritorialityWeakVsStrong", "Territoriality", effects$Variable)
      effects$Variable <- gsub("Territory_num", "Territory (1-3)", effects$Variable)
      effects$Variable <- gsub("Migration_num", "Migration (1-3)", effects$Variable)
      effects$Variable <- gsub(":", " × ", effects$Variable)
      
      p <- ggplot(effects, aes(x = OddsRatio, y = reorder(Variable, OddsRatio))) +
        geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
        geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), height = 0.2) +
        geom_point(aes(color = P.Value < 0.05), size = 3) +
        scale_x_log10() +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50"),
                          labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
        labs(title = paste("Expanded Model:", name),
             x = "Odds Ratio (95% CI)", y = "",
             color = "Significant") +
        theme_minimal()
      
      forest_plots[[name]] <- p
    }
  }
  
  if (length(forest_plots) > 0) {
    combined_forest <- wrap_plots(forest_plots, ncol = 1)
    ggsave(file.path(output_dir, "expanded_model_effects.png"),
           combined_forest, width = 10, height = 6 * length(forest_plots), dpi = 300)
  }
  
  cat("\n\nResults saved to:", output_dir, "\n")
  
  return(expansion_results)
}

# Run the analysis
results <- run_stepwise_expansion_simple()