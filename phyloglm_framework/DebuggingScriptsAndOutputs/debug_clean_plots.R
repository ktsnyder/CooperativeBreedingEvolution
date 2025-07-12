# Debug script for clean plots

# Load the results
batch_results <- readRDS('./Outputs/PhyloglmResults/PhyloGLM_Batch_20250710_235142_boot100/all_results.rds')

# Source the clean plotting functions
source("phyloglm_framework/phyloglm_summary_plots_clean.R")

# Debug the forest plot function
debug_forest_plot <- function(batch_results) {
  predictor_var = "HighConfidence_Coop"
  response_var = "FemaleSong_Agg01"
  predictor_label = "Cooperative Breeding"
  response_label = "Female Song"
  
  # Collect relevant effects
  plot_data <- data.frame()
  
  # Process first analysis for debugging
  for (name in names(batch_results$results)[1:2]) {
    cat("\nProcessing:", name, "\n")
    
    # Skip unwanted analyses
    if (grepl("Terr3", name)) {
      cat("  Skipping - has Terr3\n")
      next
    }
    if (grepl("Migration|Region", name)) {
      cat("  Skipping - has Migration/Region\n") 
      next
    }
    
    result <- batch_results$results[[name]]
    if (!result$success) {
      cat("  Skipping - not successful\n")
      next
    }
    
    # Check if this analysis has our variables
    config <- result$config
    has_forward <- config$response == response_var && predictor_var %in% config$predictors
    has_reverse <- config$response == predictor_var && response_var %in% config$predictors
    
    cat("  Has forward:", has_forward, "\n")
    cat("  Has reverse:", has_reverse, "\n")
    
    if (!has_forward && !has_reverse) {
      cat("  Skipping - no bidirectional relationship\n")
      next
    }
    
    # Must have territoriality
    if (!"TerritorialityWeakVsStrong" %in% config$predictors) {
      cat("  Skipping - no territoriality\n")
      next
    }
    
    # Get best model
    best_model <- result$comparison$comparison[1, ]
    cat("  Best model:", best_model$Model, "\n")
    
    # Get coefficient for the relevant parameter
    param_name <- if (has_forward) predictor_var else response_var
    cat("  Looking for parameter:", param_name, "\n")
    
    coef_data <- result$coefficients %>%
      filter(Model == best_model$Model, 
             Parameter == param_name)
    
    cat("  Found", nrow(coef_data), "rows\n")
    
    if (nrow(coef_data) == 0) {
      cat("  Skipping - no coefficient data\n")
      next
    }
    
    # Print coefficient data structure
    cat("  Coefficient data structure:\n")
    str(coef_data)
    
    # Extract control variable
    control <- if ("logMass_AVONET" %in% config$controls) {
      "Body Mass"
    } else if ("PercentAbsLogWingDimorphism" %in% config$controls) {
      "Wing Dimorphism"  
    } else if ("logMaleFemalePlumageDiffAbs" %in% config$controls) {
      "Plumage Dimorphism"
    } else {
      "Other"
    }
    
    # Try to create entry
    cat("  Creating entry...\n")
    
    # Debug CI values
    cat("  CI_lower in coef_data:", "CI_lower" %in% names(coef_data), "\n")
    if ("CI_lower" %in% names(coef_data)) {
      cat("  CI_lower value:", coef_data$CI_lower[1], "\n")
      cat("  CI_lower class:", class(coef_data$CI_lower), "\n")
    }
    
    # Create CI values safely
    ci_lower <- tryCatch({
      if ("CI_lower" %in% names(coef_data) && length(coef_data$CI_lower) > 0 && !is.na(coef_data$CI_lower[1])) {
        as.numeric(coef_data$CI_lower[1])
      } else {
        as.numeric(coef_data$Estimate[1]) - 1.96 * as.numeric(coef_data$StdErr[1])
      }
    }, error = function(e) {
      cat("  Error calculating CI_lower:", e$message, "\n")
      NA
    })
    
    ci_upper <- tryCatch({
      if ("CI_upper" %in% names(coef_data) && length(coef_data$CI_upper) > 0 && !is.na(coef_data$CI_upper[1])) {
        as.numeric(coef_data$CI_upper[1])
      } else {
        as.numeric(coef_data$Estimate[1]) + 1.96 * as.numeric(coef_data$StdErr[1])
      }
    }, error = function(e) {
      cat("  Error calculating CI_upper:", e$message, "\n")
      NA
    })
    
    cat("  Final CI_lower:", ci_lower, "\n")
    cat("  Final CI_upper:", ci_upper, "\n")
    
    # Create plot entry
    entry <- data.frame(
      Analysis = name,
      Direction = if (has_forward) {
        paste0(predictor_label, " → ", response_label)
      } else {
        paste0(response_label, " → ", predictor_label)
      },
      Control = control,
      Model = best_model$Model,
      Estimate = as.numeric(coef_data$Estimate[1]),
      StdErr = as.numeric(coef_data$StdErr[1]),
      p_value = as.numeric(coef_data$p_value[1]),
      OddsRatio = exp(as.numeric(coef_data$Estimate[1])),
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      stringsAsFactors = FALSE
    )
    
    plot_data <- rbind(plot_data, entry)
  }
  
  return(plot_data)
}

# Run the debug
debug_data <- debug_forest_plot(batch_results)
cat("\n\nFinal plot_data:\n")
print(debug_data)