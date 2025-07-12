# Script to reload saved PhyloGLM results and create/modify plots
# This allows you to work with results without re-running analyses

library(ggplot2)
library(dplyr)
library(patchwork)

# Source visualization framework
source("phyloglm_framework/visualization_framework.R")

#' Load results from a single analysis directory
#' 
#' @param analysis_dir Path to analysis directory
#' @return List with all saved results
load_single_analysis <- function(analysis_dir) {
  
  results <- list()
  
  # Load CSV files
  csv_files <- c(
    "model_comparison.csv",
    "coefficients.csv", 
    "effect_sizes.csv",
    "model_averaged_effects.csv",
    "data_summary.csv",
    "convergence_info.csv",
    "model_formulas.csv"
  )
  
  for (file in csv_files) {
    file_path <- file.path(analysis_dir, file)
    if (file.exists(file_path)) {
      name <- gsub("\\.csv$", "", file)
      results[[name]] <- read.csv(file_path, stringsAsFactors = FALSE)
    }
  }
  
  # Load RDS files
  rds_files <- c(
    "fitted_models.rds",
    "bootstrap_results.rds",
    "prepared_data.rds",
    "complete_analysis_results.rds"
  )
  
  for (file in rds_files) {
    file_path <- file.path(analysis_dir, file)
    if (file.exists(file_path)) {
      name <- gsub("\\.rds$", "", file)
      results[[name]] <- readRDS(file_path)
    }
  }
  
  return(results)
}

#' Load batch analysis results
#' 
#' @param batch_dir Path to batch analysis directory
#' @return List with all results
load_batch_results <- function(batch_dir) {
  
  # Load main results
  all_results_file <- file.path(batch_dir, "all_results.rds")
  if (!file.exists(all_results_file)) {
    stop("Cannot find all_results.rds in batch directory")
  }
  
  all_results <- readRDS(all_results_file)
  
  # Also load integrated results if available
  integrated_dir <- file.path(batch_dir, "integrated_results")
  if (dir.exists(integrated_dir)) {
    integrated <- list(
      comparisons = read.csv(file.path(integrated_dir, "all_model_comparisons.csv")),
      coefficients = read.csv(file.path(integrated_dir, "all_coefficients.csv")),
      effects = read.csv(file.path(integrated_dir, "effect_comparison.csv")),
      summary = read.csv(file.path(integrated_dir, "analysis_summary_stats.csv"))
    )
  } else {
    integrated <- NULL
  }
  
  return(list(
    results = all_results,
    integration = integrated,
    output_dir = batch_dir
  ))
}

#' Create custom forest plot from saved results
#' 
#' @param coefficients Coefficients data frame
#' @param parameter_name Parameter to plot
#' @param model_filter Optional model name filter
#' @return ggplot object
create_custom_forest_plot <- function(coefficients, 
                                     parameter_name,
                                     model_filter = NULL) {
  
  # Filter data
  plot_data <- coefficients[coefficients$Parameter == parameter_name, ]
  
  if (!is.null(model_filter)) {
    plot_data <- plot_data[plot_data$Model %in% model_filter, ]
  }
  
  # Calculate odds ratios
  plot_data$OddsRatio <- exp(plot_data$Estimate)
  plot_data$OR_CI_lower <- exp(plot_data$Estimate - 1.96 * plot_data$StdErr)
  plot_data$OR_CI_upper <- exp(plot_data$Estimate + 1.96 * plot_data$StdErr)
  
  # Use bootstrap CIs if available
  if ("CI_lower" %in% names(plot_data) && !all(is.na(plot_data$CI_lower))) {
    plot_data$OR_CI_lower <- ifelse(!is.na(plot_data$CI_lower),
                                     exp(plot_data$CI_lower),
                                     plot_data$OR_CI_lower)
    plot_data$OR_CI_upper <- ifelse(!is.na(plot_data$CI_upper),
                                     exp(plot_data$CI_upper),
                                     plot_data$OR_CI_upper)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = OddsRatio, y = reorder(Model, OddsRatio))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), 
                   height = 0.2, size = 0.5) +
    geom_point(aes(color = p_value < 0.05), size = 3) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "gray50"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "") +
    labs(x = "Odds Ratio (95% CI)",
         y = "",
         title = paste("Forest Plot:", parameter_name)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Create comparison table from saved results
#' 
#' @param batch_results Loaded batch results
#' @param parameter_name Parameter to compare
#' @return Data frame with comparison
create_comparison_table <- function(batch_results, parameter_name) {
  
  comparison <- data.frame()
  
  for (name in names(batch_results$results)) {
    result <- batch_results$results[[name]]
    
    if (result$success && !is.null(result$effects)) {
      # Get best model effects
      best_model <- result$comparison$comparison$Model[1]
      param_effects <- result$effects[
        result$effects$Model == best_model & 
        result$effects$Parameter == parameter_name, 
      ]
      
      if (nrow(param_effects) > 0) {
        param_effects$Analysis <- name
        param_effects$N_species <- result$prepared_data$n_species
        param_effects$Best_AIC <- result$comparison$comparison$AIC[1]
        comparison <- rbind(comparison, param_effects[1, ])
      }
    }
  }
  
  # Add significance indicator
  comparison$Significant <- ifelse(comparison$p_value < 0.05, "*", "")
  
  # Sort by odds ratio
  comparison <- comparison[order(comparison$OddsRatio), ]
  
  return(comparison)
}

# Example usage:
if (interactive()) {
  cat("PhyloGLM Results Reloader\n")
  cat("========================\n\n")
  
  cat("Example 1: Load and plot single analysis results\n")
  cat("------------------------------------------------\n")
  cat('results <- load_single_analysis("Outputs/PhyloglmResults/test_run")\n')
  cat('p <- create_custom_forest_plot(results$coefficients, "HighConfidence_Coop")\n')
  cat('ggsave("custom_forest_plot.png", p, width = 8, height = 6)\n\n')
  
  cat("Example 2: Load and work with batch results\n")
  cat("-------------------------------------------\n")
  cat('batch <- load_batch_results("Outputs/PhyloglmResults/PhyloGLM_Batch_[timestamp]")\n')
  cat('comparison <- create_comparison_table(batch, "HighConfidence_Coop")\n')
  cat('write.csv(comparison, "cb_effects_comparison.csv", row.names = FALSE)\n\n')
  
  cat("Example 3: Create custom visualizations\n")
  cat("---------------------------------------\n")
  cat('# Modify and recreate any plot\n')
  cat('p <- create_comparative_forest_plot(batch, "HighConfidence_Coop") +\n')
  cat('  theme(text = element_text(size = 14)) +\n')
  cat('  labs(title = "My Custom Title")\n')
  cat('ggsave("modified_forest_plot.png", p, width = 10, height = 8)\n')
}