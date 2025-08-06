# Simple interface for running stepwise model expansion
# This script provides an easy way to expand your best phyloglm models

#' Run stepwise expansion with simple interface
#' 
#' This function takes your existing phyloglm results and tests adding
#' additional predictors to see if they improve the model
#' 
#' @param include_familial Include Familial Living predictor (has fewer species)
#' @param bootstrap_n Number of bootstrap iterations for final model (0 = skip)
#' @return List of expansion results
run_stepwise_simple <- function(include_familial = FALSE, 
                               bootstrap_n = 100) {
  
  cat("PhyloGLM Stepwise Model Expansion\n")
  cat("=================================\n\n")
  
  # Load required functions
  source("phyloglm_framework/stepwise_model_expansion.R")
  source("phyloglm_framework/stepwise_visualization.R")
  
  # Check if combined results exist
  results_path <- "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds"
  if (!file.exists(results_path)) {
    # Try alternate paths
    alt_paths <- c(
      "Outputs/PhyloglmResults/combined_all_analyses/all_results.rds",
      "Outputs/PhyloglmResults/combined_FS_response/all_results.rds"
    )
    
    for (path in alt_paths) {
      if (file.exists(path)) {
        results_path <- path
        break
      }
    }
    
    if (!file.exists(results_path)) {
      stop("Could not find combined results file. Please run batch analyses first.")
    }
  }
  
  cat("Using results from:", results_path, "\n\n")
  
  # Define candidate predictors
  candidates <- get_candidate_predictors()
  
  # Remove familial living if not requested
  if (!include_familial) {
    candidates$binary$FamilialLiving <- NULL
    cat("Note: Excluding Familial Living predictor (use include_familial=TRUE to include)\n\n")
  }
  
  # Run expansion
  expansion_results <- run_stepwise_batch(
    results_file = results_path,
    data_path = "Data_R_2025-06-09.csv",
    tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
    target_analyses = c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass"),
    output_dir = paste0("Outputs/PhyloglmResults/stepwise_expansion_",
                       format(Sys.Date(), "%Y%m%d")),
    candidates = candidates,
    aic_threshold = 2,
    n_bootstrap = bootstrap_n,
    verbose = TRUE
  )
  
  # Create and save visualizations
  cat("\n\nCreating summary visualizations...\n")
  
  output_dir <- paste0("Outputs/PhyloglmResults/stepwise_expansion_",
                      format(Sys.Date(), "%Y%m%d"))
  
  # Summary figure
  create_stepwise_summary_figure(
    expansion_results,
    output_file = file.path(output_dir, "stepwise_summary")
  )
  
  # Model comparison table
  comparison <- create_model_comparison_table(
    expansion_results,
    output_file = file.path(output_dir, "model_comparison.csv")
  )
  
  # Predictor importance
  if (length(expansion_results) > 1) {
    importance_plot <- plot_predictor_importance(expansion_results)
    ggsave(file.path(output_dir, "predictor_importance.png"),
           importance_plot, width = 10, height = 8, dpi = 300)
  }
  
  # Print summary
  cat("\n\nSUMMARY OF EXPANSIONS\n")
  cat("=====================\n\n")
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    
    cat(name, "\n")
    cat(rep("-", nchar(name)), "\n", sep = "")
    
    if (result$n_predictors_added > 0) {
      cat("Added", result$n_predictors_added, "predictor(s):\n")
      for (step in result$expansion_history) {
        cat("  -", step$added, "(AIC improved by", 
            round(step$improvement, 2), ")\n")
      }
      cat("Total AIC improvement:", round(result$total_improvement, 2), "\n")
    } else {
      cat("No predictors improved the model\n")
    }
    cat("\n")
  }
  
  cat("Results saved to:", output_dir, "\n\n")
  
  return(invisible(expansion_results))
}

# Example usage:
# Run with default settings (excluding Familial Living)
# results <- run_stepwise_simple()

# Run including Familial Living predictor
# results <- run_stepwise_simple(include_familial = TRUE)

# Run without bootstrap (faster)
# results <- run_stepwise_simple(bootstrap_n = 0)