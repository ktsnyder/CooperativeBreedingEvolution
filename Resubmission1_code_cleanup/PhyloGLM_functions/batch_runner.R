# Batch Analysis Runner for PhyloGLM Framework
# Runs multiple phyloglm analyses with robust error handling
# 
# 9/9/25 - Removed attempted parallel processing functionality

library(phylolm)
#source("phyloglm_framework/variable_classification.R")
source(file.path("PhyloGLM_functions", "formula_builder_comprehensive.R"))
source(file.path("PhyloGLM_functions", "data_preparation_simple.R"))
source(file.path("PhyloGLM_functions", "batch_runner_helpers.R"))

#' Run a batch of phyloglm analyses
#' 
#' @param analysis_configs List of analysis configurations
#' @param data Original data frame
#' @param tree Original phylogenetic tree
#' @param output_dir Base directory for outputs
#' @param bootstrap_n Number of bootstrap replicates for top models
#' @param save_intermediate Whether to save intermediate results
#' @return List of analysis results
run_phyloglm_batch <- function(analysis_configs, 
                              data, 
                              tree,
                              output_dir,
                              bootstrap_n = 100,
                              save_intermediate = TRUE) {
  
  # Create output directory with timestamp and bootstrap info
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  batch_dir <- file.path(output_dir, paste0("PhyloGLM_Batch_", timestamp, "_boot", bootstrap_n))
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save configurations
  config_dir <- file.path(batch_dir, "configs")
  dir.create(config_dir, showWarnings = FALSE)
  saveRDS(analysis_configs, file.path(config_dir, "analysis_configs.rds"))
  
  # Create metadata file
  metadata_file <- file.path(batch_dir, "analysis_metadata.txt")
  write_analysis_metadata(
    file = metadata_file,
    bootstrap_n = bootstrap_n,
    n_analyses = length(analysis_configs),
    data_file = deparse(substitute(data)),
    tree_file = deparse(substitute(tree))
  )
  
  # Initialize results storage
  all_results <- list()
  
  # Progress tracking
  n_analyses <- length(analysis_configs)
  cat(paste("Starting batch analysis with", n_analyses, "configurations\n"))
  
  # Run each analysis
  for (i in seq_along(analysis_configs)) {
    config <- analysis_configs[[i]]
    
    # Create analysis name
    analysis_name <- create_analysis_name(config)
    cat(paste("\n[", i, "/", n_analyses, "] Running:", analysis_name, "\n"))
    
    # Create analysis-specific directory
    analysis_dir <- file.path(batch_dir, "individual_analyses", analysis_name)
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Run analysis with error handling
    result <- tryCatch({
      run_single_analysis(
        config = config,
        data = data,
        tree = tree,
        output_dir = analysis_dir,
        bootstrap_n = bootstrap_n,
        save_outputs = save_intermediate
      )
    }, error = function(e) {
      cat(paste("ERROR in analysis", analysis_name, ":", e$message, "\n"))
      list(
        success = FALSE,
        error = e$message,
        config = config,
        analysis_name = analysis_name
      )
    })
    
    # Store result
    all_results[[analysis_name]] <- result
    
    # Save intermediate results if requested
    if (save_intermediate) {
      saveRDS(result, file.path(analysis_dir, "analysis_result.rds"))
    }
  }
  
  # Create integrated results
  cat("\nCreating integrated results...\n")
  integrated_dir <- file.path(batch_dir, "integrated_results")
  dir.create(integrated_dir, showWarnings = FALSE)
  
  integration_result <- integrate_batch_results(
    all_results, 
    output_dir = integrated_dir
  )
  
  # Save complete results
  saveRDS(all_results, file.path(batch_dir, "all_results.rds"))
  
  # Create summary report
  create_batch_summary_report(
    all_results,
    output_file = file.path(batch_dir, "batch_summary.txt")
  )
  
  cat("\nBatch analysis complete. Results saved to:", batch_dir, "\n")
  
  return(list(
    results = all_results,
    integration = integration_result,
    output_dir = batch_dir,
    bootstrap_n = bootstrap_n
  ))
}

#' Write analysis metadata file
#' 
#' @param file Path to metadata file
#' @param bootstrap_n Number of bootstrap iterations
#' @param n_analyses Number of analyses
#' @param data_file Data file name
#' @param tree_file Tree file name
write_analysis_metadata <- function(file, bootstrap_n, n_analyses, 
                                   data_file, tree_file) {
  
  # Get package versions
  pkg_versions <- c(
    phylolm = as.character(packageVersion("phylolm")),
    ape = as.character(packageVersion("ape")),
    ggplot2 = as.character(packageVersion("ggplot2")),
    dplyr = as.character(packageVersion("dplyr"))
  )
  
  # Write metadata
  con <- file(file, open = "w")
  
  writeLines("PhyloGLM Analysis Metadata", con)
  writeLines("=========================", con)
  writeLines(paste("Date:", Sys.Date()), con)
  writeLines(paste("Time:", format(Sys.time(), "%H:%M:%S")), con)
  writeLines(paste("R Version:", R.version.string), con)
  writeLines("", con)
  
  writeLines("Analysis Parameters:", con)
  writeLines(paste("  Bootstrap iterations:", bootstrap_n), con)
  writeLines(paste("  Number of analyses:", n_analyses), con)
  writeLines("", con)
  
  writeLines("Input Files:", con)
  writeLines(paste("  Data file:", data_file), con)
  writeLines(paste("  Tree file:", tree_file), con)
  writeLines("", con)
  
  writeLines("Package Versions:", con)
  for (pkg in names(pkg_versions)) {
    writeLines(paste("  ", pkg, ":", pkg_versions[pkg]), con)
  }
  writeLines("", con)
  
  writeLines("System Information:", con)
  writeLines(paste("  Platform:", Sys.info()["sysname"]), con)
  writeLines(paste("  Machine:", Sys.info()["machine"]), con)
  
  close(con)
}

#' Run a single phyloglm analysis
#' 
#' @param config Analysis configuration
#' @param data Original data
#' @param tree Original tree
#' @param output_dir Output directory for this analysis
#' @param bootstrap_n Number of bootstrap replicates
#' @param save_outputs Whether to save outputs
#' @return Analysis results
run_single_analysis <- function(config, data, tree, output_dir, 
                               bootstrap_n = 100, save_outputs = TRUE) {
  
  # Create output directory if it doesn't exist
  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Prepare data
  cat("  Preparing data...\n")
  prepared <- prepare_analysis_data(data, tree, config)
  
  # Check data quality
  quality <- check_data_quality(prepared)
  if (!quality$is_acceptable) {
    stop(paste("Data quality issues:", 
               paste(names(quality$issues), collapse = ", ")))
  }
  
  # Save data summary
  if (save_outputs) {
    data_summary <- create_data_summary(prepared)
    write.csv(data_summary, 
              file.path(output_dir, "data_summary.csv"),
              row.names = FALSE)
  }
  
  # Build model formulas
  cat("  Building model set...\n")
  
  # Use 15-model builder for standard analyses (2 predictors + 1 control)
  if (length(config$predictors) == 2 && length(config$controls) == 1) {
    formulas <- build_comprehensive_15_models(
      response_var = config$response,
      pred_var = config$predictors[1],
      terr_var = config$predictors[2], 
      control_var = config$controls[1]
    )
  } else {
    # Use comprehensive builder for other cases
    formulas <- build_comprehensive_model_set(config, prepared$data)
  }
  
  cat(paste("  ", length(formulas), "models to fit\n"))
  
  # Fit all models
  cat("  Fitting models...\n")
  model_results <- fit_all_models(
    formulas = formulas,
    data = prepared$data,
    tree = prepared$tree,
    config = config
  )
  
  # Model comparison
  cat("  Comparing models...\n")
  comparison <- compare_models(model_results)
  
  # Identify top models (deltaAIC < 2)
  top_model_names <- comparison$comparison$Model[comparison$comparison$deltaAIC < 2]
  top_models <- comparison$models[top_model_names]
  cat(paste("  Found", length(top_models), "top models (ΔAIC < 2)\n"))
  
  # Refit top models with bootstrap
  if (bootstrap_n > 0 && length(top_models) > 0) {
    cat(paste("  Bootstrapping top models (", bootstrap_n, "replicates)...\n"))
    boot_results <- bootstrap_top_models(
      top_models = top_models,
      data = prepared$data,
      tree = prepared$tree,
      n_boot = bootstrap_n
    )
  } else {
    boot_results <- NULL
  }
  
  # Extract results
  cat("  Extracting coefficients and effects...\n")
  coefficients <- extract_all_coefficients(model_results, boot_results)
  effects <- calculate_effect_sizes(coefficients)
  
  # Create best model effects summary
  best_effects <- create_best_model_effects(effects, comparison, boot_results)
  
  # Save outputs
  if (save_outputs) {
    cat("  Saving outputs...\n")
    
    # Model comparison with bootstrap info
    comparison_with_info <- comparison$comparison
    comparison_with_info$bootstrap_n <- bootstrap_n
    write.csv(comparison_with_info, 
              file.path(output_dir, "model_comparison.csv"),
              row.names = FALSE)
    
    # Coefficients
    write.csv(coefficients, 
              file.path(output_dir, "coefficients.csv"),
              row.names = FALSE)
    
    # Effect sizes
    write.csv(effects,
              file.path(output_dir, "effect_sizes.csv"),
              row.names = FALSE)
    
    # Best model effects
    if (nrow(best_effects) > 0) {
      write.csv(best_effects,
                file.path(output_dir, "best_model_effects.csv"),
                row.names = FALSE)
    }
    
    # Save models
    saveRDS(model_results, file.path(output_dir, "fitted_models.rds"))
    
    # Save bootstrap results if available
    if (!is.null(boot_results)) {
      saveRDS(boot_results, file.path(output_dir, "bootstrap_results.rds"))
    }
    
    # Save prepared data
    saveRDS(prepared, file.path(output_dir, "prepared_data.rds"))
    
    # Save convergence information
    convergence_df <- data.frame(
      Model = names(model_results$convergence),
      Converged = sapply(model_results$convergence, function(x) x$converged),
      stringsAsFactors = FALSE
    )
    write.csv(convergence_df, 
              file.path(output_dir, "convergence_info.csv"),
              row.names = FALSE)
    
    # Save model formulas
    formula_df <- data.frame(
      Model = names(model_results$formulas),
      Formula = sapply(model_results$formulas, function(f) as.character(f)[3]),
      stringsAsFactors = FALSE
    )
    write.csv(formula_df,
              file.path(output_dir, "model_formulas.csv"),
              row.names = FALSE)
    
    # Save complete results object
    complete_results <- list(
      config = config,
      prepared_data = prepared,
      quality = quality,
      models = model_results,
      comparison = comparison,
      coefficients = coefficients,
      effects = effects,
      best_effects = best_effects,
      bootstrap = boot_results,
      bootstrap_n = bootstrap_n,
      analysis_date = Sys.Date()
    )
    saveRDS(complete_results, file.path(output_dir, "complete_analysis_results.rds"))
  }
  
  # Return results
  return(list(
    success = TRUE,
    config = config,
    prepared_data = prepared,
    quality = quality,
    models = model_results,
    comparison = comparison,
    coefficients = coefficients,
    effects = effects,
    best_effects = best_effects,
    bootstrap = boot_results,
    output_dir = output_dir
  ))
}

#' Fit all models in a formula list
#' 
#' @param formulas List of formulas
#' @param data Prepared data
#' @param tree Prepared tree
#' @param config Analysis configuration
#' @return List of fitted models and metadata
fit_all_models <- function(formulas, data, tree, config) {
  
  models <- list()
  fit_times <- list()
  convergence <- list()
  
  for (i in seq_along(formulas)) {
    model_name <- names(formulas)[i]
    formula <- formulas[[i]]
    
    # Time the fit
    start_time <- Sys.time()
    
    # Fit model with error handling
    fit_result <- tryCatch({
      phyloglm(
        formula = formula,
        data = data,
        phy = tree,
        method = "logistic_MPLE",
        btol = 30,  # Increased tolerance
        log.alpha.bound = 4,
        start.beta = NULL,
        start.alpha = NULL,
        boot = 0  # No bootstrap in initial fit
      )
    }, error = function(e) {
      list(
        converged = FALSE,
        error = e$message
      )
    }, warning = function(w) {
      # Still try to return the model even with warnings
      fit <- phyloglm(
        formula = formula,
        data = data,
        phy = tree,
        method = "logistic_MPLE",
        btol = 30,
        log.alpha.bound = 4,
        boot = 0
      )
      fit$warning <- w$message
      fit
    })
    
    end_time <- Sys.time()
    fit_times[[model_name]] <- as.numeric(end_time - start_time, units = "secs")
    
    # Store model if successful
    if (!is.null(fit_result) && !is.null(fit_result$coefficients)) {
      models[[model_name]] <- fit_result
      convergence[[model_name]] <- list(
        converged = ifelse(is.null(fit_result$converged), TRUE, fit_result$converged),
        iterations = fit_result$iterations,
        log_likelihood = logLik(fit_result),
        alpha = fit_result$alpha
      )
    } else {
      convergence[[model_name]] <- list(
        converged = FALSE,
        error = ifelse(is.null(fit_result$error), "Unknown error", fit_result$error)
      )
    }
  }
  
  # Remove failed models
  successful_models <- models[!sapply(models, is.null)]
  
  cat(paste("    Successfully fit", length(successful_models), "of", 
            length(formulas), "models\n"))
  
  return(list(
    models = successful_models,
    formulas = formulas,
    convergence = convergence,
    fit_times = fit_times,
    n_failed = length(formulas) - length(successful_models)
  ))
}

#' Compare fitted models
#' 
#' @param model_results Output from fit_all_models
#' @return Model comparison data frame
compare_models <- function(model_results) {
  
  models <- model_results$models
  
  # Extract AIC values
  aic_values <- sapply(models, AIC)
  
  # Create comparison table
  comparison <- data.frame(
    Model = names(models),
    AIC = aic_values,
    stringsAsFactors = FALSE
  )
  
  # Add other metrics
  comparison$LogLik <- sapply(models, logLik)
  comparison$df <- sapply(models, function(m) length(coef(m)))
  
  # Calculate deltaAIC and weights
  comparison <- comparison[order(comparison$AIC), ]
  comparison$deltaAIC <- comparison$AIC - min(comparison$AIC)
  comparison$weight <- exp(-0.5 * comparison$deltaAIC)
  comparison$weight <- comparison$weight / sum(comparison$weight)
  comparison$cum_weight <- cumsum(comparison$weight)
  
  # Add R-squared approximations
  comparison$R2_McFadden <- sapply(models[comparison$Model], function(m) {
    null_ll <- logLik(models[["Null"]])
    1 - (logLik(m) / null_ll)
  })
  
  return(list(
    comparison = comparison,
    models = models,
    best_model = models[[comparison$Model[1]]]
  ))
}

#' Create analysis name from configuration
#' 
#' @param config Analysis configuration
#' @return Character string for directory/file naming
create_analysis_name <- function(config) {
  if (!is.null(config$name)) {
    return(config$name)
  }
  
  response_short <- simplify_var_name(config$response)
  
  predictor_names <- sapply(config$predictors, simplify_var_name)
  predictor_str <- paste(predictor_names, collapse = "_")
  
  if (!is.null(config$controls) && length(config$controls) > 0) {
    control_names <- sapply(config$controls, simplify_var_name)
    control_str <- paste0("_ctrl", paste(control_names, collapse = ""))
  } else {
    control_str <- ""
  }
  
  return(paste0(response_short, "_vs_", predictor_str, control_str))
}

#' Simplify variable name for file/directory naming
#' 
#' @param var_name Full variable name
#' @return Simplified name
simplify_var_name <- function(var_name) {
  # Remove common suffixes and special characters
  simple <- gsub("_Agg01|_AVONET|_Jetz|2017", "", var_name)
  simple <- gsub("HighConfidence_", "", simple)
  simple <- gsub("TerritorialityWeakVsStrong", "TerrWS", simple)
  simple <- gsub("logMass_normalized", "Mass", simple)
  simple <- gsub("logMass", "Mass", simple)
  simple <- gsub("PercentAbsLog", "", simple)
  simple <- gsub("Dimorphism", "Dim", simple)
  simple <- gsub("GeographicRegion", "Region", simple)
  simple <- gsub("Migration", "Migr", simple)
  simple <- gsub("Centroid.Latitude", "Lat", simple)
  simple <- gsub("FemaleSong", "FS", simple)
  simple <- gsub("Coop", "CB", simple)
  simple <- gsub("FamilialLiving", "Fam", simple)
  
  # Remove any remaining special characters
  simple <- gsub("[^A-Za-z0-9]", "", simple)
  
  return(simple)
}

#' Create batch summary report
#' 
#' @param all_results List of all analysis results
#' @param output_file Path to output text file
create_batch_summary_report <- function(all_results, output_file) {
  
  # Open connection
  con <- file(output_file, open = "w")
  
  # Header
  writeLines("PhyloGLM Batch Analysis Summary", con)
  writeLines(paste("Generated:", Sys.time()), con)
  writeLines(paste("\nTotal analyses:", length(all_results)), con)
  
  # Add bootstrap info if available
  if (exists("bootstrap_n", where = parent.frame())) {
    writeLines(paste("Bootstrap iterations:", get("bootstrap_n", parent.frame())), con)
  }
  
  # Success/failure summary
  n_success <- sum(sapply(all_results, function(r) r$success))
  n_failed <- length(all_results) - n_success
  writeLines(paste("Successful:", n_success), con)
  writeLines(paste("Failed:", n_failed), con)
  
  # Individual analysis summaries
  writeLines("\n\nIndividual Analysis Summaries:", con)
  writeLines(strrep("-", 60), con)
  
  for (name in names(all_results)) {
    result <- all_results[[name]]
    writeLines(paste("\nAnalysis:", name), con)
    
    if (result$success) {
      writeLines(paste("  Status: Success"), con)
      writeLines(paste("  N species:", result$prepared_data$n_species), con)
      writeLines(paste("  N models:", length(result$models$models)), con)
      
      # Best model
      if (!is.null(result$comparison)) {
        best <- result$comparison$comparison[1, ]
        writeLines(paste("  Best model:", best$Model), con)
        writeLines(paste("  Best AIC:", round(best$AIC, 2)), con)
      }
    } else {
      writeLines(paste("  Status: Failed"), con)
      writeLines(paste("  Error:", result$error), con)
    }
  }
  
  close(con)
}