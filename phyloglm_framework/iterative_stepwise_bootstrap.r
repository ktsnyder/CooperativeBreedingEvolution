# Iterative Stepwise Model Expansion with Bootstrap
# This script implements true stepwise forward selection with bootstrap CIs

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)

# Source the bootstrap functions from the single-parameter script
# source("phyloglm_framework/stepwise_expansion_bootstrap_updated.R")

# Main function for iterative stepwise selection
run_iterative_stepwise_bootstrap <- function(
  single_param_results_dir,  # e.g., "Outputs/PhyloglmResults/stepwise_boot500_20250716"
  output_subdir = "Iterative_Stepwise_Model_Tests",
  aic_threshold = 2,
  n_bootstrap = NULL,  # If NULL, detect from input dir name
  max_iterations = 7,  # Safety limit
  save_boot_matrices = TRUE,
  verbose = TRUE
) {
  
  # Extract n_bootstrap from directory name if not provided
  if (is.null(n_bootstrap)) {
    dir_name <- basename(single_param_results_dir)
    boot_match <- regmatches(dir_name, regexpr("boot[0-9]+", dir_name))
    if (length(boot_match) > 0) {
      n_bootstrap <- as.numeric(sub("boot", "", boot_match))
      message("Detected n_bootstrap = ", n_bootstrap, " from directory name")
    } else {
      stop("Could not detect n_bootstrap from directory name. Please specify explicitly.")
    }
  }
  
  # Create output directory structure
  output_dir <- file.path(single_param_results_dir, output_subdir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create analysis-specific subdirectories with bootstrap count
  for (analysis in c("CB_vs_FS_TerrWS_Mass", "FS_vs_CB_TerrWS_Mass")) {
    analysis_dir <- file.path(output_dir, paste0(analysis, "_boot", n_bootstrap))
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
    if (save_boot_matrices) {
      dir.create(file.path(analysis_dir, "bootstrap_matrices"), recursive = TRUE, showWarnings = FALSE)
    }
    dir.create(file.path(analysis_dir, "iteration_summaries"), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Load necessary data and results
  message("\nLoading data and single-parameter results...")
  
  # Load single-parameter results
  results_file <- list.files(single_param_results_dir, pattern = "^expansion_results_boot.*\\.rds$", full.names = TRUE)[1]
  if (!file.exists(results_file)) {
    stop("Could not find results file: ", results_file)
  }
  single_param_results <- readRDS(results_file)
  
  # Load original data files (get paths from the single param results directory)
  # First, find the parent directory that contains the data files
  parent_dir <- dirname(dirname(dirname(single_param_results_dir)))
  data_path <- file.path(parent_dir, "Data_R_2025-06-09.csv")
  tree_path <- file.path(parent_dir, "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  
  if (!file.exists(data_path)) {
    stop("Could not find data file: ", data_path)
  }
  if (!file.exists(tree_path)) {
    stop("Could not find tree file: ", tree_path)
  }
  
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Initialize comprehensive results table
  comprehensive_results <- data.frame()
  
  # Initialize iteration paths table
  iteration_paths <- data.frame()
  
  # Define available predictors
  available_predictors <- list(
    Territory = list(
      name = "Territory (1-3)",
      var = "Territory_num",
      values_func = function(data) as.numeric(data$Territory),
      skip_conditions = "Terr3|Territory_num"
    ),
    Migration = list(
      name = "Migration (1-3)",
      var = "Migration_num", 
      values_func = function(data) as.numeric(data$Migration_AVONET),
      skip_conditions = "Migration_num"
    ),
    Latitude = list(
      name = "Absolute Latitude",
      var = "abs_Latitude",
      values_func = function(data) abs(data$Centroid.Latitude_AVONET),
      skip_conditions = "abs_Latitude|absLat"
    ),
    Wing = list(
      name = "Wing Dimorphism",
      var = "PercentAbsLogWingDimorphism",
      values_func = function(data) data$PercentAbsLogWingDimorphism,
      skip_conditions = "WingDimorphism|PercentAbsLogWingDimorphism"
    ),
    Plumage = list(
      name = "Plumage Dimorphism",
      var = "logMaleFemalePlumageDiffAbs",
      values_func = function(data) log(data$MaleFemalePlumageDiffAbs + 1),
      skip_conditions = "PlumageDiff|logMaleFemalePlumageDiff"
    ),
    FamilialLiving = list(
      name = "Familial Living",
      var = "Griesser2017FamilialLiving",
      values_func = function(data) data$Griesser2017FamilialLiving,
      skip_conditions = "FamilialLiving|Griesser2017FamilialLiving"
    ),
    Region = list(
      name = "Geographic Region",
      var = "GeographicRegion_Jetz",
      values_func = function(data) as.factor(data$GeographicRegion_Jetz),
      skip_conditions = "GeographicRegion|GeographicRegion_Jetz"
    )
  )
  
  # Process each analysis
  for (analysis_name in c("CB_vs_FS_TerrWS_Mass", "FS_vs_CB_TerrWS_Mass")) {
    
    message("\n", paste(rep("=", 60), collapse=""))
    message("ITERATIVE STEPWISE SELECTION: ", analysis_name)
    message(paste(rep("=", 60), collapse=""), "\n")
    
    # Get base model info from single-parameter results
    base_info <- single_param_results[[analysis_name]]
    if (is.null(base_info)) {
      message("No results found for ", analysis_name, ". Skipping.")
      next
    }
    
    # Initialize tracking for this analysis
    current_formula <- formula(base_info$base_model)
    current_data <- base_info$base_data
    current_tree <- base_info$base_tree
    current_aic <- base_info$original_base_aic
    base_formula_str <- paste(deparse(current_formula), collapse = " ")
    added_predictors <- character()
    
    message("Starting base model: ", base_formula_str)
    message("Base AIC: ", round(current_aic, 2))
    message("N species: ", nrow(current_data))
    
    # ITERATION 1: Use results from single-parameter analysis
    message("\n--- ITERATION 1 ---")
    message("Using results from single-parameter analysis...")
    
    # Get improvements from single-parameter analysis, sorted by AIC improvement
    iter1_improvements <- base_info$improvements %>%
      filter(Improved_Model == TRUE) %>%  # Only those that improved by threshold
      arrange(desc(AIC_Improvement))
    
    if (nrow(iter1_improvements) == 0) {
      message("No predictors improved model by threshold. Analysis complete.")
      next
    }
    
    # Process all tested models for iteration 1 (including non-improvements)
    all_iter1_results <- base_info$improvements
    
    for (i in 1:nrow(all_iter1_results)) {
      pred_row <- all_iter1_results[i,]
      pred_result <- base_info$predictor_results[[pred_row$Predictor_Name]]
      
      if (!is.null(pred_result)) {
        # Add all coefficients to comprehensive results
        comp_rows <- process_model_coefficients(
          model_fit = pred_result$new_bootstrap$fit,
          bootstrap_results = pred_result$new_bootstrap,
          analysis_name = analysis_name,
          iteration = 1,
          predictor_added = pred_row$Predictor_Name,
          predictor_variable = pred_row$Predictor_Variable,
          base_formula = base_formula_str,
          previous_aic = current_aic,
          expanded_aic = pred_row$Expanded_Model_AIC,
          aic_improvement = pred_row$AIC_Improvement,
          model_selected = pred_row$Improved_Model,
          n_species = pred_row$N_Species
        )
        
        comprehensive_results <- rbind(comprehensive_results, comp_rows)
      }
    }
    
    # Select best predictor from iteration 1
    best_pred <- iter1_improvements[1,]
    best_pred_result <- base_info$predictor_results[[best_pred$Predictor_Name]]
    
    message("\nBest predictor: ", best_pred$Predictor_Name)
    message("AIC improvement: ", round(best_pred$AIC_Improvement, 2))
    
    # Update current model state
    current_formula <- best_pred_result$formula
    current_data <- best_pred_result$data
    current_tree <- best_pred_result$tree
    current_aic <- best_pred$Expanded_Model_AIC
    added_predictors <- c(added_predictors, best_pred$Predictor_Name)
    
    # Add to iteration paths
    iteration_paths <- rbind(iteration_paths, data.frame(
      Analysis_Name = analysis_name,
      Iteration = 1,
      Predictor_Added = best_pred$Predictor_Name,
      Predictor_Variable = best_pred$Predictor_Variable,
      Previous_AIC = base_info$original_base_aic,
      Current_AIC = current_aic,
      AIC_Improvement = best_pred$AIC_Improvement,
      Cumulative_AIC_Improvement = base_info$original_base_aic - current_aic,
      N_Species = nrow(current_data),
      Current_Formula = paste(deparse(current_formula), collapse = " "),
      stringsAsFactors = FALSE
    ))
    
    # ITERATIONS 2+: Test remaining predictors
    for (iteration in 2:max_iterations) {
      
      message("\n--- ITERATION ", iteration, " ---")
      
      # Determine which predictors remain to be tested
      remaining_predictors <- character()
      for (pred_name in names(available_predictors)) {
        pred <- available_predictors[[pred_name]]
        # Check if already in model
        formula_str <- paste(deparse(current_formula), collapse = " ")
        if (!grepl(pred$skip_conditions, formula_str)) {
          remaining_predictors <- c(remaining_predictors, pred_name)
        }
      }
      
      if (length(remaining_predictors) == 0) {
        message("No more predictors to test. Analysis complete.")
        break
      }
      
      message("Testing ", length(remaining_predictors), " remaining predictors...")
      
      # Test each remaining predictor
      iteration_results <- list()
      
      for (pred_id in remaining_predictors) {
        pred <- available_predictors[[pred_id]]
        
        message("\nTesting: ", pred$name)
        
        # Test this predictor
          test_result <- test_predictor_iterative(
            base_formula = current_formula,
            base_data = current_data,
            base_tree = current_tree,
            base_aic = current_aic,
            predictor_name = pred$name,
            predictor_var = pred$var,
            predictor_values = pred$values_func(full_data),
            full_data = full_data,
            n_bootstrap = n_bootstrap,
            analysis_name = analysis_name,
            iteration = iteration,
            save_matrices = save_boot_matrices,
            matrix_dir = file.path(output_dir, paste0(analysis_name, "_boot", n_bootstrap), "bootstrap_matrices")
          )
        
        if (!is.null(test_result)) {
          iteration_results[[pred_id]] <- test_result
          
          # Add to comprehensive results
          comp_rows <- process_model_coefficients(
            model_fit = test_result$expanded_fit,
            bootstrap_results = test_result$expanded_bootstrap,
            analysis_name = analysis_name,
            iteration = iteration,
            predictor_added = pred$name,
            predictor_variable = pred$var,
            base_formula = paste(deparse(current_formula), collapse = " "),
            previous_aic = current_aic,
            expanded_aic = test_result$expanded_aic,
            aic_improvement = test_result$aic_improvement,
            model_selected = test_result$aic_improvement >= aic_threshold,
            n_species = test_result$n_species
          )
          
          comprehensive_results <- rbind(comprehensive_results, comp_rows)
        }
      }
      
      # Find best improvement
      best_improvement <- NULL
      best_improvement_value <- 0
      
      for (pred_id in names(iteration_results)) {
        result <- iteration_results[[pred_id]]
        if (result$aic_improvement >= aic_threshold && result$aic_improvement > best_improvement_value) {
          best_improvement <- result
          best_improvement_value <- result$aic_improvement
          best_improvement$predictor_id <- pred_id
        }
      }
      
      if (is.null(best_improvement)) {
        message("\nNo predictor improves model by threshold. Analysis complete.")
        break
      }
      
      # Update model with best predictor
      pred_info <- available_predictors[[best_improvement$predictor_id]]
      message("\nBest predictor: ", pred_info$name)
      message("AIC improvement: ", round(best_improvement$aic_improvement, 2))
      
      current_formula <- best_improvement$formula
      current_data <- best_improvement$data
      current_tree <- best_improvement$tree
      current_aic <- best_improvement$expanded_aic
      added_predictors <- c(added_predictors, pred_info$name)
      
      # Add to iteration paths
      iteration_paths <- rbind(iteration_paths, data.frame(
        Analysis_Name = analysis_name,
        Iteration = iteration,
        Predictor_Added = pred_info$name,
        Predictor_Variable = pred_info$var,
        Previous_AIC = best_improvement$base_aic,
        Current_AIC = current_aic,
        AIC_Improvement = best_improvement$aic_improvement,
        Cumulative_AIC_Improvement = base_info$original_base_aic - current_aic,
        N_Species = nrow(current_data),
        Current_Formula = paste(deparse(current_formula), collapse = " "),
        stringsAsFactors = FALSE
      ))
      
      # Save intermediate results
      saveRDS(list(
        comprehensive_results = comprehensive_results,
        iteration_paths = iteration_paths
      ), file.path(output_dir, paste0(analysis_name, "_boot", n_bootstrap), 
                   "iteration_summaries", 
                   paste0("iteration", iteration, "_results.rds")))
      
      # Also save iteration summary as CSV
      iter_summary <- comprehensive_results %>%
        filter(Analysis_Name == analysis_name, Iteration <= iteration)
      
      write.csv(iter_summary,
                file.path(output_dir, paste0(analysis_name, "_boot", n_bootstrap), 
                         "iteration_summaries",
                         paste0("iteration", iteration, "_comprehensive.csv")),
                row.names = FALSE)
    }
    
    # Save final results for this analysis
    analysis_summary <- list(
      analysis_name = analysis_name,
      base_model = base_formula_str,
      base_aic = base_info$original_base_aic,
      final_model = paste(deparse(current_formula), collapse = " "),
      final_aic = current_aic,
      total_improvement = base_info$original_base_aic - current_aic,
      predictors_added = added_predictors,
      n_iterations = length(added_predictors)
    )
    
    analysis_dir <- file.path(output_dir, paste0(analysis_name, "_boot", n_bootstrap))
    
    saveRDS(analysis_summary, 
            file.path(analysis_dir, paste0(analysis_name, "_boot", n_bootstrap, "_summary.rds")))
    
    # Also save analysis-specific CSV files
    analysis_comp_results <- comprehensive_results %>%
      filter(Analysis_Name == analysis_name)
    
    write.csv(analysis_comp_results,
              file.path(analysis_dir, paste0(analysis_name, "_boot", n_bootstrap, "_comprehensive_results.csv")),
              row.names = FALSE)
    
    analysis_iter_paths <- iteration_paths %>%
      filter(Analysis_Name == analysis_name)
    
    write.csv(analysis_iter_paths,
              file.path(analysis_dir, paste0(analysis_name, "_boot", n_bootstrap, "_iteration_paths.csv")),
              row.names = FALSE)
  }
  
  # Save comprehensive results - ENSURE THIS HAPPENS
  message("\nSaving comprehensive results CSV...")
  write.csv(comprehensive_results,
            file.path(output_dir, paste0("comprehensive_results_boot", n_bootstrap, ".csv")),
            row.names = FALSE)
  message("Saved: ", file.path(output_dir, paste0("comprehensive_results_boot", n_bootstrap, ".csv")))
  
  # Save iteration paths
  message("Saving iteration paths CSV...")
  write.csv(iteration_paths,
            file.path(output_dir, paste0("iteration_paths_boot", n_bootstrap, ".csv")),
            row.names = FALSE)
  message("Saved: ", file.path(output_dir, paste0("iteration_paths_boot", n_bootstrap, ".csv")))
  
  # Create final models comparison
  message("Creating final models comparison...")
  final_comparison <- create_final_comparison(single_param_results, iteration_paths, comprehensive_results)
  write.csv(final_comparison,
            file.path(output_dir, paste0("final_models_comparison_boot", n_bootstrap, ".csv")),
            row.names = FALSE)
  message("Saved: ", file.path(output_dir, paste0("final_models_comparison_boot", n_bootstrap, ".csv")))
  
  message("\n", paste(rep("=", 60), collapse=""))
  message("ITERATIVE STEPWISE SELECTION COMPLETE")
  message("Results saved to: ", output_dir)
  message(paste(rep("=", 60), collapse=""))
  
  return(list(
    comprehensive_results = comprehensive_results,
    iteration_paths = iteration_paths,
    final_comparison = final_comparison,
    output_dir = output_dir
  ))
}

# Function to test a predictor in the iterative context
test_predictor_iterative <- function(
  base_formula, base_data, base_tree, base_aic,
  predictor_name, predictor_var, predictor_values,
  full_data, n_bootstrap, analysis_name, iteration,
  save_matrices = TRUE, matrix_dir = NULL
) {
  
  # Add predictor to data
  test_data <- base_data
  test_data[[predictor_var]] <- predictor_values[match(test_data$species, full_data$species)]
  
  # Convert categorical variables to factors
  if (predictor_var == "GeographicRegion_Jetz") {
    test_data[[predictor_var]] <- as.factor(test_data[[predictor_var]])
  }
  
  # Remove NAs
  test_data <- test_data[!is.na(test_data[[predictor_var]]), ]
  
  if (nrow(test_data) < 100) {
    message("   Not enough species with data (", nrow(test_data), ")")
    return(NULL)
  }
  
  # Match tree
  test_tree <- keep.tip(base_tree, test_data$species)
  test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
  
  message("   N species with ", predictor_name, " data: ", nrow(test_data))
  
  # Refit base model to subset (without bootstrap for efficiency)
  base_subset_fit <- phyloglm(
    formula = base_formula,
    data = test_data,
    phy = test_tree,
    method = "logistic_MPLE",
    btol = 50,
    log.alpha.bound = 4
  )
  
  base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
  
  # Fit expanded model with bootstrap
  new_formula <- update(base_formula, paste("~ . +", predictor_var))
  
  # Create save prefix
  save_prefix <- NULL
  if (save_matrices && !is.null(matrix_dir)) {
    save_prefix <- paste0("iter", iteration, "_", gsub("[^A-Za-z0-9]", "", predictor_name))
  }
  
  # Run bootstrap for expanded model
  expanded_boot <- run_bootstrap_model(
    formula = new_formula,
    data = test_data,
    tree = test_tree,
    n_boot = n_bootstrap,
    save_prefix = if(!is.null(save_prefix)) paste0(save_prefix, "_expanded") else NULL,
    save_matrices = save_matrices,
    matrix_dir = matrix_dir
  )
  
  expanded_aic <- -2 * expanded_boot$fit$logLik + 2 * expanded_boot$fit$d
  aic_improvement <- base_subset_aic - expanded_aic
  
  message("   Model with ", predictor_name, " AIC: ", round(expanded_aic, 2))
  message("   AIC improvement: ", round(aic_improvement, 2))
  
  return(list(
    predictor = predictor_name,
    variable = predictor_var,
    n_species = nrow(test_data),
    base_aic = base_subset_aic,
    expanded_aic = expanded_aic,
    aic_improvement = aic_improvement,
    expanded_fit = expanded_boot$fit,
    expanded_bootstrap = expanded_boot,
    data = test_data,
    tree = test_tree,
    formula = new_formula
  ))
}

# Function to process model coefficients into comprehensive table format
process_model_coefficients <- function(
  model_fit, bootstrap_results, analysis_name, iteration,
  predictor_added, predictor_variable, base_formula,
  previous_aic, expanded_aic, aic_improvement,
  model_selected, n_species
) {
  
  # Get coefficient summary
  coef_summary <- summary(model_fit)$coefficients
  
  # Initialize results list
  results_list <- list()
  
  for (i in 1:nrow(coef_summary)) {
    param_name <- rownames(coef_summary)[i]
    
    # Get bootstrap info for this parameter
    boot_row <- bootstrap_results$coefficients[bootstrap_results$coefficients$Parameter == param_name, ]
    
    if (nrow(boot_row) == 0) {
      # If no bootstrap info (shouldn't happen), use fallback
      boot_mean <- coef_summary[i, "Estimate"]
      boot_sd <- coef_summary[i, "StdErr"]
      ci_lower <- boot_mean - 1.96 * boot_sd
      ci_upper <- boot_mean + 1.96 * boot_sd
    } else {
      boot_mean <- boot_row$Boot_Mean[1]
      boot_sd <- boot_row$Boot_SD[1]
      ci_lower <- boot_row$CI_Lower[1]
      ci_upper <- boot_row$CI_Upper[1]
    }
    
    # Calculate odds ratios
    or <- exp(coef_summary[i, "Estimate"])
    or_lower <- exp(ci_lower)
    or_upper <- exp(ci_upper)
    
    results_list[[i]] <- data.frame(
      Analysis_Name = analysis_name,
      Iteration = iteration,
      Base_Model_Formula = base_formula,
      Predictor_Added = predictor_added,
      Predictor_Variable = predictor_variable,
      Previous_Model_AIC = previous_aic,
      Expanded_Model_AIC = expanded_aic,
      AIC_Improvement = aic_improvement,
      Model_Selected = model_selected,
      N_Species_In_Model = n_species,
      Coefficient_Name = param_name,
      Coefficient_Estimate = coef_summary[i, "Estimate"],
      Coefficient_SE = coef_summary[i, "StdErr"],
      Coefficient_pvalue = coef_summary[i, "p.value"],
      Coefficient_Bootstrap_Mean_Estimate = boot_mean,
      Coefficient_Bootstrap_SD = boot_sd,
      Coefficient_CI_Lower_2.5 = ci_lower,
      Coefficient_CI_Upper_97.5 = ci_upper,
      Parameter_Odds_Ratio = or,
      Parameter_Odds_Ratio_CI_Lower_2.5 = or_lower,
      Parameter_Odds_Ratio_CI_Upper_97.5 = or_upper,
      OR_pvalue = coef_summary[i, "p.value"],
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, results_list))
}

# Function to create final comparison table
create_final_comparison <- function(single_param_results, iteration_paths, comprehensive_results) {
  
  final_models <- data.frame()
  
  for (analysis_name in c("CB_vs_FS_TerrWS_Mass", "FS_vs_CB_TerrWS_Mass")) {
    
    # Get original base info
    base_info <- single_param_results[[analysis_name]]
    if (is.null(base_info)) next
    
    # Get final iteration info
    analysis_paths <- iteration_paths[iteration_paths$Analysis_Name == analysis_name,]
    
    if (nrow(analysis_paths) == 0) {
      # No improvements were made
      final_models <- rbind(final_models, data.frame(
        Analysis_Name = analysis_name,
        Original_Base_AIC = base_info$original_base_aic,
        Final_Model_AIC = base_info$original_base_aic,
        Total_AIC_Improvement = 0,
        N_Predictors_Added = 0,
        Predictors_Added = "None",
        Original_MainPred_Effect = get_mainpred_effect(base_info$base_model, analysis_name),
        Final_MainPred_Effect = get_mainpred_effect(base_info$base_model, analysis_name),
        Total_MainPred_Change = 0,
        stringsAsFactors = FALSE
      ))
    } else {
      # Get final model info
      final_iter <- analysis_paths[nrow(analysis_paths),]
      
      # Get MainPred effects
      orig_mainpred <- get_mainpred_from_comprehensive(
        comprehensive_results, analysis_name, iteration = 0
      )
      final_mainpred <- get_mainpred_from_comprehensive(
        comprehensive_results, analysis_name, iteration = max(analysis_paths$Iteration)
      )
      
      final_models <- rbind(final_models, data.frame(
        Analysis_Name = analysis_name,
        Original_Base_AIC = base_info$original_base_aic,
        Final_Model_AIC = final_iter$Current_AIC,
        Total_AIC_Improvement = final_iter$Cumulative_AIC_Improvement,
        N_Predictors_Added = nrow(analysis_paths),
        Predictors_Added = paste(analysis_paths$Predictor_Added, collapse = " + "),
        Original_MainPred_Effect = orig_mainpred,
        Final_MainPred_Effect = final_mainpred,
        Total_MainPred_Change = final_mainpred - orig_mainpred,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(final_models)
}

# Helper function to get MainPred effect from model
get_mainpred_effect <- function(model, analysis_name) {
  if (grepl("^FS_vs_CB", analysis_name)) {
    main_pred <- "HighConfidence_Coop"
  } else {
    main_pred <- "FemaleSong_Agg01"
  }
  
  coefs <- coef(model)
  if (main_pred %in% names(coefs)) {
    return(coefs[main_pred])
  }
  return(NA)
}

# Helper function to get MainPred from comprehensive results
get_mainpred_from_comprehensive <- function(comp_results, analysis_name, iteration) {
  if (grepl("^FS_vs_CB", analysis_name)) {
    main_pred <- "HighConfidence_Coop"
  } else {
    main_pred <- "FemaleSong_Agg01"
  }
  
  # For iteration 0, need to look at the base model
  if (iteration == 0) {
    # This would need to be extracted from the original base model
    # For now, return NA
    return(NA)
  }
  
  pred_row <- comp_results %>%
    filter(Analysis_Name == analysis_name,
           Iteration == iteration,
           Coefficient_Name == main_pred,
           Model_Selected == TRUE) %>%
    slice(1)
  
  if (nrow(pred_row) > 0) {
    return(pred_row$Coefficient_Estimate[1])
  }
  return(NA)
}

# Load bootstrap model function from the single-parameter script
# This should be available from sourcing stepwise_expansion_bootstrap_updated.R
# If not, include it here:

run_bootstrap_model <- function(formula, data, tree, n_boot = 1000, 
                               method = "logistic_MPLE", save_prefix = NULL,
                               save_matrices = TRUE, matrix_dir = NULL) {
  
  message("      Running ", n_boot, " bootstrap iterations...")
  
  # Fit original model
  original_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4
  )
  
  # Get bootstrap results
  message("      Fitting model with bootstrap=", n_boot, "...")
  boot_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4,
    boot = n_boot
  )
  
  # Extract bootstrap results
  if (!is.null(boot_fit$bootstrap)) {
    boot_matrix <- boot_fit$bootstrap
    
    # Save matrix if requested
    matrix_file <- NULL
    if (save_matrices && !is.null(save_prefix) && !is.null(matrix_dir)) {
      matrix_file <- file.path(matrix_dir, paste0(save_prefix, "_boot", n_boot, "_matrix.rds"))
      saveRDS(boot_matrix, matrix_file, compress = TRUE)
      message("      Saved bootstrap matrix to: ", basename(matrix_file))
    }
    
    # Calculate statistics
    original_coef <- coef(original_fit)
    coef_names <- names(original_coef)
    n_coef <- length(coef_names)
    
    # Handle potential mismatch with bootstrap matrix
    if (ncol(boot_matrix) > n_coef) {
      boot_matrix_coef <- boot_matrix[, 1:n_coef, drop = FALSE]
    } else {
      boot_matrix_coef <- boot_matrix
    }
    
    boot_means <- colMeans(boot_matrix_coef, na.rm = TRUE)
    boot_sds <- apply(boot_matrix_coef, 2, sd, na.rm = TRUE)
    boot_lower <- apply(boot_matrix_coef, 2, quantile, probs = 0.025, na.rm = TRUE)
    boot_upper <- apply(boot_matrix_coef, 2, quantile, probs = 0.975, na.rm = TRUE)
    
    coef_summary <- data.frame(
      Parameter = coef_names,
      Estimate = original_coef,
      Boot_Mean = boot_means,
      Boot_SD = boot_sds,
      CI_Lower = boot_lower,
      CI_Upper = boot_upper,
      stringsAsFactors = FALSE
    )
    
    n_successful <- sum(complete.cases(boot_matrix_coef))
    n_converged <- sum(!is.na(boot_matrix_coef[,1]))
    convergence_rate <- n_converged / n_boot
    
  } else {
    # Fallback
    message("      Warning: Bootstrap results not available, using standard errors")
    coef_summary <- summary(original_fit)$coefficients
    coef_summary <- data.frame(
      Parameter = rownames(coef_summary),
      Estimate = coef_summary[, "Estimate"],
      Boot_Mean = coef_summary[, "Estimate"],
      Boot_SD = coef_summary[, "StdErr"],
      CI_Lower = coef_summary[, "Estimate"] - 1.96 * coef_summary[, "StdErr"],
      CI_Upper = coef_summary[, "Estimate"] + 1.96 * coef_summary[, "StdErr"],
      stringsAsFactors = FALSE
    )
    matrix_file <- NULL
    n_successful <- 1
    n_converged <- 1
    convergence_rate <- 1
  }
  
  return(list(
    fit = original_fit,
    coefficients = coef_summary,
    bootstrap_matrix_file = matrix_file,
    n_successful_boots = n_successful,
    n_converged = n_converged,
    convergence_rate = convergence_rate
  ))
}

# Usage example:
results <- run_iterative_stepwise_bootstrap(
  single_param_results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  n_bootstrap = 500,
  aic_threshold = 2
)
