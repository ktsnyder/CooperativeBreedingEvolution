# Phylogenetic GLM Stepwise Expansion Script
# Combines functionality of stepwise_expansion_bootstrap_updated.R and iterative_stepwise_bootstrap_KTS.R
# Then split from "phyloglm_unified (2).R"
# Supports both RDS input (from previous analyses) and direct formula input
# Performs fast iterations with bootstrap only for initial and final models

library(phylolm)
library(ape)

# Define available predictors for stepwise expansion
define_available_predictors <- function() {
  list(
    Migration = list(
      name = "Migration (1-3)",
      var = "Migration_num", 
      values_func = function(data) data$Migration_num,
      skip_conditions = "Migration_num"
    ),
    Latitude = list(
      name = "Absolute Latitude",
      var = "abs_Latitude_normalized",
      values_func = function(data) data$abs_Latitude_normalized,
      skip_conditions = "abs_Latitude_normalized|abs_Latitude|absLat|Latitude"
    ),
    Wing = list(
      name = "Wing Dimorphism",
      var = "PercentAbsLogWingDimorphism_normalized",
      values_func = function(data) data$PercentAbsLogWingDimorphism_normalized,
      skip_conditions = "PercentAbsLogWingDimorphism_normalized|WingDimorphism|PercentAbsLogWingDimorphism"
    ),
    Plumage = list(
      name = "Plumage Dimorphism",
      var = "logMaleFemalePlumageDiffAbs_normalized",
      values_func = function(data) data$logMaleFemalePlumageDiffAbs_normalized,
      skip_conditions = "logMaleFemalePlumageDiffAbs_normalized|PlumageDiff|logMaleFemalePlumageDiff"
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
    ),
    Coop = list(
      name = "Cooperative Breeding",
      var = "HighConfidence_Coop",
      values_func = function(data) as.factor(data$HighConfidence_Coop),
      skip_conditions = "HighConfidence_Coop|Coop"
    ),
    Mass = list(
      name = "Body Mass",
      var = "logMass_normalized",
      values_func = function(data) data$logMass_normalized,
      skip_conditions = "logMass_normalized|logMass|Mass_AVONET|Mass"
    ),
    Polygyny = list(
      name = "Polygyny",
      var = "Final.polygyny",
      values_func = function(data) data$Final.polygyny,
      skip_conditions = "Polygyny"
    )
  )
}

# Get remaining predictors not yet in current formula
get_remaining_predictors <- function(current_formula, available_predictors, full_data) {
  formula_str <- paste(deparse(current_formula), collapse = " ")
  remaining <- list()
  
  for (pred_name in names(available_predictors)) {
    pred <- available_predictors[[pred_name]]
    if (!grepl(pred$skip_conditions, formula_str)) {
      remaining[[pred_name]] <- pred
    }
  }
  
  return(remaining)
}

# Prepare data for a specific formula
prepare_data_for_formula <- function(formula, full_data, tree) {
  all_vars <- all.vars(formula)
  data_subset <- full_data[, c("species", all_vars)]
  data_clean <- na.omit(data_subset)
  tree_subset <- keep.tip(tree, data_clean$species)
  data_final <- data_clean[match(tree_subset$tip.label, data_clean$species), ]
  
  return(list(data = data_final, tree = tree_subset))
}

# Test predictors without bootstrap (for speed during selection)
test_predictors_fast <- function(base_formula, base_data, base_tree, base_aic,
                                remaining_predictors, full_data, verbose = TRUE) {
  
  results <- list()
  
  for (pred_name in names(remaining_predictors)) {
    pred <- remaining_predictors[[pred_name]]
    
    if (verbose) cat("  Testing:", pred$name, "\n")
    
    # Add predictor values to data
    test_data <- base_data
    test_data[[pred$var]] <- pred$values_func(full_data)[match(test_data$species, full_data$species)]
    
    # Handle categorical variables
    if (pred$var == "GeographicRegion_Jetz") {
      test_data[[pred$var]] <- as.factor(test_data[[pred$var]])
    }
    
    # Remove NAs
    test_data <- test_data[!is.na(test_data[[pred$var]]), ]
    
    if (nrow(test_data) < 100) {
      if (verbose) cat("    Not enough species with data (", nrow(test_data), ")\n")
      next
    }
    
    # Match tree
    test_tree <- keep.tip(base_tree, test_data$species)
    test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
    
    # Fit base model to subset
    tryCatch({
      base_subset_fit <- phyloglm(
        formula = base_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
    }, error = function(e) {
      if (verbose) cat("    Base model failed to converge\n")
      return(NULL)
    })
    
    # Fit expanded model
    new_formula <- update(base_formula, paste("~ . +", pred$var))
    
    tryCatch({
      expanded_fit <- phyloglm(
        formula = new_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      expanded_aic <- -2 * expanded_fit$logLik + 2 * expanded_fit$d
      aic_improvement <- base_subset_aic - expanded_aic
      
      if (verbose) cat("    AIC improvement:", round(aic_improvement, 2), "\n")
      
      results[[pred_name]] <- list(
        predictor_name = pred$name,
        predictor_var = pred$var,
        formula = new_formula,
        data = test_data,
        tree = test_tree,
        model = expanded_fit,
        aic_improvement = aic_improvement,
        n_species = nrow(test_data)
      )
    }, error = function(e) {
      if (verbose) cat("    Expanded model failed to converge\n")
    })
  }
  
  return(results)
}

# Find best improvement from test results
find_best_improvement <- function(iteration_results, aic_threshold) {
  best_improvement <- NULL
  best_aic_improvement <- 0
  
  for (result in iteration_results) {
    if (!is.null(result$aic_improvement) &&
        result$aic_improvement >= aic_threshold && 
        result$aic_improvement > best_aic_improvement) {
      best_improvement <- result
      best_aic_improvement <- result$aic_improvement
    }
  }
  
  return(best_improvement)
}

# Main stepwise expansion function
run_stepwise_expansion <- function(
  formulas,                    # List of starting formulas
  data_path,                   # Path to data file
  tree_path,                   # Path to tree file
  aic_threshold = 2,           # AIC improvement threshold
  n_bootstrap = 100,           # Bootstrap iterations for final models
  max_iterations = 10,         # Safety limit
  output_dir = "Outputs/PhyloGLM_outputs/Stepwise",
  verbose = TRUE
) {
  
  # Create output directory
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  main_output_dir <- file.path(output_dir, paste0("stepwise_", timestamp))
  dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data and tree
  if (verbose) cat("Loading data and tree...\n")
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Initialize results storage
  all_results <- list()
  summary_table <- data.frame()
  
  # Define available predictors
  available_predictors <- define_available_predictors()
  
  # Process each starting formula
  for (i in seq_along(formulas)) {
    formula_input <- formulas[[i]]
    model_name <- paste0("model_", i)
    
    if (verbose) {
      cat("\n", paste(rep("=", 50), collapse=""), "\n")
      cat("Processing", model_name, ":", deparse(formula_input), "\n")
      cat(paste(rep("=", 50), collapse=""), "\n\n")
    }
    
    # Prepare initial data
    prepared <- prepare_data_for_formula(formula_input, full_data, tree)
    rownames(prepared$data) <- prepared$data$species
    
    # Fit initial model
    initial_fit <- phyloglm(
      formula = formula_input,
      data = prepared$data,
      phy = prepared$tree,
      method = "logistic_MPLE",
      btol = 50,
      log.alpha.bound = 4
    )
    
    initial_aic <- -2 * initial_fit$logLik + 2 * initial_fit$d
    
    # Initialize tracking for this model
    current_formula <- formula_input
    current_data <- prepared$data
    current_tree <- prepared$tree
    current_aic <- initial_aic
    iteration_history <- data.frame()
    
    # Iterative expansion
    for (iteration in 1:max_iterations) {
      if (verbose) cat("\n--- Iteration", iteration, "---\n")
      
      # Get remaining predictors
      remaining_predictors <- get_remaining_predictors(current_formula, available_predictors, full_data)
      
      if (length(remaining_predictors) == 0) {
        if (verbose) cat("No more predictors to test\n")
        break
      }
      
      # Test each remaining predictor
      iteration_results <- test_predictors_fast(
        base_formula = current_formula,
        base_data = current_data,
        base_tree = current_tree,
        base_aic = current_aic,
        remaining_predictors = remaining_predictors,
        full_data = full_data,
        verbose = verbose
      )
      
      # Find best improvement
      best_improvement <- find_best_improvement(iteration_results, aic_threshold)
      
      if (is.null(best_improvement)) {
        if (verbose) cat("No predictor improves model by threshold\n")
        break
      }
      
      # Update model with best predictor
      if (verbose) {
        cat("\nBest predictor:", best_improvement$predictor_name, "\n")
        cat("AIC improvement:", round(best_improvement$aic_improvement, 2), "\n")
      }
      
      # Record this iteration
      iteration_history <- rbind(iteration_history, data.frame(
        Model = model_name,
        Iteration = iteration,
        Predictor_Added = best_improvement$predictor_name,
        AIC_Improvement = best_improvement$aic_improvement,
        N_Species = best_improvement$n_species,
        stringsAsFactors = FALSE
      ))
      
      # Update current model state
      current_formula <- best_improvement$formula
      current_data <- best_improvement$data
      current_tree <- best_improvement$tree
      current_aic <- -2 * best_improvement$model$logLik + 2 * best_improvement$model$d
    }
    
    # Run final model with bootstrap
    if (verbose) cat("\nRunning final model with bootstrap...\n")
    
    final_bootstrap <- run_bootstrap_model(
      formula = current_formula,
      data = current_data,
      tree = current_tree,
      n_boot = n_bootstrap,
      save_prefix = paste0(model_name, "_final"),
      save_matrices = FALSE,
      matrix_dir = main_output_dir,
      save_coefficient_csv = TRUE,
      use_bootstrap_pvalues = FALSE
    )
    
    # Store results
    all_results[[model_name]] <- list(
      initial_formula = formula_input,
      initial_aic = initial_aic,
      final_formula = current_formula,
      final_aic = current_aic,
      total_improvement = initial_aic - current_aic,
      n_iterations = nrow(iteration_history),
      iteration_history = iteration_history,
      final_bootstrap = final_bootstrap
    )
    
    # Add to summary table
    summary_table <- rbind(summary_table, data.frame(
      Model = model_name,
      Initial_Formula = paste(deparse(formula_input), collapse = " "),
      Final_Formula = paste(deparse(current_formula), collapse = " "),
      Initial_AIC = initial_aic,
      Final_AIC = current_aic,
      Total_AIC_Improvement = initial_aic - current_aic,
      N_Iterations = nrow(iteration_history),
      stringsAsFactors = FALSE
    ))
  }
  
  # Save results
  write.csv(summary_table, 
            file.path(main_output_dir, paste0("stepwise_summary_boot", n_bootstrap, ".csv")),
            row.names = FALSE)
  
  # Save detailed iteration history
  all_iterations <- data.frame()
  for (model_name in names(all_results)) {
    if (nrow(all_results[[model_name]]$iteration_history) > 0) {
      all_iterations <- rbind(all_iterations, all_results[[model_name]]$iteration_history)
    }
  }
  
  if (nrow(all_iterations) > 0) {
    write.csv(all_iterations,
              file.path(main_output_dir, paste0("iteration_history_boot", n_bootstrap, ".csv")),
              row.names = FALSE)
  }
  
  # Save complete results
  saveRDS(all_results, file.path(main_output_dir, paste0("stepwise_results_boot", n_bootstrap, ".rds")))
  
  if (verbose) {
    cat("\n", paste(rep("=", 50), collapse=""), "\n")
    cat("Stepwise expansion complete\n")
    cat("Results saved to:", main_output_dir, "\n")
    cat(paste(rep("=", 50), collapse=""), "\n")
  }
  
  return(all_results)
}