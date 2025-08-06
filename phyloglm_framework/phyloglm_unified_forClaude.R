# Unified Phylogenetic GLM Stepwise Expansion Script
# Combines functionality of stepwise_expansion_bootstrap_updated.R and iterative_stepwise_bootstrap_KTS.R
# Then split from "phyloglm_unified (2).R"
# Supports both RDS input (from previous analyses) and direct formula input
# Performs fast iterations with bootstrap only for initial and final models
# KTS edited 7/29/2025
#
# Key features:
# - Accepts either RDS file with previous results or list of formulas
# - Expands top N models independently (configurable)
# - Fast AIC-based selection during iterations (no bootstrap)
# - Full bootstrap only for initial and final models
# - Options for bootstrap-based or parametric p-values
# - Saves both parametric and bootstrap p-values for comparison
# - Detects and warns about bootstrap bias
# - Ensures proper parameter alignment across all calculations
# - Comprehensive output with expansion history and convergence analysis
# - Can replace logMass_AVONET with logMass_normalized throughout (use_normalized_mass parameter)

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

# Source helper functions (ensure this path is correct)
# source("phyloglm_framework/batch_runner_helpers.R")

# Main unified function for stepwise expansion
run_unified_stepwise_expansion <- function(
  input_source = NULL,           # Either RDS file path or list of formulas
  data_path = NULL,              # Required if input_source is formulas
  tree_path = NULL,              # Required if input_source is formulas
  top_n_models = 1,              # Number of top models to expand (ignored if formulas provided)
  output_base_dir = "Outputs/PhyloglmResults/",
  aic_threshold = 2,
  n_bootstrap = 500,
  max_iterations = 10,           # Safety limit
  save_boot_matrices = TRUE,
  save_coefficient_csv = TRUE,   # Whether to save coefficient CSVs
  use_bootstrap_pvalues = FALSE, # Whether to use bootstrap-based p-values
  bias_threshold_sd = 1.0,       # Threshold for bootstrap bias warning
  verbose = TRUE,
  detailed_logging = TRUE,       # Whether to save detailed expansion logs
  use_normalized_mass = TRUE,    # Whether to use logMass_normalized instead of logMass_AVONET
  ncores = 1                     # For potential parallelization
) {
  
  # Create output directory
  timestamp <- format(Sys.Date(), "%Y%m%d")
  main_output_dir <- file.path(output_base_dir, paste0("stepwise_unified_", timestamp))
  dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize results storage
  all_expansion_paths <- list()
  
  # Step 1: Handle input and prepare base models
  if (verbose) message("\n=== UNIFIED STEPWISE EXPANSION ===\n")
  
  base_models <- prepare_base_models(
    input_source = input_source,
    data_path = data_path,
    tree_path = tree_path,
    top_n_models = top_n_models,
    n_bootstrap = n_bootstrap,
    main_output_dir = main_output_dir,
    save_boot_matrices = save_boot_matrices,
    use_bootstrap_pvalues = use_bootstrap_pvalues,
    use_normalized_mass = use_normalized_mass,
    verbose = verbose
  )
  
  if (length(base_models) == 0) {
    stop("No valid base models found")
  }
  
  # Step 2: Define available predictors
  available_predictors <- define_available_predictors()
  
  # Step 3: Perform stepwise expansion for each base model
  model_counter <- 1
  
  for (model_name in names(base_models)) {
    if (verbose) {
      message("\n", paste(rep("=", 60), collapse=""))
      message("Processing model ", model_counter, " of ", length(base_models), ": ", model_name)
      message(paste(rep("=", 60), collapse=""))
    }
    
    # Initialize expansion path for this model
    expansion_path <- list(
      model_name = model_name,
      analysis_type = base_models[[model_name]]$analysis_type,
      initial_formula = base_models[[model_name]]$formula,
      initial_aic = base_models[[model_name]]$aic,
      initial_n_species = nrow(base_models[[model_name]]$data),
      initial_bootstrap = base_models[[model_name]]$bootstrap_results,
      iterations = list(),
      final_formula = NULL,
      final_aic = NULL,
      final_bootstrap = NULL,
      status = "active"
    )
    
    # Current model state
    current_formula <- expansion_path$initial_formula
    current_data <- base_models[[model_name]]$data
    current_tree <- base_models[[model_name]]$tree
    current_aic <- expansion_path$initial_aic
    
    # Perform iterative expansion
    for (iteration in 1:max_iterations) {
      if (verbose) message("\n--- Iteration ", iteration, " ---")
      
      # Get remaining predictors
      remaining_predictors <- get_remaining_predictors(
        current_formula = current_formula,
        available_predictors = available_predictors,
        full_data = base_models[[model_name]]$full_data
      )
      
      if (length(remaining_predictors) == 0) {
        if (verbose) message("No more predictors to test")
        expansion_path$status <- "completed_no_predictors"
        break
      }
      
      # Test each remaining predictor
      iteration_results <- test_predictors_fast(
        base_formula = current_formula,
        base_data = current_data,
        base_tree = current_tree,
        base_aic = current_aic,
        remaining_predictors = remaining_predictors,
        full_data = base_models[[model_name]]$full_data,
        verbose = verbose,
        use_normalized_mass = use_normalized_mass
      )
      
      # Find best improvement
      best_improvement <- find_best_improvement(
        iteration_results = iteration_results,
        aic_threshold = aic_threshold
      )
      
      # Always store iteration results, even if no predictor is selected
      expansion_path$iterations[[iteration]] <- list(
        predictors_tested = names(iteration_results),
        selected_predictor = if (!is.null(best_improvement)) best_improvement$predictor_name else NA,
        aic_improvement = if (!is.null(best_improvement)) best_improvement$aic_improvement else NA,
        new_aic = if (!is.null(best_improvement)) best_improvement$new_aic else current_aic,
        n_species = if (!is.null(best_improvement)) best_improvement$n_species_final else nrow(current_data),
        # Store all test results if detailed logging is enabled
        all_test_results = if (detailed_logging) iteration_results else NULL,
        # Summary of all tests
        test_summary = lapply(iteration_results, function(res) {
          list(
            predictor = res$predictor_name,
            status = res$status,
            n_species_final = res$n_species_final,
            aic_improvement = if (!is.null(res$aic_improvement)) res$aic_improvement else NA,
            formula = if (!is.null(res$expanded_formula)) 
              paste(deparse(res$expanded_formula), collapse = " ") else NA
          )
        })
      )
      
      if (is.null(best_improvement)) {
        if (verbose) message("No predictor improves model by threshold")
        expansion_path$status <- "completed_no_improvement"
        break
      }
      
      # Update model with best predictor
      if (verbose) {
        message("\nBest predictor: ", best_improvement$predictor_name)
        message("AIC improvement: ", round(best_improvement$aic_improvement, 2))
      }
      
      # Update current model state
      current_formula <- best_improvement$formula
      current_data <- best_improvement$data
      current_tree <- best_improvement$tree
      current_aic <- best_improvement$new_aic
    }
    
    # Store final state
    expansion_path$final_formula <- current_formula
    expansion_path$final_aic <- current_aic
    expansion_path$n_iterations <- length(expansion_path$iterations)
    
    # Set final status if not already set
    if (expansion_path$status == "active") {
      if (length(expansion_path$iterations) > 0) {
        expansion_path$status <- "completed_with_expansions"
      } else {
        expansion_path$status <- "completed_no_expansions"
      }
    }
    
    # Run final model with bootstrap
    if (verbose) message("\nRunning final model with bootstrap...")
    
    final_bootstrap <- run_bootstrap_model(
      formula = current_formula,
      data = current_data,
      tree = current_tree,
      n_boot = n_bootstrap,
      save_prefix = paste0(model_name, "_final"),
      save_matrices = save_boot_matrices,
      matrix_dir = file.path(main_output_dir, 
                           paste0(expansion_path$analysis_type, "_", 
                                 round(expansion_path$initial_aic, 1)),
                           "bootstrap_matrices"),
      save_coefficient_csv = save_coefficient_csv,
      use_bootstrap_pvalues = use_bootstrap_pvalues,
      bias_threshold_sd = bias_threshold_sd
    )
    
    expansion_path$final_bootstrap <- final_bootstrap
    expansion_path$final_n_species <- nrow(current_data)
    
    # Save expansion path
    all_expansion_paths[[model_name]] <- expansion_path
    
    # Create model-specific output directory using INITIAL AIC
    model_output_dir <- file.path(main_output_dir, 
                                paste0(expansion_path$analysis_type, "_", 
                                      round(expansion_path$initial_aic, 1)))
    dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save detailed results for this model
    save_model_results(
      expansion_path = expansion_path,
      output_dir = model_output_dir,
      n_bootstrap = n_bootstrap
    )
    
    model_counter <- model_counter + 1
  }
  
  # Create summary across all paths
  create_summary_results(
    all_expansion_paths = all_expansion_paths,
    output_dir = main_output_dir,
    n_bootstrap = n_bootstrap
  )
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse=""))
    message("UNIFIED STEPWISE EXPANSION COMPLETE")
    message("Results saved to: ", main_output_dir)
    message(paste(rep("=", 60), collapse=""))
  }
  
  return(list(
    expansion_paths = all_expansion_paths,
    output_dir = main_output_dir
  ))
}

# Function to replace mass variable in formula if needed
replace_mass_in_formula <- function(formula, use_normalized_mass = TRUE) {
  if (!use_normalized_mass) {
    return(formula)
  }
  
  # Convert formula to character
  formula_str <- paste(deparse(formula), collapse = " ")
  
  # Replace various mass variable names with logMass_normalized
  formula_str <- gsub("logMass_AVONET", "logMass_normalized", formula_str)
  formula_str <- gsub("\\blogMass\\b", "logMass_normalized", formula_str)
  formula_str <- gsub("Mass_AVONET", "logMass_normalized", formula_str)
  
  # Convert back to formula
  return(as.formula(formula_str))
}

# Function to ensure data has the correct mass column
ensure_mass_column <- function(data, use_normalized_mass = TRUE) {
  if (!use_normalized_mass) {
    return(data)
  }
  
  # Check if we need to rename the column
  if ("logMass_AVONET" %in% names(data) && !"logMass_normalized" %in% names(data)) {
    data$logMass_normalized <- data$logMass_AVONET
    # Optionally remove the old column to avoid confusion
    # data$logMass_AVONET <- NULL
  } else if ("Mass_AVONET" %in% names(data) && !"logMass_normalized" %in% names(data)) {
    # If it's not already log-transformed, we might need to log it
    # But for now, just rename
    data$logMass_normalized <- data$Mass_AVONET
  }
  
  return(data)
}

# Function to prepare base models from input
prepare_base_models <- function(input_source, data_path, tree_path, top_n_models, 
                               n_bootstrap, main_output_dir, save_boot_matrices, use_bootstrap_pvalues = TRUE, 
                               bias_threshold_sd = 1.0, save_coefficient_csv = TRUE, use_normalized_mass = TRUE, verbose) {
  
  base_models <- list()
  
  if (is.character(input_source) && file.exists(input_source)) {
    # Input is an RDS file with previous results
    if (verbose) message("Loading results from: ", input_source)
    
    all_results <- readRDS(input_source)
    
    # Extract top n models for each analysis type
    analysis_types <- c("CB_vs_FS_Terr3_Mass", "FS_vs_CB_Terr3_Mass", 
                       "CB_vs_FS_TerrWS_Mass", "FS_vs_CB_TerrWS_Mass")
    
    for (analysis_type in analysis_types) {
      if (!is.null(all_results[[analysis_type]])) {
        result <- all_results[[analysis_type]]
        
        # Get top n models
        if (!is.null(result$comparison$comparison)) {
          top_models <- head(result$comparison$comparison, top_n_models)
          
          for (i in 1:nrow(top_models)) {
            model_name <- paste0(analysis_type, "_model", i)
            model_info <- result$models$models[[top_models$Model[i]]]
            
            # Run bootstrap for this base model if not already done
            if (verbose) message("Preparing base model: ", model_name)
            
            # Get full data (may need to load separately if not in results)
            if (!is.null(result$prepared_data$full_data)) {
              full_data_for_model <- result$prepared_data$full_data
            } else {
              # If full data not in results, try to load it
              if (!is.null(data_path) && file.exists(data_path)) {
                full_data_for_model <- read.csv(data_path)
              } else {
                # Try default location based on other scripts
                default_data_path <- "Data_R_2025-06-09.csv"
                if (file.exists(default_data_path)) {
                  if (verbose) message("  Loading full data from default location: ", default_data_path)
                  full_data_for_model <- read.csv(default_data_path)
                } else {
                  stop("Full data not found in results and data_path not provided or not found")
                }
              }
            }
            
            # Replace mass variable in formula if needed
            model_formula <- replace_mass_in_formula(formula(model_info), use_normalized_mass)
            
            # Ensure data has the correct mass column
            prepared_data <- ensure_mass_column(result$prepared_data$data, use_normalized_mass)
            full_data_for_model <- ensure_mass_column(full_data_for_model, use_normalized_mass)
            
            base_models[[model_name]] <- list(
              analysis_type = analysis_type,
              formula = model_formula,
              model = model_info,
              aic = top_models$AIC[i],
              data = prepared_data,
              tree = result$prepared_data$tree,
              full_data = full_data_for_model,
              bootstrap_results = NULL  # Will be filled below
            )
            
            # Run bootstrap for this base model
            if (verbose) message("  Running bootstrap for base model...")
            
            boot_results <- run_bootstrap_model(
              formula = model_formula,
              data = prepared_data,
              tree = result$prepared_data$tree,
              n_boot = n_bootstrap,
              save_prefix = paste0(model_name, "_initial"),
              save_matrices = save_boot_matrices,
              matrix_dir = file.path(main_output_dir, "initial_bootstraps"),
              save_coefficient_csv = save_coefficient_csv,
              use_bootstrap_pvalues = use_bootstrap_pvalues,
              bias_threshold_sd = bias_threshold_sd
            )
            
            base_models[[model_name]]$bootstrap_results <- boot_results
          }
        }
      }
    }
    
  } else if (is.list(input_source)) {
    # Input is a list of formulas
    if (is.null(data_path) || is.null(tree_path)) {
      stop("data_path and tree_path must be provided when using formula input")
    }
    
    if (verbose) message("Loading data and tree...")
    full_data <- read.csv(data_path)
    # Ensure full data has the correct mass column
    full_data <- ensure_mass_column(full_data, use_normalized_mass)
    tree <- read.nexus(tree_path)
    
    formula_counter <- 1
    for (formula_input in input_source) {
      if (verbose) message("\nProcessing formula ", formula_counter, "...")
      
      # Replace mass variable in formula if needed
      formula_input <- replace_mass_in_formula(formula_input, use_normalized_mass)
      
      # Infer analysis type
      analysis_type <- infer_analysis_type(formula_input)
      model_name <- paste0(analysis_type, "_custom", formula_counter)
      
      # Prepare data for this formula
      prepared <- prepare_data_for_formula(
        formula = formula_input,
        full_data = full_data,
        tree = tree
      )
      
      rownames(prepared$data) <- prepared$data$species
      
      # Run initial model with bootstrap
      if (verbose) message("Running initial model with bootstrap...")
      
      initial_fit <- phyloglm(
        formula = formula_input,
        data = prepared$data,
        phy = prepared$tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      
      initial_aic <- -2 * initial_fit$logLik + 2 * initial_fit$d
      
      # Run bootstrap
      boot_results <- run_bootstrap_model(
        formula = formula_input,
        data = prepared$data,
        tree = prepared$tree,
        n_boot = n_bootstrap,
        save_prefix = paste0(model_name, "_initial"),
        save_matrices = save_boot_matrices,
        matrix_dir = file.path(main_output_dir, "initial_bootstraps"),
        save_coefficient_csv = save_coefficient_csv,
        use_bootstrap_pvalues = use_bootstrap_pvalues,
        bias_threshold_sd = bias_threshold_sd
      )
      
      base_models[[model_name]] <- list(
        analysis_type = analysis_type,
        formula = formula_input,
        model = initial_fit,
        aic = initial_aic,
        data = prepared$data,
        tree = prepared$tree,
        full_data = full_data,
        bootstrap_results = boot_results
      )
      
      formula_counter <- formula_counter + 1
    }
  } else {
    stop("input_source must be either an RDS file path or a list of formulas")
  }
  
  return(base_models)
}

# Function to infer analysis type from formula
infer_analysis_type <- function(formula) {
  formula_str <- paste(deparse(formula), collapse = " ")
  
  # Determine response variable
  if (grepl("^HighConfidence_Coop", formula_str)) {
    response_type <- "CB_vs_FS"
  } else if (grepl("^FemaleSong_Agg01", formula_str)) {
    response_type <- "FS_vs_CB"
  } else {
    response_type <- "Unknown"
  }
  
  # Determine territory type
  if (grepl("Territory_12vs3|Territory_num", formula_str)) {
    terr_type <- "Terr3"
  } else if (grepl("TerritorialityWeakVsStrong", formula_str)) {
    terr_type <- "TerrWS"
  } else {
    terr_type = "NoTerr"
  }
  
  # Check for mass
  mass_type <- ifelse(grepl("logMass_normalized|logMass|Mass_AVONET", formula_str), "Mass", "NoMass")
  
  return(paste(response_type, terr_type, mass_type, sep = "_"))
}

# Function to prepare data for a specific formula
prepare_data_for_formula <- function(formula, full_data, tree) {
  # Extract variables from formula
  all_vars <- all.vars(formula)
  
  # Select relevant columns
  data_subset <- full_data[, c("species", all_vars)]
  
  # Remove rows with NAs
  data_clean <- na.omit(data_subset)
  
  # Match tree
  tree_subset <- keep.tip(tree, data_clean$species)
  data_final <- data_clean[match(tree_subset$tip.label, data_clean$species), ]
  
  return(list(data = data_final, tree = tree_subset))
}

# Function to define available predictors
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
    ),
    Colonial = list(
      name = "Colonial",
      var = "colonial_Griesser2023",
      values_func = function(data) data$colonial_Griesser2023,
      skip_conditions = "colonial|Colonial"
    )
  )
}

# Function to get remaining predictors for current formula
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

# Function to test predictors without bootstrap (fast)
test_predictors_fast <- function(base_formula, base_data, base_tree, base_aic,
                                remaining_predictors, full_data, verbose, use_normalized_mass = TRUE) {
  
  results <- list()
  
  for (pred_name in names(remaining_predictors)) {
    pred <- remaining_predictors[[pred_name]]
    
    if (verbose) message("  Testing: ", pred$name)
    
    # Track detailed information for this test
    test_details <- list(
      predictor_name = pred$name,
      predictor_var = pred$var,
      base_formula = base_formula,
      n_species_base = nrow(base_data),
      start_time = Sys.time()
    )
    
    # Add predictor values to data
    test_data <- base_data
    test_data[[pred$var]] <- pred$values_func(full_data)[match(test_data$species, full_data$species)]
    
    # Handle categorical variables
    if (pred$var == "GeographicRegion_Jetz") {
      test_data[[pred$var]] <- as.factor(test_data[[pred$var]])
    }
    
    # Track species before NA removal
    n_species_before_na <- nrow(test_data)
    
    # Remove NAs
    test_data <- test_data[!is.na(test_data[[pred$var]]), ]
    n_species_after_na <- nrow(test_data)
    
    test_details$n_species_before_na_removal <- n_species_before_na
    test_details$n_species_after_na_removal <- n_species_after_na
    test_details$n_species_lost_to_na <- n_species_before_na - n_species_after_na
    
    if (nrow(test_data) < 100) {
      if (verbose) message("    Not enough species with data (", nrow(test_data), ")")
      test_details$status <- "insufficient_data"
      test_details$n_species_final <- nrow(test_data)
      test_details$end_time <- Sys.time()
      test_details$elapsed_seconds <- as.numeric(difftime(test_details$end_time, 
                                                          test_details$start_time, 
                                                          units = "secs"))
      results[[pred_name]] <- test_details
      next
    }
    
    # Match tree
    test_tree <- keep.tip(base_tree, test_data$species)
    test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
    test_details$n_species_final <- nrow(test_data)
    
    # Try to refit base model to subset
    base_subset_converged <- TRUE
    base_subset_warnings <- NULL
    
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
    }, warning = function(w) {
      base_subset_warnings <<- c(base_subset_warnings, conditionMessage(w))
    }, error = function(e) {
      base_subset_converged <<- FALSE
      test_details$status <<- "base_model_failed"
      test_details$error_message <<- conditionMessage(e)
    })
    
    if (!base_subset_converged) {
      test_details$end_time <- Sys.time()
      test_details$elapsed_seconds <- as.numeric(difftime(test_details$end_time, 
                                                          test_details$start_time, 
                                                          units = "secs"))
      results[[pred_name]] <- test_details
      next
    }
    
    test_details$base_subset_aic <- base_subset_aic
    test_details$base_subset_warnings <- base_subset_warnings
    
    # Fit expanded model
    new_formula <- update(base_formula, paste("~ . +", pred$var))
    test_details$expanded_formula <- new_formula
    
    expanded_converged <- TRUE
    expanded_warnings <- NULL
    
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
    }, warning = function(w) {
      expanded_warnings <<- c(expanded_warnings, conditionMessage(w))
    }, error = function(e) {
      expanded_converged <<- FALSE
      test_details$status <<- "expanded_model_failed"
      test_details$error_message <<- conditionMessage(e)
    })
    
    if (!expanded_converged) {
      test_details$end_time <- Sys.time()
      test_details$elapsed_seconds <- as.numeric(difftime(test_details$end_time, 
                                                          test_details$start_time, 
                                                          units = "secs"))
      results[[pred_name]] <- test_details
      next
    }
    
    if (verbose) {
      message("    AIC improvement: ", round(aic_improvement, 2))
    }
    
    # Store all results
    test_details$status <- "successful"
    test_details$new_aic <- expanded_aic
    test_details$aic_improvement <- aic_improvement
    test_details$expanded_warnings <- expanded_warnings
    test_details$formula <- new_formula
    test_details$data <- test_data
    test_details$tree <- test_tree
    test_details$model <- expanded_fit
    test_details$coefficients <- coef(expanded_fit)
    test_details$n_parameters <- expanded_fit$d
    test_details$end_time <- Sys.time()
    test_details$elapsed_seconds <- as.numeric(difftime(test_details$end_time, 
                                                        test_details$start_time, 
                                                        units = "secs"))
    
    results[[pred_name]] <- test_details
  }
  
  return(results)
}

# Function to find best improvement
find_best_improvement <- function(iteration_results, aic_threshold) {
  best_improvement <- NULL
  best_aic_improvement <- 0
  
  for (result in iteration_results) {
    # Only consider successful tests
    if (!is.null(result$status) && result$status == "successful" &&
        !is.null(result$aic_improvement) &&
        result$aic_improvement >= aic_threshold && 
        result$aic_improvement > best_aic_improvement) {
      best_improvement <- result
      best_aic_improvement <- result$aic_improvement
    }
  }
  
  return(best_improvement)
}

# Function to run bootstrap model
run_bootstrap_model <- function(formula, data, tree, n_boot = 500, 
                               method = "logistic_MPLE", save_prefix = NULL,
                               save_matrices = TRUE, matrix_dir = NULL,
                               save_coefficient_csv = TRUE,
                               use_bootstrap_pvalues = FALSE,
                               bias_threshold_sd = 1.0) {
  
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
  boot_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4,
    boot = n_boot
  )
  
  # Get coefficient information with GUARANTEED alignment
  original_coef <- coef(original_fit)
  coef_names <- names(original_coef)
  n_coef <- length(coef_names)
  
  # Get summary information - extract by name to ensure alignment
  fit_summary <- summary(original_fit)$coefficients
  
  # CRITICAL: Extract p-values and SEs by matching names
  param_pvals <- numeric(n_coef)
  std_errors <- numeric(n_coef)
  
  for (i in 1:n_coef) {
    param_name <- coef_names[i]
    if (param_name %in% rownames(fit_summary)) {
      param_pvals[i] <- fit_summary[param_name, "p.value"]
      std_errors[i] <- fit_summary[param_name, "StdErr"]
    } else {
      warning(paste("Parameter", param_name, "not found in summary output"))
      param_pvals[i] <- NA
      std_errors[i] <- NA
    }
  }
  
  # Initialize bootstrap p-values as NA
  boot_pvals <- rep(NA, n_coef)
  
  # Extract bootstrap results
  if (!is.null(boot_fit$bootstrap)) {
    boot_matrix <- boot_fit$bootstrap
    
    # CRITICAL: Check if bootstrap matrix has column names
    # If not, we need to be very careful about order
    if (!is.null(colnames(boot_matrix))) {
      # Reorder bootstrap matrix to match coefficient order
      boot_matrix_ordered <- matrix(NA, nrow = nrow(boot_matrix), ncol = n_coef)
      colnames(boot_matrix_ordered) <- coef_names
      
      for (i in 1:n_coef) {
        param_name <- coef_names[i]
        if (param_name %in% colnames(boot_matrix)) {
          boot_matrix_ordered[, i] <- boot_matrix[, param_name]
        } else {
          # Try to match by position if names don't match
          if (i <= ncol(boot_matrix)) {
            boot_matrix_ordered[, i] <- boot_matrix[, i]
            warning(paste("Parameter", param_name, 
                         "not found in bootstrap matrix by name, using position", i))
          }
        }
      }
      boot_matrix_coef <- boot_matrix_ordered
    } else {
      # No column names - assume order matches (but warn)
      if (ncol(boot_matrix) >= n_coef) {
        boot_matrix_coef <- boot_matrix[, 1:n_coef, drop = FALSE]
        warning("Bootstrap matrix has no column names - assuming parameter order matches coefficient order")
      } else {
        stop("Bootstrap matrix has fewer columns than coefficients")
      }
    }
    
    # Save matrix if requested
    matrix_file <- NULL
    if (save_matrices && !is.null(save_prefix) && !is.null(matrix_dir)) {
      dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
      matrix_file <- file.path(matrix_dir, paste0(save_prefix, "_boot", n_boot, "_matrix.rds"))
      # Save with column names for future reference
      colnames(boot_matrix_coef) <- coef_names
      saveRDS(boot_matrix_coef, matrix_file, compress = TRUE)
    }
    
    # Calculate bootstrap statistics
    boot_means <- colMeans(boot_matrix_coef, na.rm = TRUE)
    boot_sds <- apply(boot_matrix_coef, 2, sd, na.rm = TRUE)
    boot_lower <- apply(boot_matrix_coef, 2, quantile, probs = 0.025, na.rm = TRUE)
    boot_upper <- apply(boot_matrix_coef, 2, quantile, probs = 0.975, na.rm = TRUE)
    
    # ALWAYS calculate bootstrap-based p-values (for comparison)
    for (i in 1:n_coef) {
      boot_samples <- boot_matrix_coef[, i]
      boot_samples <- boot_samples[!is.na(boot_samples)]
      if (length(boot_samples) > 0) {
        if (original_coef[i] > 0) {
          boot_pvals[i] <- 2 * min(mean(boot_samples <= 0), mean(boot_samples >= 0))
        } else if (original_coef[i] < 0) {
          boot_pvals[i] <- 2 * min(mean(boot_samples >= 0), mean(boot_samples <= 0))
        } else {
          boot_pvals[i] <- 1
        }
      } else {
        boot_pvals[i] <- NA
      }
    }
    
    # Decide which p-values to use for significance stars
    if (use_bootstrap_pvalues) {
      pvals_for_sig <- boot_pvals
    } else {
      pvals_for_sig <- param_pvals
    }
    
    # Add significance stars based on selected p-values
    Significance <- character(n_coef)
    for (i in 1:n_coef) {
      if (is.na(pvals_for_sig[i])) {
        Significance[i] <- ""
      } else if (pvals_for_sig[i] < 0.001) {
        Significance[i] <- "***"
      } else if (pvals_for_sig[i] < 0.01) {
        Significance[i] <- "**"
      } else if (pvals_for_sig[i] < 0.05) {
        Significance[i] <- "*"
      } else if (pvals_for_sig[i] < 0.1) {
        Significance[i] <- "."
      } else {
        Significance[i] <- ""
      }
    }
    
    # Calculate odds ratios - these are aligned with coefficients
    Odds_Ratio_OG <- exp(original_coef)
    Odds_Ratio_Boot <- exp(boot_means)
    OR_CI_Lower <- exp(boot_lower)
    OR_CI_Upper <- exp(boot_upper)
    
    # Calculate bias metrics
    Bias <- boot_means - original_coef
    Bias_SE_Units <- Bias / boot_sds
    Relative_Bias <- Bias / abs(original_coef)
    Relative_Bias[is.infinite(Relative_Bias)] <- NA  # Handle division by zero
    
    # Create coefficient summary with all information
    coef_summary <- data.frame(
      Parameter = coef_names,
      Estimate = original_coef,
      Boot_Mean = boot_means,
      Boot_SD = boot_sds,
      Bias = Bias,
      Bias_SE_Units = Bias_SE_Units,
      Relative_Bias = Relative_Bias,
      CI_Lower = boot_lower,
      CI_Upper = boot_upper,
      Odds_Ratio_OG = Odds_Ratio_OG,
      Odds_Ratio = Odds_Ratio_Boot,
      OR_CI_Lower = OR_CI_Lower,
      OR_CI_Upper = OR_CI_Upper,
      p_value_param = param_pvals,
      p_value_boot = boot_pvals,
      p_value = pvals_for_sig,  # The one used for significance
      Significance = Significance,
      row.names = NULL,  # Avoid row names to prevent confusion
      stringsAsFactors = FALSE
    )
    
    n_successful <- sum(complete.cases(boot_matrix_coef))
    n_converged <- sum(!is.na(boot_matrix_coef[,1]))
    convergence_rate <- n_converged / n_boot
    
  } else {
    # Fallback when bootstrap fails
    message("Note: boot_fit$bootstrap was NULL; using parametric standard errors for CIs")
    
    # Add significance stars
    Significance <- character(n_coef)
    for (i in 1:n_coef) {
      if (is.na(param_pvals[i])) {
        Significance[i] <- ""
      } else if (param_pvals[i] < 0.001) {
        Significance[i] <- "***"
      } else if (param_pvals[i] < 0.01) {
        Significance[i] <- "**"
      } else if (param_pvals[i] < 0.05) {
        Significance[i] <- "*"
      } else if (param_pvals[i] < 0.1) {
        Significance[i] <- "."
      } else {
        Significance[i] <- ""
      }
    }
    
    # Calculate CIs and odds ratios using parametric estimates
    param_ci_lower <- original_coef - 1.96 * std_errors
    param_ci_upper <- original_coef + 1.96 * std_errors
    
    coef_summary <- data.frame(
      Parameter = coef_names,
      Estimate = original_coef,
      Boot_Mean = original_coef,
      Boot_SD = std_errors,
      Bias = 0,  # No bias if no bootstrap
      Bias_SE_Units = 0,
      Relative_Bias = 0,
      CI_Lower = param_ci_lower,
      CI_Upper = param_ci_upper,
      Odds_Ratio_OG = exp(original_coef),
      Odds_Ratio = exp(original_coef),
      OR_CI_Lower = exp(param_ci_lower),
      OR_CI_Upper = exp(param_ci_upper),
      p_value_param = param_pvals,
      p_value_boot = boot_pvals,  # Will be NA
      p_value = param_pvals,
      Significance = Significance,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    matrix_file <- NULL
    n_successful <- 1
    n_converged <- 1
    convergence_rate <- 1
  }
  
  # Check for bias issues
  if (!is.null(boot_fit$bootstrap)) {
    bias_issues <- abs(coef_summary$Bias_SE_Units) > bias_threshold_sd & !is.na(coef_summary$Bias_SE_Units)
    if (any(bias_issues)) {
      message("\nWARNING: Large bootstrap bias detected:")
      for (i in which(bias_issues)) {
        message(sprintf("  %s: Estimate = %.3f, Boot_Mean = %.3f (bias = %.1f SDs)",
                       coef_summary$Parameter[i],
                       coef_summary$Estimate[i],
                       coef_summary$Boot_Mean[i],
                       coef_summary$Bias_SE_Units[i]))
      }
    }
  }
  
  # Diagnostic check: Print warning if p-value and CI disagree substantially
  for (i in 1:nrow(coef_summary)) {
    param <- coef_summary$Parameter[i]
    p_val <- coef_summary$p_value[i]
    or_lower <- coef_summary$OR_CI_Lower[i]
    or_upper <- coef_summary$OR_CI_Upper[i]
    
    if (!is.na(p_val) && !is.na(or_lower) && !is.na(or_upper)) {
      ci_excludes_1 <- (or_lower > 1) || (or_upper < 1)
      p_significant <- p_val < 0.05
      
      if (ci_excludes_1 != p_significant) {
        # Also check if parametric and bootstrap p-values disagree
        param_sig <- coef_summary$p_value_param[i] < 0.05
        boot_sig <- !is.na(coef_summary$p_value_boot[i]) && coef_summary$p_value_boot[i] < 0.05
        
        message(paste("\nWARNING: Parameter", param, "has inconsistent inference:"))
        message(sprintf("  Parametric p-value: %.4f %s", 
                       coef_summary$p_value_param[i],
                       ifelse(param_sig, "(significant)", "(not significant)")))
        if (!is.na(coef_summary$p_value_boot[i])) {
          message(sprintf("  Bootstrap p-value: %.4f %s", 
                         coef_summary$p_value_boot[i],
                         ifelse(boot_sig, "(significant)", "(not significant)")))
        }
        message(sprintf("  OR 95%% CI: [%.3f, %.3f] %s",
                       or_lower, or_upper,
                       ifelse(ci_excludes_1, "(excludes 1)", "(includes 1)")))
      }
    }
  }
  
  # Save coefficient summary if requested
  if (save_coefficient_csv && !is.null(save_prefix) && !is.null(matrix_dir)) {
    csv_filename <- paste0("coefficients_OddsRatios_", save_prefix, "_boot", n_boot, ".csv")
    write.csv(coef_summary, 
              file.path(matrix_dir, csv_filename), 
              row.names = FALSE)
  }
  
  return(list(
    fit = original_fit,
    bootstrap_fit = boot_fit,
    coefficients = coef_summary,
    bootstrap_matrix_file = matrix_file,
    n_successful_boots = n_successful,
    n_converged = n_converged,
    convergence_rate = convergence_rate,
    bias_threshold_sd = bias_threshold_sd,
    inference_method = ifelse(use_bootstrap_pvalues, "bootstrap", "parametric")
  ))
}

# Function to save detailed expansion log
save_detailed_expansion_log <- function(expansion_path, output_dir) {
  
  # Create a detailed log of all tests performed
  detailed_log <- data.frame()
  
  # Check if there are any iterations
  if (length(expansion_path$iterations) > 0) {
    for (iter_num in 1:length(expansion_path$iterations)) {
      iter <- expansion_path$iterations[[iter_num]]
      
      # Get base formula for this iteration
      if (iter_num == 1) {
        base_formula_str <- paste(deparse(expansion_path$initial_formula), collapse = " ")
      } else {
        # Get the formula after the previous iteration's selected predictor
        prev_iter <- expansion_path$iterations[[iter_num - 1]]
        # Find the selected test summary
        for (test_name in names(prev_iter$test_summary)) {
          if (prev_iter$test_summary[[test_name]]$predictor == prev_iter$selected_predictor) {
            base_formula_str <- prev_iter$test_summary[[test_name]]$formula
            break
          }
        }
      }
      
      # Process each test in this iteration (test_summary is a list)
      for (test_name in names(iter$test_summary)) {
        test_summary <- iter$test_summary[[test_name]]
        row_data <- data.frame(
          Model_Name = expansion_path$model_name,
          Analysis_Type = expansion_path$analysis_type,
          Iteration = iter_num,
          Base_Formula = base_formula_str,
          Predictor_Tested = test_summary$predictor,
          Expanded_Formula = test_summary$formula,
          Status = test_summary$status,
          N_Species = test_summary$n_species_final,
          AIC_Improvement = test_summary$aic_improvement,
          Selected = !is.na(iter$selected_predictor) && test_summary$predictor == iter$selected_predictor,
          stringsAsFactors = FALSE
        )
        
        detailed_log <- rbind(detailed_log, row_data)
      }
    }
  }
  
  # Save the detailed log
  write.csv(detailed_log,
            file.path(output_dir, "detailed_expansion_log.csv"),
            row.names = FALSE)
  
  # If detailed results are available, save additional information
  if (length(expansion_path$iterations) > 0 && 
      !is.null(expansion_path$iterations[[1]]$all_test_results)) {
    # Create an even more detailed log with species counts at each stage
    ultra_detailed_log <- data.frame()
    
    for (iter_num in 1:length(expansion_path$iterations)) {
      iter <- expansion_path$iterations[[iter_num]]
      
      if (!is.null(iter$all_test_results)) {
        for (pred_name in names(iter$all_test_results)) {
          test_result <- iter$all_test_results[[pred_name]]
          
          row_data <- data.frame(
            Iteration = iter_num,
            Predictor = test_result$predictor_name,
            Variable = test_result$predictor_var,
            N_Species_Base = test_result$n_species_base,
            N_Species_Before_NA = test_result$n_species_before_na_removal,
            N_Species_After_NA = test_result$n_species_after_na_removal,
            N_Species_Lost_To_NA = test_result$n_species_lost_to_na,
            N_Species_Final = test_result$n_species_final,
            Status = test_result$status,
            Base_Subset_AIC = if (!is.null(test_result$base_subset_aic)) 
              test_result$base_subset_aic else NA,
            Expanded_AIC = if (!is.null(test_result$new_aic)) 
              test_result$new_aic else NA,
            AIC_Improvement = if (!is.null(test_result$aic_improvement)) 
              test_result$aic_improvement else NA,
            N_Parameters = if (!is.null(test_result$n_parameters)) 
              test_result$n_parameters else NA,
            Selected = !is.na(iter$selected_predictor) && test_result$predictor_name == iter$selected_predictor,
            Elapsed_Seconds = if (!is.null(test_result$elapsed_seconds)) 
              test_result$elapsed_seconds else NA,
            Error_Message = if (!is.null(test_result$error_message)) 
              test_result$error_message else "",
            stringsAsFactors = FALSE
          )
          
          ultra_detailed_log <- rbind(ultra_detailed_log, row_data)
        }
      }
    }
    
    if (nrow(ultra_detailed_log) > 0) {
      write.csv(ultra_detailed_log,
                file.path(output_dir, "ultra_detailed_expansion_log.csv"),
                row.names = FALSE)
    }
    
    # Save complete test results as RDS for future reference
    all_test_results <- lapply(expansion_path$iterations, function(iter) {
      iter$all_test_results
    })
    saveRDS(all_test_results,
            file.path(output_dir, "all_test_results.rds"))
  }
}

# Function to save model-specific results
save_model_results <- function(expansion_path, output_dir, n_bootstrap) {
  
  # Save complete expansion path as RDS
  saveRDS(expansion_path, file.path(output_dir, "expansion_path.rds"))
  
  # Save detailed expansion logs
  save_detailed_expansion_log(expansion_path, output_dir)
  
  # Create expansion history CSV
  if (length(expansion_path$iterations) > 0) {
    expansion_history <- data.frame()
    
    for (i in 1:length(expansion_path$iterations)) {
      iter <- expansion_path$iterations[[i]]
      expansion_history <- rbind(expansion_history, data.frame(
        Iteration = i,
        Selected_Predictor = iter$selected_predictor,
        AIC_Improvement = iter$aic_improvement,
        New_AIC = iter$new_aic,
        N_Species = iter$n_species,
        N_Predictors_Tested = length(iter$predictors_tested),
        stringsAsFactors = FALSE
      ))
    }
    
    write.csv(expansion_history, 
              file.path(output_dir, "expansion_history.csv"),
              row.names = FALSE)
  }
  
  # Save final model coefficients with CIs
  if (!is.null(expansion_path$final_bootstrap)) {
    final_coefs <- expansion_path$final_bootstrap$coefficients
    
    # The bootstrap results already include odds ratios, p-values, bias metrics, and significance
    # Just save as is
    write.csv(final_coefs,
              file.path(output_dir, paste0("final_coefficients_boot", n_bootstrap, ".csv")),
              row.names = FALSE)
    
    # Add inference method and bias info
    inference_info <- data.frame(
      Inference_Method = expansion_path$final_bootstrap$inference_method,
      Bootstrap_Iterations = n_bootstrap,
      Convergence_Rate = expansion_path$final_bootstrap$convergence_rate,
      Bias_Threshold_SD = expansion_path$final_bootstrap$bias_threshold_sd,
      Parameters_With_High_Bias = sum(abs(final_coefs$Bias_SE_Units) > expansion_path$final_bootstrap$bias_threshold_sd, na.rm = TRUE)
    )
    write.csv(inference_info,
              file.path(output_dir, "inference_method.csv"),
              row.names = FALSE)
  }
  
  # Save model comparison with all iterations
  model_comparison <- data.frame()
  
  # Add initial model
  model_comparison <- rbind(model_comparison, data.frame(
    Iteration = 0,
    Model = "Initial",
    Formula = paste(deparse(expansion_path$initial_formula), collapse = " "),
    AIC = expansion_path$initial_aic,
    N_Species = expansion_path$initial_n_species,
    N_Parameters = length(coef(expansion_path$initial_bootstrap$fit)),
    AIC_Improvement_From_Previous = 0,
    Total_AIC_Improvement = 0,
    Predictor_Added = NA,
    stringsAsFactors = FALSE
  ))
  
  # Add each iteration
  current_formula <- expansion_path$initial_formula
  cumulative_improvement <- 0
  
  if (length(expansion_path$iterations) > 0) {
    for (i in 1:length(expansion_path$iterations)) {
    iter <- expansion_path$iterations[[i]]
    if (!is.na(iter$aic_improvement)) {
      cumulative_improvement <- cumulative_improvement + iter$aic_improvement
    }
    
    # Find the formula for this iteration if a predictor was selected
    if (!is.na(iter$selected_predictor)) {
      for (test_name in names(iter$test_summary)) {
        if (iter$test_summary[[test_name]]$predictor == iter$selected_predictor) {
          current_formula <- iter$test_summary[[test_name]]$formula
          break
        }
      }
    }
    # If no predictor was selected, formula stays the same
    
    model_comparison <- rbind(model_comparison, data.frame(
      Iteration = i,
      Model = paste0("Iteration_", i),
      Formula = if (!is.na(iter$selected_predictor)) current_formula else paste(deparse(current_formula), collapse = " "),
      AIC = iter$new_aic,
      N_Species = iter$n_species,
      N_Parameters = NA,  # Could extract from test results if needed
      AIC_Improvement_From_Previous = if (!is.na(iter$aic_improvement)) iter$aic_improvement else 0,
      Total_AIC_Improvement = cumulative_improvement,
      Predictor_Added = if (!is.na(iter$selected_predictor)) iter$selected_predictor else "None - No improvement",
      stringsAsFactors = FALSE
    ))
    }
  }
  
  # Add final model if different from last iteration
  if (length(expansion_path$iterations) == 0 || 
      expansion_path$final_aic != expansion_path$iterations[[length(expansion_path$iterations)]]$new_aic) {
    model_comparison <- rbind(model_comparison, data.frame(
      Iteration = length(expansion_path$iterations) + 1,
      Model = "Final",
      Formula = paste(deparse(expansion_path$final_formula), collapse = " "),
      AIC = expansion_path$final_aic,
      N_Species = expansion_path$final_n_species,
      N_Parameters = length(coef(expansion_path$final_bootstrap$fit)),
      AIC_Improvement_From_Previous = 0,
      Total_AIC_Improvement = expansion_path$initial_aic - expansion_path$final_aic,
      Predictor_Added = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(model_comparison,
            file.path(output_dir, "model_comparison.csv"),
            row.names = FALSE)
}

# Function to create summary across all expansion paths
create_summary_results <- function(all_expansion_paths, output_dir, n_bootstrap) {
  
  # Create summary table
  summary_table <- data.frame()
  
  for (path_name in names(all_expansion_paths)) {
    path <- all_expansion_paths[[path_name]]
    
    summary_table <- rbind(summary_table, data.frame(
      Model_Name = path_name,
      Analysis_Type = path$analysis_type,
      Initial_Formula = paste(deparse(path$initial_formula), collapse = " "),
      Final_Formula = paste(deparse(path$final_formula), collapse = " "),
      Initial_AIC = path$initial_aic,
      Final_AIC = path$final_aic,
      Total_AIC_Improvement = path$initial_aic - path$final_aic,
      N_Iterations = path$n_iterations,
      N_Predictors_Added = length(path$iterations),
      Status = path$status,
      Initial_N_Species = path$initial_n_species,
      Final_N_Species = path$final_n_species,
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(summary_table,
            file.path(output_dir, paste0("all_paths_summary_boot", n_bootstrap, ".csv")),
            row.names = FALSE)
  
  # Create convergence analysis
  convergence_analysis <- analyze_convergence(all_expansion_paths)
  if (!is.null(convergence_analysis)) {
    write.csv(convergence_analysis,
              file.path(output_dir, "convergence_analysis.csv"),
              row.names = FALSE)
  }
  
  # Create predictor importance summary
  predictor_importance <- calculate_predictor_importance(all_expansion_paths)
  write.csv(predictor_importance,
            file.path(output_dir, "predictor_importance.csv"),
            row.names = FALSE)
  
  # Save complete results object
  saveRDS(all_expansion_paths,
          file.path(output_dir, paste0("all_expansion_paths_boot", n_bootstrap, ".rds")))
}

# Function to analyze convergence of paths
analyze_convergence <- function(all_expansion_paths) {
  # Group paths by analysis type
  by_analysis <- list()
  
  for (path_name in names(all_expansion_paths)) {
    path <- all_expansion_paths[[path_name]]
    analysis_type <- path$analysis_type
    
    if (!analysis_type %in% names(by_analysis)) {
      by_analysis[[analysis_type]] <- list()
    }
    
    by_analysis[[analysis_type]][[path_name]] <- path
  }
  
  # Check for convergence within each analysis type
  convergence_results <- data.frame()
  
  for (analysis_type in names(by_analysis)) {
    paths <- by_analysis[[analysis_type]]
    
    if (length(paths) > 1) {
      # Compare final formulas
      final_formulas <- sapply(paths, function(p) paste(deparse(p$final_formula), collapse = " "))
      unique_formulas <- unique(final_formulas)
      
      convergence_results <- rbind(convergence_results, data.frame(
        Analysis_Type = analysis_type,
        N_Paths = length(paths),
        N_Unique_Final_Models = length(unique_formulas),
        Converged = length(unique_formulas) == 1,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(if(nrow(convergence_results) > 0) convergence_results else NULL)
}

# Function to calculate predictor importance
calculate_predictor_importance <- function(all_expansion_paths) {
  predictor_counts <- list()
  
  for (path_name in names(all_expansion_paths)) {
    path <- all_expansion_paths[[path_name]]
    
    for (iter in path$iterations) {
      pred <- iter$selected_predictor
      # Skip if no predictor was selected (NA)
      if (!is.na(pred)) {
        if (!pred %in% names(predictor_counts)) {
          predictor_counts[[pred]] <- 0
        }
        predictor_counts[[pred]] <- predictor_counts[[pred]] + 1
      }
    }
  }
  
  # Handle case where no predictors were selected
  if (length(predictor_counts) == 0) {
    return(data.frame(
      Predictor = character(0),
      Times_Selected = numeric(0),
      Proportion = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  importance_df <- data.frame(
    Predictor = names(predictor_counts),
    Times_Selected = unlist(predictor_counts),
    Proportion = unlist(predictor_counts) / sum(unlist(predictor_counts)),
    stringsAsFactors = FALSE
  )
  
  return(importance_df[order(importance_df$Times_Selected, decreasing = TRUE), ])
}

# # Example usage:
# # Using RDS input:
# results <- run_unified_stepwise_expansion(
#   input_source = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
#   top_n_models = 2,
#   n_bootstrap = 500,
#   use_bootstrap_pvalues = TRUE,  # For consistent p-values and CIs
#   bias_threshold_sd = 1.0,         # Warn if bias > 1 SD
#   use_normalized_mass = TRUE,      # Replace logMass_AVONET with logMass_normalized
#   output_base_dir = "Outputs/PhyloglmResults/StepwiseIterative_Top2Models_boot500_aicThresh0.5",save_boot_matrices = T, save_coefficient_csv = T,verbose = T, data_path = "Data_R_2025-07-23.csv", tree_path = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex", detailed_logging = TRUE, aic_threshold = 0.5
# )

# Using formula input:
formulas <- list(
  FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized,
  FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized,
  HighConfidence_Coop ~ FemaleSong_Agg01 * Territory_12vs3,
  HighConfidence_Coop ~ FemaleSong_Agg01 * TerritorialityWeakVsStrong
)
results <- run_unified_stepwise_expansion(
  input_source = formulas,
  data_path = "Data_R_2025-07-23.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  n_bootstrap = 500,
  save_coefficient_csv = TRUE,
  use_bootstrap_pvalues = TRUE,  # Use parametric p-values
  bias_threshold_sd = 1,         # More lenient bias threshold
  use_normalized_mass = TRUE, top_n_models =  1, output_base_dir = "Outputs/PhyloglmResults/StepwiseIterative_Top1Model_MassNorm_AddPolygyny_boot500_aicThresh2.0", aic_threshold = 2.0, save_boot_matrices = T, verbose = T
)

