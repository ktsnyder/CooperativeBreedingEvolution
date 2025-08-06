# Stepwise Model Expansion Framework for PhyloGLM
# This framework takes the best models from existing analyses and iteratively
# adds predictors to find optimal "megamodels"

library(phylolm)
library(dplyr)
library(ape)

source("phyloglm_framework/batch_runner_helpers.R")  # Contains calculate_effect_sizes and bootstrap functions

#' Define candidate predictors for stepwise addition
#' @return List of predictor information
get_candidate_predictors <- function() {
  list(
    # Continuous predictors
    continuous = list(
      Territory_numeric = list(
        name = "Territory",
        type = "continuous",
        description = "Territory strength (1-3)",
        transform = function(x) as.numeric(x)
      ),
      Migration_numeric = list(
        name = "Migration_AVONET", 
        type = "continuous",
        description = "Migration (1-3)",
        transform = function(x) as.numeric(x)
      ),
      Latitude_abs = list(
        name = "Centroid.Latitude_AVONET",
        type = "continuous", 
        description = "Absolute latitude",
        transform = function(x) abs(x),
        check_redundancy = "absLat"  # Check if already in model as abs(Lat)
      ),
      PlumDim = list(
        name = "logMaleFemalePlumageDiffAbs",
        type = "continuous",
        description = "Plumage dimorphism",
        transform = NULL,
        check_redundancy = "PlumDim"
      ),
      WingDim = list(
        name = "PercentAbsLogWingDimorphism",
        type = "continuous",
        description = "Wing dimorphism", 
        transform = NULL,
        check_redundancy = "WingDim"
      )
    ),
    
    # Binary predictors
    binary = list(
      GeographicRegion = list(
        name = "GeographicRegion_Jetz",
        type = "binary",
        description = "Geographic region (Holarctic=0, Tropical=1)",
        transform = NULL,
        check_redundancy = "Region"
      ),
      FamilialLiving = list(
        name = "Griesser2017FamilialLiving",
        type = "binary",
        description = "Familial living",
        transform = NULL,
        check_redundancy = c("Fam", "FL"),
        note = "Reduced species coverage"
      )
    )
  )
}

#' Extract predictors already in a model formula
#' @param formula Model formula
#' @return Character vector of predictor names
extract_model_predictors <- function(formula) {
  # Convert formula to character and parse
  form_char <- as.character(formula)[3]  # Right-hand side
  
  # Remove interactions and spaces
  form_char <- gsub("\\*|:", "+", form_char)
  form_char <- gsub("\\s+", "", form_char)
  
  # Split by + and clean
  predictors <- unlist(strsplit(form_char, "\\+"))
  predictors <- unique(predictors)
  predictors <- predictors[predictors != ""]
  
  # Also check for transformed variables
  # e.g., abs(Centroid.Latitude_AVONET) -> Centroid.Latitude_AVONET
  predictors_clean <- gsub("^abs\\((.+)\\)$", "\\1", predictors)
  predictors_clean <- gsub("^log\\((.+)\\)$", "\\1", predictors_clean)
  
  return(unique(c(predictors, predictors_clean)))
}

#' Check if a predictor is redundant with existing model
#' @param candidate Candidate predictor info
#' @param existing_predictors Vector of existing predictor names
#' @param model_name Name of the model being expanded
#' @return Logical indicating if redundant
is_redundant_predictor <- function(candidate, existing_predictors, model_name) {
  # Check if already in model
  if (candidate$name %in% existing_predictors) {
    return(TRUE)
  }
  
  # Check redundancy patterns
  if (!is.null(candidate$check_redundancy)) {
    patterns <- candidate$check_redundancy
    for (pattern in patterns) {
      if (any(grepl(pattern, existing_predictors)) || grepl(pattern, model_name)) {
        return(TRUE)
      }
    }
  }
  
  return(FALSE)
}

#' Perform stepwise model expansion
#' @param base_result Result object from original analysis
#' @param data Full dataset
#' @param tree Phylogenetic tree
#' @param candidates List of candidate predictors
#' @param aic_threshold Minimum AIC improvement to keep predictor (default 2)
#' @param n_bootstrap Number of bootstrap iterations for final model
#' @param verbose Print progress messages
#' @return List with expansion results
stepwise_expand_model <- function(base_result, 
                                 data,
                                 tree,
                                 candidates = get_candidate_predictors(),
                                 aic_threshold = 2,
                                 n_bootstrap = 100,
                                 verbose = TRUE) {
  
  # Get best model info
  best_model_name <- base_result$comparison$comparison$Model[1]
  best_model <- base_result$models$models[[best_model_name]]
  base_formula <- formula(best_model)
  base_aic <- base_result$comparison$comparison$AIC[1]
  
  if (verbose) {
    cat("\nStarting stepwise expansion from base model:\n")
    cat("Formula:", deparse(base_formula), "\n")
    cat("AIC:", base_aic, "\n")
    cat("N species:", nrow(best_model$X), "\n\n")
  }
  
  # Extract current predictors
  current_predictors <- extract_model_predictors(base_formula)
  response_var <- as.character(base_formula)[2]
  
  # Get ALL variables from the base model's data to ensure we have everything
  base_data_vars <- names(best_model$data)
  
  # Prepare data with all potential predictors
  all_vars <- unique(c(response_var, current_predictors, base_data_vars,
                      sapply(c(candidates$continuous, candidates$binary), 
                             function(x) x$name)))
  
  # Check which variables exist
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    warning("Variables not found in data: ", paste(missing_vars, collapse = ", "))
    all_vars <- intersect(all_vars, names(data))
  }
  
  # Subset and match with tree
  data_subset <- data[, all_vars, drop = FALSE]
  
  # Handle species column
  if (!"species" %in% names(data)) {
    if ("Species" %in% names(data)) {
      data_subset$species <- data$Species
    } else if ("X" %in% names(data)) {
      # Sometimes the first column is species names
      data_subset$species <- data$X
    } else if (!is.null(rownames(data)) && rownames(data)[1] != "1") {
      data_subset$species <- rownames(data)
    } else {
      stop("Cannot find species names in data")
    }
  } else {
    data_subset$species <- data$species
  }
  
  # Match species between data and tree
  shared_species <- intersect(data_subset$species, tree$tip.label)
  if (verbose) {
    cat("Species in data:", length(unique(data_subset$species)), "\n")
    cat("Species in tree:", length(tree$tip.label), "\n")
    cat("Shared species:", length(shared_species), "\n")
  }
  
  data_subset <- data_subset[data_subset$species %in% shared_species, ]
  tree_subset <- keep.tip(tree, shared_species)
  
  # Initialize tracking
  expansion_history <- list()
  current_formula <- base_formula
  current_aic <- base_aic
  current_model <- best_model
  improvements <- data.frame()
  
  # Test each candidate predictor
  all_candidates <- c(candidates$continuous, candidates$binary)
  
  for (i in seq_along(all_candidates)) {
    candidate <- all_candidates[[i]]
    
    # Check redundancy
    if (is_redundant_predictor(candidate, current_predictors, base_result$config$name)) {
      if (verbose) {
        cat("Skipping", candidate$name, "- redundant with existing predictors\n")
      }
      next
    }
    
    # Apply transformation if needed
    if (!is.null(candidate$transform)) {
      data_subset[[paste0(candidate$name, "_transformed")]] <- 
        candidate$transform(data_subset[[candidate$name]])
      pred_name <- paste0(candidate$name, "_transformed")
    } else {
      pred_name <- candidate$name
    }
    
    # Check if variable exists and has variation
    if (!(pred_name %in% names(data_subset))) {
      if (verbose) cat("Skipping", candidate$name, "- not in data\n")
      next
    }
    
    # Create new formula by adding predictor
    new_formula <- update(current_formula, paste("~ . +", pred_name))
    
    # Prepare data for this model
    model_vars <- all.vars(new_formula)
    model_data <- data_subset[, model_vars, drop = FALSE]
    model_data$species <- data_subset$species
    
    # Remove rows with any missing values
    complete_rows <- complete.cases(model_data)
    model_data <- model_data[complete_rows, ]
    model_species <- model_data$species
    
    # Match tree
    model_tree <- keep.tip(tree_subset, intersect(tree_subset$tip.label, model_species))
    model_data <- model_data[model_data$species %in% model_tree$tip.label, ]
    
    # Order data to match tree - CRITICAL for phyloglm
    model_data <- model_data[match(model_tree$tip.label, model_data$species), ]
    
    # Remove any NA rows from failed matching
    if (any(is.na(model_data$species))) {
      model_data <- model_data[!is.na(model_data$species), ]
    }
    
    # Final check that row names match tree tips
    rownames(model_data) <- model_data$species
    
    if (verbose) {
      cat("\nTesting addition of", candidate$name, 
          "( n =", nrow(model_data), "species )\n")
    }
    
    # Double-check data/tree alignment
    if (!all(rownames(model_data) == model_tree$tip.label)) {
      warning("Data and tree tips are not aligned properly")
      next
    }
    
    # Fit new model
    tryCatch({
      new_model <- phyloglm(
        formula = new_formula,
        data = model_data,
        phy = model_tree,
        method = "logistic_MPLE",
        btol = 30,  # Reduced tolerance
        log.alpha.bound = 4  # Reduced bound
      )
      
      new_aic <- -2 * new_model$logLik + 2 * new_model$d
      aic_improvement <- current_aic - new_aic
      
      # Store result
      test_result <- data.frame(
        predictor = candidate$name,
        description = candidate$description,
        type = candidate$type,
        n_species = nrow(model_data),
        base_aic = current_aic,
        new_aic = new_aic,
        aic_improvement = aic_improvement,
        kept = aic_improvement > aic_threshold,
        stringsAsFactors = FALSE
      )
      
      improvements <- rbind(improvements, test_result)
      
      if (verbose) {
        cat("  AIC:", round(new_aic, 2), 
            "| Improvement:", round(aic_improvement, 2),
            if(aic_improvement > aic_threshold) "*** KEPT ***" else "", "\n")
      }
      
      # Keep if improvement exceeds threshold
      if (aic_improvement > aic_threshold) {
        current_formula <- new_formula
        current_aic <- new_aic
        current_model <- new_model
        current_predictors <- c(current_predictors, pred_name)
        
        expansion_history[[length(expansion_history) + 1]] <- list(
          step = length(expansion_history) + 1,
          added = candidate$name,
          formula = new_formula,
          aic = new_aic,
          improvement = aic_improvement,
          n_species = nrow(model_data)
        )
      }
      
    }, error = function(e) {
      if (verbose) {
        cat("  Error fitting model:", e$message, "\n")
      }
      
      test_result <- data.frame(
        predictor = candidate$name,
        description = candidate$description,
        type = candidate$type,
        n_species = nrow(model_data),
        base_aic = current_aic,
        new_aic = NA,
        aic_improvement = NA,
        kept = FALSE,
        error = e$message,
        stringsAsFactors = FALSE
      )
      
      improvements <- rbind(improvements, test_result)
    })
  }
  
  # If we found improvements, bootstrap the final model
  final_results <- list(
    base_model = best_model,
    base_formula = base_formula,
    base_aic = base_aic,
    final_model = current_model,
    final_formula = current_formula,
    final_aic = current_aic,
    total_improvement = base_aic - current_aic,
    expansion_history = expansion_history,
    improvement_tests = improvements,
    n_predictors_added = length(expansion_history)
  )
  
  if (length(expansion_history) > 0 && n_bootstrap > 0) {
    if (verbose) {
      cat("\nBootstrapping final model with", n_bootstrap, "iterations...\n")
    }
    
    # Get data for final model
    final_vars <- all.vars(current_formula)
    final_data <- data_subset[, c(final_vars, "species"), drop = FALSE]
    final_data <- final_data[complete.cases(final_data), ]
    final_tree <- keep.tip(tree_subset, intersect(tree_subset$tip.label, final_data$species))
    final_data <- final_data[final_data$species %in% final_tree$tip.label, ]
    final_data <- final_data[match(final_tree$tip.label, final_data$species), ]
    
    # For now, skip bootstrap and just calculate effects
    # Bootstrap functionality can be added later
    final_results$bootstrap <- NULL
    
    # Calculate effects without bootstrap CIs
    final_results$effects <- calculate_effect_sizes(
      current_model,
      model_name = "Expanded_Model"
    )
  }
  
  class(final_results) <- c("stepwise_phyloglm", "list")
  return(final_results)
}

#' Print method for stepwise results
#' @param x Stepwise result object
#' @param ... Additional arguments
print.stepwise_phyloglm <- function(x, ...) {
  cat("Stepwise PhyloGLM Expansion Results\n")
  cat("===================================\n\n")
  
  cat("Base model:\n")
  cat("  Formula:", deparse(x$base_formula), "\n")
  cat("  AIC:", round(x$base_aic, 2), "\n\n")
  
  if (x$n_predictors_added > 0) {
    cat("Final model:\n")
    cat("  Formula:", deparse(x$final_formula), "\n")
    cat("  AIC:", round(x$final_aic, 2), "\n")
    cat("  Total AIC improvement:", round(x$total_improvement, 2), "\n\n")
    
    cat("Predictors added:\n")
    for (step in x$expansion_history) {
      cat("  Step", step$step, ":", step$added, 
          "(AIC improvement:", round(step$improvement, 2), ")\n")
    }
  } else {
    cat("No predictors improved the model beyond threshold.\n")
  }
  
  cat("\nAll predictors tested:\n")
  if (nrow(x$improvement_tests) > 0) {
    # Check which columns exist
    available_cols <- intersect(c("predictor", "type", "aic_improvement", "kept", "error"), 
                               names(x$improvement_tests))
    if (length(available_cols) > 0) {
      print(x$improvement_tests[, available_cols])
    } else {
      print(x$improvement_tests)
    }
  } else {
    cat("No predictors were tested.\n")
  }
}

#' Run stepwise expansion for multiple analyses
#' @param results_file Path to combined results RDS file
#' @param data_path Path to data CSV
#' @param tree_path Path to tree file
#' @param target_analyses Specific analyses to expand (NULL for all)
#' @param output_dir Output directory for results
#' @param ... Additional arguments passed to stepwise_expand_model
run_stepwise_batch <- function(results_file,
                              data_path,
                              tree_path, 
                              target_analyses = NULL,
                              output_dir = "Outputs/PhyloglmResults/stepwise_expansion",
                              ...) {
  
  # Load data
  cat("Loading data and results...\n")
  results <- readRDS(results_file)
  data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Determine which analyses to run
  if (is.null(target_analyses)) {
    # Default: one FS->CB and one CB->FS analysis
    target_analyses <- c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass")
  }
  
  # Run stepwise expansion for each target
  expansion_results <- list()
  
  for (analysis_name in target_analyses) {
    if (!(analysis_name %in% names(results))) {
      warning("Analysis not found:", analysis_name)
      next
    }
    
    cat("\n\n========================================\n")
    cat("Expanding:", analysis_name, "\n")
    cat("========================================\n")
    
    base_result <- results[[analysis_name]]
    
    # Run expansion
    expanded <- stepwise_expand_model(
      base_result = base_result,
      data = data,
      tree = tree,
      ...
    )
    
    # Save individual result
    saveRDS(expanded, file.path(output_dir, paste0(analysis_name, "_expanded.rds")))
    
    # Save summary
    sink(file.path(output_dir, paste0(analysis_name, "_summary.txt")))
    print(expanded)
    sink()
    
    expansion_results[[analysis_name]] <- expanded
  }
  
  # Save all results
  saveRDS(expansion_results, file.path(output_dir, "all_expansion_results.rds"))
  
  cat("\n\nStepwise expansion complete. Results saved to:", output_dir, "\n")
  
  return(expansion_results)
}