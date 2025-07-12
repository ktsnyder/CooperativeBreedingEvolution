# Enhanced Data Preparation for PhyloGLM Framework
# Handles variable type conversions, missing data, and tree matching

library(ape)
library(dplyr)
source("phyloglm_framework/variable_classification.R")

#' Prepare data for phyloglm analysis with automatic variable handling
#' 
#' @param data Original data frame
#' @param tree Original phylogenetic tree
#' @param config Analysis configuration list containing:
#'   - response: response variable name
#'   - predictors: vector of predictor names
#'   - controls: vector of control variable names
#'   - transformations: list of transformations to apply
#' @param min_species Minimum number of species required for analysis
#' @param remove_outliers Whether to remove outliers (beyond 4 SD)
#' @return List containing prepared data, tree, and metadata
prepare_analysis_data <- function(data, tree, config, 
                                 min_species = 50,
                                 remove_outliers = FALSE) {
  
  # Extract variables from config
  response_var <- config$response
  predictor_vars <- config$predictors
  control_vars <- config$controls
  all_vars <- unique(c(response_var, predictor_vars, control_vars))
  
  # Initial checks
  missing_vars <- setdiff(all_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Ensure species column exists
  if (!"species" %in% colnames(data)) {
    if ("Species" %in% colnames(data)) {
      data$species <- data$Species
    } else if (rownames(data)[1] != "1") {
      data$species <- rownames(data)
    } else {
      stop("No species identifier found in data")
    }
  }
  
  # Remove duplicate species
  if (any(duplicated(data$species))) {
    message(paste("Removing", sum(duplicated(data$species)), "duplicate species"))
    data <- data[!duplicated(data$species), ]
  }
  
  # Get variable classifications
  classifications <- get_variable_classifications()
  
  # Detect and prepare variables
  variable_info <- list()
  for (var in all_vars) {
    var_type <- detect_variable_type(data, var, classifications)
    variable_info[[var]] <- list(
      type = var_type,
      original_name = var,
      label = get_variable_label(var)
    )
    
    # Apply transformations if specified
    if (!is.null(config$transformations) && var %in% names(config$transformations)) {
      trans <- config$transformations[[var]]
      if (trans == "abs") {
        data[[var]] <- abs(data[[var]])
        variable_info[[var]]$transformation <- "absolute value"
      } else if (trans == "log") {
        # Check for non-positive values
        if (any(data[[var]] <= 0, na.rm = TRUE)) {
          data[[var]] <- log(data[[var]] + 1)
          variable_info[[var]]$transformation <- "log(x + 1)"
        } else {
          data[[var]] <- log(data[[var]])
          variable_info[[var]]$transformation <- "log"
        }
      } else if (trans == "sqrt") {
        data[[var]] <- sqrt(data[[var]])
        variable_info[[var]]$transformation <- "square root"
      } else if (trans == "scale") {
        data[[var]] <- scale(data[[var]])[, 1]
        variable_info[[var]]$transformation <- "standardized"
      }
    }
    
    # Prepare variable based on type
    data[[var]] <- prepare_variable(data, var, var_type)
  }
  
  # Filter to complete cases
  data_complete <- data[complete.cases(data[, all_vars]), ]
  n_removed <- nrow(data) - nrow(data_complete)
  
  if (n_removed > 0) {
    message(paste("Removed", n_removed, "species with missing data"))
  }
  
  # Remove outliers if requested
  if (remove_outliers) {
    continuous_vars <- all_vars[sapply(all_vars, function(v) {
      variable_info[[v]]$type %in% c("continuous", "continuous_special")
    })]
    
    outliers <- rep(FALSE, nrow(data_complete))
    for (var in continuous_vars) {
      if (is.numeric(data_complete[[var]])) {
        var_mean <- mean(data_complete[[var]], na.rm = TRUE)
        var_sd <- sd(data_complete[[var]], na.rm = TRUE)
        outliers <- outliers | abs(data_complete[[var]] - var_mean) > 4 * var_sd
      }
    }
    
    if (sum(outliers) > 0) {
      message(paste("Removing", sum(outliers), "outlier species"))
      data_complete <- data_complete[!outliers, ]
    }
  }
  
  # Match tree and data
  tree_species <- tree$tip.label
  data_species <- data_complete$species
  
  # Find common species
  common_species <- intersect(tree_species, data_species)
  
  if (length(common_species) < min_species) {
    warning(paste("Only", length(common_species), "species in common between tree and data"))
    if (length(common_species) < 10) {
      stop("Too few species for meaningful analysis")
    }
  }
  
  # Prune tree and filter data
  tree_pruned <- drop.tip(tree, setdiff(tree_species, common_species))
  data_final <- data_complete[data_complete$species %in% common_species, ]
  
  # Ensure row names match tree tips
  rownames(data_final) <- data_final$species
  data_final <- data_final[tree_pruned$tip.label, ]
  
  # Create summary statistics
  summary_stats <- data.frame(
    variable = all_vars,
    type = sapply(all_vars, function(v) variable_info[[v]]$type),
    n_complete = sapply(all_vars, function(v) sum(!is.na(data_final[[v]]))),
    n_missing_original = sapply(all_vars, function(v) sum(is.na(data[[v]]))),
    prop_missing_original = sapply(all_vars, function(v) mean(is.na(data[[v]]))),
    stringsAsFactors = FALSE
  )
  
  # Add variable-specific summaries
  for (var in all_vars) {
    if (variable_info[[var]]$type == "binary") {
      if (is.numeric(data_final[[var]])) {
        summary_stats[summary_stats$variable == var, "prop_1s"] <- 
          mean(data_final[[var]] == 1, na.rm = TRUE)
      }
    } else if (variable_info[[var]]$type %in% c("continuous", "continuous_special")) {
      if (is.numeric(data_final[[var]])) {
        summary_stats[summary_stats$variable == var, "mean"] <- 
          mean(data_final[[var]], na.rm = TRUE)
        summary_stats[summary_stats$variable == var, "sd"] <- 
          sd(data_final[[var]], na.rm = TRUE)
      }
    }
  }
  
  # Return prepared data and metadata
  return(list(
    data = data_final,
    tree = tree_pruned,
    n_species = nrow(data_final),
    variable_info = variable_info,
    summary_stats = summary_stats,
    removed_species = list(
      missing_data = setdiff(data$species, data_complete$species),
      not_in_tree = setdiff(data_complete$species, tree_species),
      not_in_data = setdiff(tree_species, data_complete$species)
    ),
    config = config
  ))
}

#' Create interaction variables if needed
#' 
#' @param data Data frame
#' @param interaction_terms Character vector of interaction terms (e.g., "var1:var2")
#' @return Data frame with interaction columns added
create_interaction_variables <- function(data, interaction_terms) {
  
  for (term in interaction_terms) {
    vars <- strsplit(term, ":")[[1]]
    if (length(vars) == 2 && all(vars %in% colnames(data))) {
      # Create interaction column
      int_name <- paste0(vars[1], "_X_", vars[2])
      
      # Handle different variable types
      if (all(sapply(data[vars], is.numeric))) {
        data[[int_name]] <- data[[vars[1]]] * data[[vars[2]]]
      } else {
        # For factors, create dummy interactions
        warning(paste("Complex interaction between non-numeric variables:", term))
      }
    }
  }
  
  return(data)
}

#' Check data quality for phyloglm
#' 
#' @param prepared_data Output from prepare_analysis_data
#' @return List of data quality metrics and warnings
check_data_quality <- function(prepared_data) {
  
  data <- prepared_data$data
  config <- prepared_data$config
  var_info <- prepared_data$variable_info
  
  issues <- list()
  warnings <- list()
  
  # Check sample size
  n <- prepared_data$n_species
  if (n < 50) {
    issues$sample_size <- paste("Very small sample size:", n)
  } else if (n < 100) {
    warnings$sample_size <- paste("Limited sample size:", n)
  }
  
  # Check response variable
  response_var <- config$response
  response_type <- var_info[[response_var]]$type
  
  if (response_type == "binary") {
    # Check balance
    prop_1 <- mean(data[[response_var]] == 1, na.rm = TRUE)
    if (prop_1 < 0.05 || prop_1 > 0.95) {
      issues$response_balance <- paste("Extremely unbalanced response:", 
                                      round(prop_1 * 100, 1), "% positive")
    } else if (prop_1 < 0.1 || prop_1 > 0.9) {
      warnings$response_balance <- paste("Unbalanced response:", 
                                        round(prop_1 * 100, 1), "% positive")
    }
  }
  
  # Check predictors
  for (pred in config$predictors) {
    pred_type <- var_info[[pred]]$type
    
    if (pred_type == "binary") {
      prop_1 <- mean(data[[pred]] == 1, na.rm = TRUE)
      if (prop_1 < 0.05 || prop_1 > 0.95) {
        issues[[paste0(pred, "_balance")]] <- 
          paste("Extremely unbalanced predictor", pred, ":", 
                round(prop_1 * 100, 1), "% positive")
      }
    } else if (pred_type == "categorical") {
      tab <- table(data[[pred]])
      min_cat <- min(tab)
      if (min_cat < 5) {
        issues[[paste0(pred, "_categories")]] <- 
          paste("Category with <5 observations in", pred)
      } else if (min_cat < 10) {
        warnings[[paste0(pred, "_categories")]] <- 
          paste("Small category in", pred, ":", min_cat, "observations")
      }
    }
  }
  
  # Check for multicollinearity among continuous predictors
  continuous_preds <- config$predictors[sapply(config$predictors, function(v) {
    var_info[[v]]$type %in% c("continuous", "continuous_special")
  })]
  
  if (length(continuous_preds) >= 2) {
    cor_matrix <- cor(data[continuous_preds], use = "complete.obs")
    high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
    
    if (nrow(high_cor) > 0) {
      for (i in 1:nrow(high_cor)) {
        if (high_cor[i, 1] < high_cor[i, 2]) {  # Avoid duplicates
          var1 <- rownames(cor_matrix)[high_cor[i, 1]]
          var2 <- colnames(cor_matrix)[high_cor[i, 2]]
          r <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
          warnings$collinearity <- c(warnings$collinearity,
                                    paste(var1, "and", var2, "are highly correlated: r =", 
                                          round(r, 2)))
        }
      }
    }
  }
  
  # Check tree
  if (!is.ultrametric(prepared_data$tree)) {
    warnings$tree <- "Tree is not ultrametric"
  }
  
  # Create quality report
  quality_report <- list(
    n_species = n,
    n_variables = length(c(config$response, config$predictors, config$controls)),
    issues = issues,
    warnings = warnings,
    is_acceptable = length(issues) == 0
  )
  
  return(quality_report)
}

#' Create a data summary table
#' 
#' @param prepared_data Output from prepare_analysis_data
#' @return Data frame with summary statistics
create_data_summary <- function(prepared_data) {
  
  data <- prepared_data$data
  var_info <- prepared_data$variable_info
  all_vars <- names(var_info)
  
  summaries <- data.frame(
    Variable = character(),
    Type = character(),
    N = integer(),
    Missing = integer(),
    Summary = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in all_vars) {
    var_data <- data[[var]]
    var_type <- var_info[[var]]$type
    
    n_complete <- sum(!is.na(var_data))
    n_missing <- sum(is.na(var_data))
    
    # Create appropriate summary based on type
    if (var_type == "binary") {
      n_1 <- sum(var_data == 1, na.rm = TRUE)
      n_0 <- sum(var_data == 0, na.rm = TRUE)
      summary_str <- paste0("0: ", n_0, " (", round(100 * n_0/n_complete, 1), "%); ",
                           "1: ", n_1, " (", round(100 * n_1/n_complete, 1), "%)")
    } else if (var_type == "categorical") {
      if (is.factor(var_data)) {
        tab <- table(var_data, useNA = "no")
        summary_str <- paste(paste0(names(tab), ": ", tab), collapse = "; ")
      } else {
        summary_str <- paste("Levels:", length(unique(var_data[!is.na(var_data)])))
      }
    } else if (var_type %in% c("continuous", "continuous_special")) {
      summary_str <- paste0("Mean: ", round(mean(var_data, na.rm = TRUE), 2),
                           " (SD: ", round(sd(var_data, na.rm = TRUE), 2), "); ",
                           "Range: [", round(min(var_data, na.rm = TRUE), 2),
                           ", ", round(max(var_data, na.rm = TRUE), 2), "]")
    } else {
      summary_str <- "Unknown type"
    }
    
    summaries <- rbind(summaries, data.frame(
      Variable = var,
      Type = var_type,
      N = n_complete,
      Missing = n_missing,
      Summary = summary_str,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summaries)
}