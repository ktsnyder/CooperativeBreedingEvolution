# Comprehensive Formula Builder for PhyloGLM Framework
# Generates ALL possible model combinations for a given set of predictors

source("phyloglm_framework/variable_classification.R")

#' Build ALL possible phyloglm formulas for a given set of predictors
#' 
#' For 3 predictors (A, B, C), this generates all 15 possible models:
#' 1. Null
#' 2-4. Single predictors: A, B, C
#' 5-7. Pairs additive: A+B, A+C, B+C
#' 8. All additive: A+B+C
#' 9-11. Single interactions: A*B, A*C, B*C
#' 12-14. Interaction + additive: A*B+C, A*C+B, B*C+A
#' 15. Full three-way: A*B*C
#' 
#' @param response_var Response variable name
#' @param predictor_vars Vector of predictor variable names
#' @param control_vars Vector of control variables (always included)
#' @param variable_types Named list of variable types
#' @return Named list of formula objects
build_all_phyloglm_formulas <- function(response_var, 
                                       predictor_vars, 
                                       control_vars = NULL,
                                       variable_types = NULL) {
  
  formulas <- list()
  
  # Helper function to add controls to formula string
  add_controls <- function(pred_str, control_vars) {
    if (!is.null(control_vars) && length(control_vars) > 0) {
      control_str <- paste(control_vars, collapse = " + ")
      return(paste(pred_str, "+", control_str))
    }
    return(pred_str)
  }
  
  # Helper function to create clean model names
  make_model_name <- function(preds, interaction_pairs = NULL, three_way = FALSE) {
    pred_names <- sapply(preds, simplify_var_name)
    
    if (three_way) {
      return(paste0("ThreeWay_", paste(pred_names, collapse = "x")))
    } else if (!is.null(interaction_pairs)) {
      int_names <- paste(sapply(interaction_pairs, simplify_var_name), collapse = "x")
      other_preds <- setdiff(preds, interaction_pairs)
      if (length(other_preds) > 0) {
        other_names <- paste(sapply(other_preds, simplify_var_name), collapse = "_")
        return(paste0("Int_", int_names, "_Add_", other_names))
      } else {
        return(paste0("Int_", int_names))
      }
    } else {
      if (length(preds) == 1) {
        return(paste0("Main_", pred_names))
      } else {
        return(paste0("Add_", paste(pred_names, collapse = "_")))
      }
    }
  }
  
  # 1. NULL MODEL
  if (!is.null(control_vars)) {
    control_str <- paste(control_vars, collapse = " + ")
    formulas[["Null"]] <- as.formula(paste(response_var, "~", control_str))
  } else {
    formulas[["Null"]] <- as.formula(paste(response_var, "~ 1"))
  }
  
  # 2. SINGLE PREDICTOR MODELS
  for (pred in predictor_vars) {
    model_name <- make_model_name(pred)
    pred_str <- add_controls(pred, control_vars)
    formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
  }
  
  # 3. TWO-PREDICTOR MODELS (if we have at least 2 predictors)
  if (length(predictor_vars) >= 2) {
    pred_pairs <- combn(predictor_vars, 2, simplify = FALSE)
    
    for (pair in pred_pairs) {
      # Additive model
      model_name <- make_model_name(pair)
      pred_str <- paste(pair, collapse = " + ")
      pred_str <- add_controls(pred_str, control_vars)
      formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
      
      # Interaction model
      model_name <- make_model_name(pair, interaction_pairs = pair)
      pred_str <- paste(pair[1], "*", pair[2])
      pred_str <- add_controls(pred_str, control_vars)
      formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
    }
  }
  
  # 4. THREE-PREDICTOR MODELS (if we have at least 3 predictors)
  if (length(predictor_vars) >= 3) {
    # All three additive
    model_name <- make_model_name(predictor_vars)
    pred_str <- paste(predictor_vars, collapse = " + ")
    pred_str <- add_controls(pred_str, control_vars)
    formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
    
    # Each two-way interaction plus the third as additive
    pred_pairs <- combn(predictor_vars, 2, simplify = FALSE)
    
    for (pair in pred_pairs) {
      remaining <- setdiff(predictor_vars, pair)
      model_name <- make_model_name(predictor_vars, interaction_pairs = pair)
      pred_str <- paste(paste(pair[1], "*", pair[2]), "+", paste(remaining, collapse = " + "))
      pred_str <- add_controls(pred_str, control_vars)
      formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
    }
    
    # Three-way interaction (if exactly 3 predictors)
    if (length(predictor_vars) == 3) {
      model_name <- make_model_name(predictor_vars, three_way = TRUE)
      pred_str <- paste(predictor_vars[1], "*", predictor_vars[2], "*", predictor_vars[3])
      pred_str <- add_controls(pred_str, control_vars)
      formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
    }
  }
  
  # 5. FOUR OR MORE PREDICTORS (generate systematically)
  if (length(predictor_vars) > 3) {
    # All combinations of size 3
    if (length(predictor_vars) >= 3) {
      pred_triples <- combn(predictor_vars, 3, simplify = FALSE)
      for (triple in pred_triples) {
        # Additive model
        model_name <- make_model_name(triple)
        pred_str <- paste(triple, collapse = " + ")
        pred_str <- add_controls(pred_str, control_vars)
        formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
      }
    }
    
    # All predictors additive
    model_name <- paste0("Add_All", length(predictor_vars), "Predictors")
    pred_str <- paste(predictor_vars, collapse = " + ")
    pred_str <- add_controls(pred_str, control_vars)
    formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
  }
  
  return(formulas)
}

#' Build model set from configuration using comprehensive formula generation
#' 
#' @param config Analysis configuration
#' @param data Data frame to check variable types
#' @return Named list of formulas
build_comprehensive_model_set <- function(config, data) {
  
  # Get variable types
  classifications <- get_variable_classifications()
  all_vars <- c(config$response, config$predictors, config$controls)
  
  variable_types <- list()
  for (var in all_vars) {
    if (var %in% colnames(data)) {
      variable_types[[var]] <- detect_variable_type(data, var, classifications)
    }
  }
  
  # Generate all formulas
  formulas <- build_all_phyloglm_formulas(
    response_var = config$response,
    predictor_vars = config$predictors,
    control_vars = config$controls,
    variable_types = variable_types
  )
  
  return(formulas)
}

#' Simplify variable name for model naming
#' 
#' @param var_name Full variable name
#' @return Simplified name
simplify_var_name <- function(var_name) {
  # Remove common suffixes and special characters
  simple <- gsub("_Agg01|_AVONET|_Jetz|2017", "", var_name)
  simple <- gsub("HighConfidence_", "", simple)
  simple <- gsub("TerritorialityWeakVsStrong", "TerrWS", simple)
  simple <- gsub("Territory", "Terr", simple)
  simple <- gsub("logMass", "Mass", simple)
  simple <- gsub("PercentAbsLog", "", simple)
  simple <- gsub("Dimorphism", "Dim", simple)
  simple <- gsub("logMaleFemalePlumage", "Plum", simple)
  simple <- gsub("DiffAbs", "", simple)
  simple <- gsub("GeographicRegion", "Region", simple)
  simple <- gsub("Migration", "Migr", simple)
  simple <- gsub("Centroid.Latitude", "Lat", simple)
  simple <- gsub("FemaleSong", "FS", simple)
  simple <- gsub("Griesser", "", simple)
  simple <- gsub("FamilialLiving", "FL", simple)
  simple <- gsub("Coop", "CB", simple)
  
  # Remove any remaining special characters
  simple <- gsub("[^A-Za-z0-9]", "", simple)
  
  return(simple)
}