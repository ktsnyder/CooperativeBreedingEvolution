# Automated Formula Building System for PhyloGLM Framework
# Handles different variable types and creates appropriate model sets
# 8/26/2025 - KTS commented out auto-building the 15 standard formulas whenever 2 predictors plus 1 control are present, so it should build just based on the complexity_levels now in build_model_set_from_config.

source("phyloglm_framework/variable_classification.R")

#' Build a comprehensive set of phyloglm formulas
#' 
#' @param response_var Character string of response variable name
#' @param predictor_vars Character vector of predictor variable names
#' @param control_vars Character vector of control variables (always included)
#' @param variable_types Named list of variable types (from detect_variable_type)
#' @param complexity_levels Character vector specifying which complexity levels to include
#' @param max_interactions Maximum order of interactions to consider
#' @return Named list of formula objects
build_phyloglm_formulas <- function(response_var, 
                                   predictor_vars, 
                                   control_vars = NULL,
                                   variable_types,
                                   complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
                                   max_interactions = 2) {
  
  formulas <- list()
  
  # Ensure response is not in predictors
  predictor_vars <- setdiff(predictor_vars, response_var)
  
  # 1. NULL MODEL
  if ("null" %in% complexity_levels) {
    formulas[["Null"]] <- as.formula(paste(response_var, "~ 1"))
  }
  
  # 2. MAIN EFFECTS MODELS (single predictors)
  if ("main" %in% complexity_levels) {
    for (pred in predictor_vars) {
      model_name <- paste0("MainPred_", simplify_var_name(pred))
      if (is.null(control_vars)) {
        formulas[[model_name]] <- as.formula(paste(response_var, "~", pred))
      } else {
        control_str <- paste(control_vars, collapse = " + ")
        formulas[[model_name]] <- as.formula(paste(response_var, "~", pred, "+", control_str))
      }
    }
    
    # Control variables only models
    if (!is.null(control_vars)) {
      for (ctrl in control_vars) {
        model_name <- paste0("Control_", simplify_var_name(ctrl))
        formulas[[model_name]] <- as.formula(paste(response_var, "~", ctrl))
      }
    }
  }
  
  # 3. ADDITIVE MODELS (multiple predictors, no interactions)
  if ("additive" %in% complexity_levels && length(predictor_vars) >= 2) {
    # Two predictor combinations
    if (length(predictor_vars) >= 2) {
      pred_combos <- combn(predictor_vars, 2, simplify = FALSE)
      for (combo in pred_combos) {
        model_name <- paste0("Add_", paste(sapply(combo, simplify_var_name), collapse = "_"))
        pred_str <- paste(combo, collapse = " + ")
        if (!is.null(control_vars)) {
          control_str <- paste(control_vars, collapse = " + ")
          formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str, "+", control_str))
        } else {
          formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
        }
      }
    }
    
    # All predictors together
    if (length(predictor_vars) > 2) {
      model_name <- "Add_AllPredictors"
      pred_str <- paste(predictor_vars, collapse = " + ")
      if (!is.null(control_vars)) {
        control_str <- paste(control_vars, collapse = " + ")
        formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str, "+", control_str))
      } else {
        formulas[[model_name]] <- as.formula(paste(response_var, "~", pred_str))
      }
    }
  }
  
  # 4. TWO-WAY INTERACTION MODELS
  if ("twoway" %in% complexity_levels && length(predictor_vars) >= 2 && max_interactions >= 2) {
    pred_combos <- combn(predictor_vars, 2, simplify = FALSE)
    
    for (combo in pred_combos) {
      var1 <- combo[1]
      var2 <- combo[2]
      type1 <- variable_types[[var1]]
      type2 <- variable_types[[var2]]
      
      # Full interaction model
      model_name <- paste0("Int_", simplify_var_name(var1), "x", simplify_var_name(var2))
      if (!is.null(control_vars)) {
        control_str <- paste(control_vars, collapse = " + ")
        formulas[[model_name]] <- as.formula(
          paste(response_var, "~", var1, "*", var2, "+", control_str)
        )
      } else {
        formulas[[model_name]] <- as.formula(paste(response_var, "~", var1, "*", var2))
      }
      
      # Additional models if there are control variables
      if (!is.null(control_vars) && length(control_vars) > 0) {
        # Interaction with first predictor and control
        for (ctrl in control_vars) {
          model_name <- paste0("Int_", simplify_var_name(var1), "x", simplify_var_name(ctrl))
          other_vars <- c(var2, setdiff(control_vars, ctrl))
          other_str <- paste(other_vars, collapse = " + ")
          formulas[[model_name]] <- as.formula(
            paste(response_var, "~", var1, "*", ctrl, "+", other_str)
          )
        }
      }
    }
  }
  
  # 5. THREE-WAY INTERACTION MODELS
  if ("threeway" %in% complexity_levels && length(predictor_vars) >= 3 && max_interactions >= 3) {
    # Only include if we have exactly 2 predictors + 1 control
    if (length(predictor_vars) == 2 && length(control_vars) == 1) {
      model_name <- "FullInteractions3"
      formulas[[model_name]] <- as.formula(
        paste(response_var, "~", predictor_vars[1], "*", predictor_vars[2], "*", control_vars[1])
      )
    }
    # Or 3 predictors without controls
    else if (length(predictor_vars) == 3 && length(control_vars) == 0) {
      model_name <- "FullInteractions3"
      formulas[[model_name]] <- as.formula(
        paste(response_var, "~", predictor_vars[1], "*", predictor_vars[2], "*", predictor_vars[3])
      )
    }
  }
  
  # 6. FOUR-WAY INTERACTION MODELS
  if ("fourway" %in% complexity_levels && length(c(predictor_vars, control_vars)) >= 4 && max_interactions >= 4) {
    # Include if we have 3 predictors + 1 control
    if (length(predictor_vars) == 3 && length(control_vars) == 1) {
      model_name <- "FullInteractions4"
      formulas[[model_name]] <- as.formula(
        paste(response_var, "~", predictor_vars[1], "*", predictor_vars[2], "*", predictor_vars[3], "*", control_vars[1])
      )
    }
  }
  
  # 7. SPECIAL MODELS for specific variable type combinations
  formulas <- add_special_models(formulas, response_var, predictor_vars, 
                                control_vars, variable_types)
  
  return(formulas)
}

#' Add special models based on variable type combinations
#' 
#' @param formulas Existing list of formulas
#' @param response_var Response variable name
#' @param predictor_vars Predictor variable names
#' @param control_vars Control variable names
#' @param variable_types Variable types
#' @return Updated formula list
add_special_models <- function(formulas, response_var, predictor_vars, 
                              control_vars, variable_types) {
  
  # Special handling for categorical predictors with >2 levels
  categorical_preds <- predictor_vars[sapply(predictor_vars, function(x) {
    variable_types[[x]] == "categorical" && x %in% c("Territory", "Migration_AVONET")
  })]
  
  if (length(categorical_preds) > 0) {
    for (cat_pred in categorical_preds) {
      # Add polynomial contrasts model if sensible
      if (cat_pred %in% c("Territory", "Migration_AVONET")) {
        model_name <- paste0("Poly_", simplify_var_name(cat_pred))
        # This will be handled in the model fitting with contrast specification
        formulas[[model_name]] <- as.formula(paste(response_var, "~", cat_pred))
      }
    }
  }
  
  # Special handling for continuous predictors that might need polynomials
  continuous_preds <- predictor_vars[sapply(predictor_vars, function(x) {
    variable_types[[x]] %in% c("continuous", "continuous_special")
  })]
  
  latitude_vars <- grep("Latitude", continuous_preds, value = TRUE)
  if (length(latitude_vars) > 0) {
    # Add absolute latitude model
    for (lat_var in latitude_vars) {
      model_name <- paste0("Abs_", simplify_var_name(lat_var))
      # Note: actual transformation happens in data prep
      if (!is.null(control_vars)) {
        control_str <- paste(control_vars, collapse = " + ")
        formulas[[model_name]] <- as.formula(
          paste(response_var, "~ abs(", lat_var, ") +", control_str)
        )
      } else {
        formulas[[model_name]] <- as.formula(paste(response_var, "~ abs(", lat_var, ")"))
      }
    }
  }
  
  return(formulas)
}

#' Simplify variable name for model naming
#' 
#' @param var_name Full variable name
#' @return Simplified name for use in model labels
simplify_var_name <- function(var_name) {
  simplifications <- list(
    "HighConfidence_Coop" = "CB",
    "FemaleSong_Agg01" = "FS", 
    "Griesser2017FamilialLiving" = "Fam",
    "TerritorialityWeakVsStrong" = "TerrWS",
    "Territory_12vs3" = "Terr123",
    "Territory" = "Terr3",
    "Migration_AVONET" = "Migr",
    "logMass_AVONET" = "Mass",
    "logMass_normalized" = "Mass",
    "PercentAbsLogWingDimorphism" = "WingDim",
    "logMaleFemalePlumageDiffAbs" = "PlumDim",
    "Centroid.Latitude_AVONET" = "Lat",
    "GeographicRegion_Jetz" = "GeoReg"
  )
  
  if (var_name %in% names(simplifications)) {
    return(simplifications[[var_name]])
  } else {
    # Generic simplification
    return(substr(gsub("[^A-Za-z0-9]", "", var_name), 1, 8))
  }
}

#' Generate model names from formulas
#' 
#' @param formula_list List of formula objects
#' @param response_var Response variable name
#' @param predictor_vars Predictor variable names
#' @return Named list with standardized model names
standardize_model_names <- function(formula_list, response_var, predictor_vars) {
  new_names <- names(formula_list)
  
  # Additional standardization if needed
  for (i in seq_along(formula_list)) {
    f <- formula_list[[i]]
    terms <- attr(terms(f), "term.labels")
    
    # Detect interaction terms
    has_interaction <- any(grepl(":", terms))
    
    if (has_interaction) {
      # Extract main effects from interaction
      main_effects <- unique(unlist(strsplit(terms[grepl(":", terms)], ":")))
      main_effects <- intersect(main_effects, predictor_vars)
      
      if (length(main_effects) == 2) {
        new_names[i] <- paste0(simplify_var_name(main_effects[1]), "x",
                              simplify_var_name(main_effects[2]))
      }
    }
  }
  
  names(formula_list) <- new_names
  return(formula_list)
}

#' Build model set for a specific analysis configuration
#' 
#' @param config List containing analysis configuration
#' @param data Data frame for variable type detection
#' @return List of formulas
build_model_set_from_config <- function(config, data) {
  
  # Check if we should use the comprehensive 15-model set
  # This applies when we have 2 predictors and 1 control (standard analysis)
  # if (length(config$predictors) == 2 && length(config$controls) == 1 && 
  #     config$controls[1] %in% c("logMass_AVONET", "logMass_normalized", 
  #                               "PercentAbsLogWingDimorphism", 
  #                               "logMaleFemalePlumageDiffAbs")) {
  #   
  #   # Use the comprehensive 15-model formula builder
  #   formulas <- build_comprehensive_15_models(
  #     response_var = config$response,
  #     pred_var = config$predictors[1],
  #     terr_var = config$predictors[2],
  #     control_var = config$controls[1]
  #   )
  #   
  #   return(formulas)
  # }
  
  # Otherwise use the standard formula builder
  # Detect variable types
  all_vars <- unique(c(config$response, config$predictors, config$controls))
  variable_types <- list()
  
  classifications <- get_variable_classifications()
  for (var in all_vars) {
    variable_types[[var]] <- detect_variable_type(data, var, classifications)
  }
  
  # Build formulas
  formulas <- build_phyloglm_formulas(
    response_var = config$response,
    predictor_vars = config$predictors,
    control_vars = config$controls,
    variable_types = variable_types,
    complexity_levels = config$complexity_levels,
    max_interactions = config$max_interactions
  )
  
  return(formulas)
}

#' Build comprehensive 15-model set following standard naming scheme
#' 
#' @param response_var Response variable name
#' @param pred_var Main predictor variable
#' @param terr_var Secondary predictor (e.g., territoriality)
#' @param control_var Control variable (e.g., mass)
#' @return List of 15 formulas with standard names
build_comprehensive_15_models <- function(response_var, pred_var, terr_var, control_var) {
  
  # Create short names for cleaner model names
  pred_short <- simplify_var_name(pred_var)
  terr_short <- simplify_var_name(terr_var)
  control_short <- simplify_var_name(control_var)
  
  # Build all 15 models with exact naming scheme
  formulas <- list(
    # 1. Null model
    Null = as.formula(paste(response_var, "~ 1")),
    
    # 2-4. Single predictor models
    MainPred = as.formula(paste(response_var, "~", pred_var)),
    Terr = as.formula(paste(response_var, "~", terr_var)),
    Mass = as.formula(paste(response_var, "~", control_var)),
    
    # 5-7. Two predictor additive models
    MainPred_Terr = as.formula(paste(response_var, "~", pred_var, "+", terr_var)),
    MainPred_Mass = as.formula(paste(response_var, "~", pred_var, "+", control_var)),
    Terr_Mass = as.formula(paste(response_var, "~", terr_var, "+", control_var)),
    
    # 8-10. Two predictor interaction models
    MainPredxTerr = as.formula(paste(response_var, "~", pred_var, "*", terr_var)),
    MainPredxMass = as.formula(paste(response_var, "~", pred_var, "*", control_var)),
    TerrxMass = as.formula(paste(response_var, "~", terr_var, "*", control_var)),
    
    # 11. Three predictor additive model
    MainPred_Terr_Mass = as.formula(paste(response_var, "~", pred_var, "+", terr_var, "+", control_var)),
    
    # 12-14. One interaction plus one additive
    MainPredxTerr_Mass = as.formula(paste(response_var, "~", pred_var, "*", terr_var, "+", control_var)),
    MainPred_TerrxMass = as.formula(paste(response_var, "~", pred_var, "+", terr_var, "*", control_var)),
    Terr_MainPredxMass = as.formula(paste(response_var, "~", terr_var, "+", pred_var, "*", control_var)),
    
    # 15. Full three-way interaction
    FullInteractions = as.formula(paste(response_var, "~", pred_var, "*", terr_var, "*", control_var))
  )
  
  return(formulas)
}

#' Create a model matrix showing which variables are in each model
#' 
#' @param formula_list List of formulas
#' @return Data frame showing model structures
create_model_matrix <- function(formula_list) {
  
  # Extract all unique variables
  all_vars <- unique(unlist(lapply(formula_list, function(f) {
    all.vars(f)[-1]  # Exclude response variable
  })))
  
  # Create matrix
  model_matrix <- matrix(0, nrow = length(formula_list), ncol = length(all_vars))
  rownames(model_matrix) <- names(formula_list)
  colnames(model_matrix) <- all_vars
  
  for (i in seq_along(formula_list)) {
    f <- formula_list[[i]]
    vars_in_model <- all.vars(f)[-1]
    terms_in_model <- attr(terms(f), "term.labels")
    
    # Mark main effects
    for (var in vars_in_model) {
      model_matrix[i, var] <- 1
    }
    
    # Mark interactions
    interaction_terms <- terms_in_model[grepl(":", terms_in_model)]
    for (int_term in interaction_terms) {
      vars_in_int <- strsplit(int_term, ":")[[1]]
      # Mark as 2 for interaction
      for (var in vars_in_int) {
        if (var %in% colnames(model_matrix)) {
          model_matrix[i, var] <- 2
        }
      }
    }
  }
  
  # Convert to data frame and add model names
  model_df <- as.data.frame(model_matrix)
  model_df$Model <- rownames(model_matrix)
  model_df <- model_df[, c("Model", all_vars)]
  
  return(model_df)
}