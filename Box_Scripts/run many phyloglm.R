#################################################
# Comprehensive Analysis of Female Song and Cooperative Breeding using Phyloglm
# Testing multiple territoriality variables and bidirectional relationships
# Date: May 14, 2025
#################################################

# Purpose:
# This script performs a comprehensive comparative analysis of the evolutionary 
# relationships between female song, cooperative breeding, and multiple measures 
# of territoriality in birds, using phylogenetic logistic regression.
#
# The script conducts bidirectional analyses (each variable as both predictor and 
# response) across six different territoriality variables, systematically testing 
# a range of models with different interaction structures, and producing comprehensive 
# visualizations and summaries.
#
# Main Functionalities:
# 1. Data preparation and cleaning specific to each analysis
# 2. Model fitting (including null, simple, and interaction models)
# 3. Model selection using AIC
# 4. Rerunning top models with bootstrapping for robust confidence intervals
# 5. Visualization of model predictions and interactions
# 6. Coefficient extraction and significance testing
# 7. Model averaging to account for model uncertainty
# 8. Effect size calculation (odds ratios, marginal effects, etc.)
# 9. Generation of comprehensive summary tables
#
# Input Requirements:
# - A dataset (dfIn) containing binary and continuous variables:
#   * FemaleSong_Agg01: Binary response variable for female song presence (0=absent, 1=present)
#   * HighConfidence_Coop: Binary predictor for cooperative breeding (0=absent, 1=present)
#   * Multiple territoriality variables (TerritorialityWeakVsStrong, etc.)
#   * Mass_AVONET: Body mass to be log-transformed
# - A phylogenetic tree (tree) in 'phylo' format
#
# Output Files:
# - all_model_comparisons.csv: AIC comparison tables for all model sets
# - all_coefficient_summaries.csv: Detailed parameter estimates for top models
# - best_models_summary.csv: Summary of the single best model for each analysis
# - all_model_averages.csv: Model-averaged parameter estimates
# - all_interaction_plots.pdf: Visualizations of effects from top models
# - effect_sizes.csv: Odds ratios and confidence intervals
# - marginal_effects.csv: Average changes in probability
# - standardized_coefficients.csv: Standardized effect sizes
# - predictor_effects.csv & territoriality_effects.csv: Effect sizes in probability scale
#
# Key Analyses:
# - Tests for interactions between cooperation and territoriality on female song
# - Tests for interactions between territoriality and body mass
# - Tests for interactions between cooperation and body mass
# - Examines bidirectional relationships (each variable as both predictor and response)
# - Compares effects across different measures of territoriality
# - Quantifies effect sizes in multiple ways (odds ratios, probability changes)
#
# Best practice (Usage):
# 1. First run the entire script to generate all outputs
# 2. Review model_comparisons.csv to understand model fit
# 3. Examine coefficient_summaries.csv for parameter significance
# 4. Use interaction_plots.pdf to visualize key relationships
# 5. Consult effect_sizes.csv to interpret biological significance
#
# Dependencies:
# - phylolm: For phylogenetic logistic regression
# - ggplot2: For visualization
# - dplyr & tidyr: For data manipulation
#
#################################################

# Load required packages
library(phylolm)    # For phylogenetic logistic regression
library(ggplot2)    # For visualization
library(gridExtra)  # For arranging multiple plots
library(dplyr)      # For data manipulation
library(tidyr)      # For reshaping data

###
#### Function to prepare data for specific variables ----
###

prepare_cleaned_data <- function(data, tree, response_var, terr_var) {
  # Identify required variables
  required_vars <- c(response_var, terr_var, "HighConfidence_Coop", "FemaleSong_Agg01", "logMass_AVONET")
  
  # Remove duplicate species (just in case)
  data <- data[!duplicated(data$species), ]
  
  # Filter to complete cases for the required variables
  data_clean <- data[complete.cases(data[, required_vars]), ]
  
  # Drop tips from the tree that aren't in the cleaned data
  tree_clean <- drop.tip(tree, setdiff(tree$tip.label, data_clean$species))
  
  # Ensure rownames match
  rownames(data_clean) <- data_clean$species
  
  # Return both cleaned data and tree
  return(list(data = data_clean, tree = tree_clean))
}

###
#### Function to run all models for a specific response and territoriality variable ----
###

run_model_set <- function(response_var, terr_var, original_data, original_tree) {
  # First, prepare clean data specific to these variables
  cleaned <- prepare_cleaned_data(original_data, original_tree, response_var, terr_var)
  data <- cleaned$data
  tree <- cleaned$tree
  
  # Report on data dimensions
  message(paste("Analysis for", response_var, "~", terr_var, "using", nrow(data), "species"))
  
  # Define the other main predictor based on response
  other_pred <- ifelse(response_var == "FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
  
  # Create formula strings by substituting the specific variables
  formulas <- list(
    f0 = paste(response_var, "~ 1"),
    f1 = paste(response_var, "~", other_pred),
    f2 = paste(response_var, "~", terr_var),
    f2b = paste(response_var, "~ logMass_AVONET"),
    f3ni = paste(response_var, "~", other_pred, "+", terr_var),
    f3i = paste(response_var, "~", other_pred, "*", terr_var),
    f1wMassNI = paste(response_var, "~", other_pred, "+ logMass_AVONET"),
    f1wMassI = paste(response_var, "~", other_pred, "* logMass_AVONET"),
    f2wMassNI = paste(response_var, "~", terr_var, "+ logMass_AVONET"),
    f2wMassI = paste(response_var, "~", terr_var, "* logMass_AVONET"),
    f3niWmassNI = paste(response_var, "~", other_pred, "+", terr_var, "+ logMass_AVONET"),
    f3iWmassI = paste(response_var, "~", other_pred, "*", terr_var, "* logMass_AVONET"),
    fInteract1 = paste(response_var, "~", other_pred, "*", terr_var, "+ logMass_AVONET"),
    fInteract2 = paste(response_var, "~", other_pred, "+", terr_var, "* logMass_AVONET"),
    fInteract3 = paste(response_var, "~", terr_var, "+", other_pred, "* logMass_AVONET")
  )
  
  # Run all models
  models <- list()
  for (i in seq_along(formulas)) {
    model_name <- names(formulas)[i]
    tryCatch({
      models[[model_name]] <- phyloglm(
        formula = as.formula(formulas[[i]]),
        data = data,
        phy = tree,
        method = c("logistic_MPLE"),
        btol = 20,  # Increased to handle boundary issues
        log.alpha.bound = 4,
        start.beta = NULL,
        start.alpha = NULL,
        boot = 0,
        full.matrix = TRUE,
        save = FALSE
      )
    }, error = function(e) {
      message(paste("Error in model", model_name, ":", e$message))
      models[[model_name]] <- NULL
    })
  }
  
  # Calculate AIC for each model
  aic_values <- sapply(models, AIC)
  
  # Create model comparison table
  model_comparison <- data.frame(
    Model = names(models),
    Formula = unlist(formulas)[names(models)],
    AIC = aic_values,
    n_species = nrow(data)  # Add sample size information
  )
  
  # Sort by AIC and calculate delta AIC and weights
  model_comparison <- model_comparison[order(model_comparison$AIC), ]
  model_comparison$deltaAIC <- model_comparison$AIC - min(model_comparison$AIC)
  model_comparison$weight <- exp(-0.5 * model_comparison$deltaAIC) / 
    sum(exp(-0.5 * model_comparison$deltaAIC))
  model_comparison$cum_weight <- cumsum(model_comparison$weight)
  
  # Round numeric columns for display
  model_comparison$AIC <- round(model_comparison$AIC, 2)
  model_comparison$deltaAIC <- round(model_comparison$deltaAIC, 2)
  model_comparison$weight <- round(model_comparison$weight, 3)
  model_comparison$cum_weight <- round(model_comparison$cum_weight, 3)
  
  # Return both models and comparison, along with the cleaned data
  return(list(
    models = models,
    comparison = model_comparison,
    response = response_var,
    territoriality = terr_var,
    data = data,
    tree = tree
  ))
}

###
#### Function to identify top models and run with bootstrapping ----
###

run_top_models_with_boot <- function(model_results, n_boot = 100) {
  # Get models with deltaAIC < 2
  top_models <- model_results$comparison$Model[model_results$comparison$deltaAIC < 2]
  
  # If no models with deltaAIC < 2, just take the best model
  if (length(top_models) == 0) {
    top_models <- model_results$comparison$Model[1]
  }
  
  # Get the cleaned data and tree from the model_results
  data <- model_results$data
  tree <- model_results$tree
  
  # Re-run top models with bootstrapping
  top_model_results <- list()
  
  for (model_name in top_models) {
    formula_str <- model_results$comparison$Formula[model_results$comparison$Model == model_name]
    
    top_model_results[[model_name]] <- tryCatch({
      phyloglm(
        formula = as.formula(formula_str),
        data = data,
        phy = tree,
        method = c("logistic_MPLE"),
        btol = 20,  # Increased to handle boundary issues
        log.alpha.bound = 4,
        start.beta = NULL,
        start.alpha = NULL,
        boot = n_boot,
        full.matrix = TRUE,
        save = TRUE
      )
    }, error = function(e) {
      message(paste("Error in bootstrap for model", model_name, ":", e$message))
      # Return the original model without bootstrapping
      model_results$models[[model_name]]
    })
  }
  
  return(top_model_results)
}

###
#### Function to create effect plots for top models ----
###

create_effect_plot <- function(model, data, response_var, terr_var) {
  tryCatch({
    # Add debugging info
    message(paste("\nCreating plot for model with formula:", deparse(formula(model))))
    message(paste("Response variable:", response_var))
    message(paste("Territoriality variable:", terr_var))
    
    # Determine the predictor variables based on the response
    pred_var <- ifelse(response_var == "FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
    message(paste("Predictor variable:", pred_var))
    
    # Extract model formula
    model_formula <- formula(model)
    formula_text <- paste(deparse(model_formula), collapse = " ")
    message(paste("Full formula:", formula_text))
    
    # Create the plot subtitle with formula and AIC
    formula_text <- deparse(formula(model))
    if (length(formula_text) > 1) {
      formula_text <- paste(formula_text, collapse = " ")
    }
    subtitle <- paste(formula_text, "| AIC:", round(AIC(model), 2))
    
    # Check what kind of plot to make
    has_interaction <- grepl(":", formula_text, fixed = TRUE) || grepl("*", formula_text, fixed = TRUE)
    has_pred_var <- grepl(pred_var, formula_text, fixed = TRUE)
    has_terr_var <- grepl(terr_var, formula_text, fixed = TRUE)
    
    message(paste("Has interaction:", has_interaction))
    message(paste("Has predictor:", has_pred_var))
    message(paste("Has territoriality:", has_terr_var))
    
    # Skip if neither predictor is in the model
    if (!has_pred_var && !has_terr_var) {
      message("Neither predictor is in the model, skipping plot")
      return(NULL)
    }
    
    # Start with mean mass value
    mean_mass <- mean(data$logMass_AVONET, na.rm = TRUE)
    message(paste("Mean logMass_AVONET:", mean_mass))
    
    # COMPLETELY DIFFERENT APPROACH: Create the dataframe piece by piece, 
    # then add variables one by one to avoid dimension issues
    
    # Start with a simple data frame with just one row
    new_data <- data.frame(intercept = 1)
    
    # Add mass first - always just the mean value
    new_data$logMass_AVONET <- mean_mass
    
    # Now add the variables for plotting
    if (has_pred_var && has_terr_var) {
      # Safest approach: Build up step by step
      message("Creating grid for both predictor and territoriality")
      
      # Create with explicit combinations
      new_data <- data.frame(
        intercept = rep(1, 4),
        logMass_AVONET = rep(mean_mass, 4)
      )
      
      # Now we'll manually add each combination
      new_data[[pred_var]] <- c(0, 0, 1, 1)
      new_data[[terr_var]] <- c(0, 1, 0, 1)
      
    } else if (has_pred_var) {
      # Only predictor variable
      message("Creating grid for predictor only")
      
      # Create with explicit rows
      new_data <- data.frame(
        intercept = c(1, 1),
        logMass_AVONET = c(mean_mass, mean_mass)
      )
      
      # Add predictor values
      new_data[[pred_var]] <- c(0, 1)
      
      # Add constant territoriality
      new_data[[terr_var]] <- 0
      
    } else if (has_terr_var) {
      # Only territoriality variable
      message("Creating grid for territoriality only")
      
      # Create with explicit rows
      new_data <- data.frame(
        intercept = c(1, 1),
        logMass_AVONET = c(mean_mass, mean_mass)
      )
      
      # Add territoriality values
      new_data[[terr_var]] <- c(0, 1)
      
      # Add constant predictor
      new_data[[pred_var]] <- 0
    }
    
    # Check and report on the dataframe dimensions
    message(paste("Prediction grid has", nrow(new_data), "rows and", ncol(new_data), "columns"))
    message("Column names: ", paste(colnames(new_data), collapse = ", "))
    
    # Extract coefficients from the model
    coefs <- coef(model)
    message("Model coefficients: ", paste(names(coefs), collapse = ", "))
    
    # Calculate the linear predictor manually
    new_data$linear_pred <- coefs["(Intercept)"]
    
    # Add effect of each term if present in the model
    if (pred_var %in% names(coefs) && has_pred_var) {
      new_data$linear_pred <- new_data$linear_pred + coefs[pred_var] * new_data[[pred_var]]
    }
    
    if (terr_var %in% names(coefs) && has_terr_var) {
      new_data$linear_pred <- new_data$linear_pred + coefs[terr_var] * new_data[[terr_var]]
    }
    
    if ("logMass_AVONET" %in% names(coefs)) {
      new_data$linear_pred <- new_data$linear_pred + coefs["logMass_AVONET"] * new_data$logMass_AVONET
    }
    
    # Add interaction terms if present
    interaction_term <- paste0(pred_var, ":", terr_var)
    if (interaction_term %in% names(coefs) && has_pred_var && has_terr_var) {
      new_data$linear_pred <- new_data$linear_pred + 
        coefs[interaction_term] * new_data[[pred_var]] * new_data[[terr_var]]
    }
    
    # Convert to probabilities
    new_data$pred_prob <- plogis(new_data$linear_pred)
    
    # We'll start with simple CIs equal to prediction
    new_data$lower_ci <- new_data$pred_prob
    new_data$upper_ci <- new_data$pred_prob
    
    # Only try bootstrap CIs if the model has them
    if (!is.null(model$bootstrap) && is.matrix(model$bootstrap)) {
      message("Model has bootstrap results, attempting to calculate CIs")
      
      tryCatch({
        # Initialize matrices to store bootstrap predictions
        n_boot <- nrow(model$bootstrap)
        boot_preds <- matrix(NA, nrow = nrow(new_data), ncol = n_boot)
        
        # For each bootstrap sample
        for (i in 1:n_boot) {
          # Extract coefficients from this bootstrap sample
          boot_coefs <- model$bootstrap[i, ]
          
          # Calculate linear predictor for this sample
          boot_linear_pred <- boot_coefs[1]  # Intercept
          
          # Add main effects if present in the model
          if (pred_var %in% names(coefs) && has_pred_var) {
            col_idx <- which(names(coefs) == pred_var)
            if (length(col_idx) > 0 && col_idx <= length(boot_coefs)) {
              boot_linear_pred <- boot_linear_pred + boot_coefs[col_idx] * new_data[[pred_var]]
            }
          }
          
          if (terr_var %in% names(coefs) && has_terr_var) {
            col_idx <- which(names(coefs) == terr_var)
            if (length(col_idx) > 0 && col_idx <= length(boot_coefs)) {
              boot_linear_pred <- boot_linear_pred + boot_coefs[col_idx] * new_data[[terr_var]]
            }
          }
          
          if ("logMass_AVONET" %in% names(coefs)) {
            col_idx <- which(names(coefs) == "logMass_AVONET")
            if (length(col_idx) > 0 && col_idx <= length(boot_coefs)) {
              boot_linear_pred <- boot_linear_pred + boot_coefs[col_idx] * new_data$logMass_AVONET
            }
          }
          
          # Add interaction term if present
          if (interaction_term %in% names(coefs) && has_pred_var && has_terr_var) {
            col_idx <- which(names(coefs) == interaction_term)
            if (length(col_idx) > 0 && col_idx <= length(boot_coefs)) {
              boot_linear_pred <- boot_linear_pred + 
                boot_coefs[col_idx] * new_data[[pred_var]] * new_data[[terr_var]]
            }
          }
          
          # Convert to probabilities
          boot_preds[, i] <- plogis(boot_linear_pred)
        }
        
        # Calculate quantiles for confidence intervals
        new_data$lower_ci <- apply(boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE)
        new_data$upper_ci <- apply(boot_preds, 1, quantile, probs = 0.975, na.rm = TRUE)
        
        message("Successfully calculated bootstrap CIs")
      }, error = function(e) {
        message(paste("Failed to calculate bootstrap CIs:", e$message))
        message("Using point estimates only for plot")
      })
    } else {
      message("Model does not have bootstrap results, using point estimates only")
    }
    
    # Create more readable variable names for the plot
    pred_label <- ifelse(pred_var == "HighConfidence_Coop", "Cooperative Breeding", "Female Song")
    terr_label <- gsub("Territoriality", "Territoriality: ", terr_var)
    terr_label <- gsub("WeakVsStrong", "Weak vs. Strong", terr_label)
    terr_label <- gsub("PermissiveExclusive", "Permissive vs. Exclusive", terr_label)
    terr_label <- gsub("HighConf", "High Confidence", terr_label)
    terr_label <- gsub("Territory_12vs3", "Territory: 1-2 vs 3", terr_label)
    
    # Determine plot title prefix based on whether there's an interaction
    title_prefix <- ifelse(has_interaction, "Interaction between", "Effects of")
    
    # Determine what type of plot to create based on variables in the model
    if (has_pred_var && has_terr_var) {
      message("Creating interaction plot with both variables")
      # Create an interaction/grouped plot
      p <- ggplot(new_data, aes(x = factor(.data[[pred_var]]), y = pred_prob, 
                           group = factor(.data[[terr_var]]), 
                           color = factor(.data[[terr_var]]),
                           fill = factor(.data[[terr_var]]))) +
        # Add confidence interval ribbons
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
        # Add lines and points
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        # Labels and styling
        labs(x = pred_label, 
            y = paste("Probability of", ifelse(response_var == "FemaleSong_Agg01", "Female Song", "Cooperative Breeding")),
            color = terr_label,
            fill = terr_label,
            subtitle = subtitle,
            title = paste(title_prefix, pred_label, "and", terr_label)) +
        scale_color_manual(values = c("darkblue", "darkred"),
                       labels = c("0", "1")) +
        scale_fill_manual(values = c("darkblue", "darkred"),
                      labels = c("0", "1")) +
        scale_x_discrete(labels = c("Absent", "Present")) +
        ylim(0, 1) +
        theme_bw() +
        theme(text = element_text(size = 12),
             plot.title = element_text(size = 12, face = "bold"),
             legend.position = "right",
             plot.subtitle = element_text(size = 6))
      
      message("Successfully created interaction plot")
      return(p)
    } else if (has_pred_var) {
      message("Creating single predictor plot")
      # Create a single predictor plot
      p <- ggplot(new_data, aes(x = factor(.data[[pred_var]]), y = pred_prob,
                           group = 1)) +
        # Add confidence interval
        geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
        # Add line and points
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        # Labels and styling
        labs(x = pred_label, 
            y = paste("Probability of", ifelse(response_var == "FemaleSong_Agg01", "Female Song", "Cooperative Breeding")),
            title = paste("Effect of", pred_label, "and", terr_label),
            subtitle = subtitle) +
        scale_x_discrete(labels = c("Absent", "Present")) +
        ylim(0, 1) +
        theme_bw() +
        theme(text = element_text(size = 12),
             plot.title = element_text(size = 12, face = "bold"))
      
      message("Successfully created predictor plot")
      return(p)
    } else if (has_terr_var) {
      message("Creating single territoriality plot")
      # Create a single territoriality plot
      p <- ggplot(new_data, aes(x = factor(.data[[terr_var]]), y = pred_prob,
                           group = 1)) +
        # Add confidence interval
        geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
        # Add line and points
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        # Labels and styling
        labs(x = terr_label, 
            y = paste("Probability of", ifelse(response_var == "FemaleSong_Agg01", "Female Song", "Cooperative Breeding")),
            title = paste("Effect of", pred_label, "and", terr_label),
            subtitle = subtitle) +
        scale_x_discrete(labels = c("0", "1")) +
        ylim(0, 1) +
        theme_bw() +
        theme(text = element_text(size = 12),
             plot.title = element_text(size = 12, face = "bold"),
             plot.subtitle = element_text(size = 6))
      
      message("Successfully created territoriality plot")
      return(p)
    }
    
    # If we get here, there's nothing to plot
    message("No suitable plot type identified")
    return(NULL)
    
  }, error = function(e) {
    # Print full error and call stack
    message("Error in create_effect_plot function:")
    message(paste(e$message))
    message(paste("Call:", deparse(e$call)))
    message("This error has been handled, but no plot will be created for this model")
    return(NULL)
  })
}


###
#### Function to create summary of top model coefficients ----
###

create_coefficient_summary <- function(top_models, response_var, terr_var) {
  # Create a data frame to store coefficient estimates
  coef_summary <- data.frame()
  
  for (model_name in names(top_models)) {
    model <- top_models[[model_name]]
    
    # Extract model formula as a string
    formula_str <- deparse(formula(model))
    # If formula is split across multiple strings, combine them
    if (length(formula_str) > 1) {
      formula_str <- paste(formula_str, collapse = " ")
    }
    
    # Extract model coefficients and their significance
    coefs <- tryCatch({
      summary(model)$coefficients
    }, error = function(e) {
      message(paste("Error extracting coefficients from model", model_name, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(coefs)) {
      # Add model information
      model_data <- data.frame(
        Response = response_var,
        Territoriality = terr_var,
        Model = model_name,
        Formula = formula_str,  # Add the formula as a string
        Parameter = rownames(coefs),
        Estimate = coefs[, "Estimate"],
        StdErr = coefs[, "StdErr"],
        z_value = coefs[, "z.value"],
        p_value = coefs[, "p.value"],
        Significance = case_when(
          coefs[, "p.value"] < 0.001 ~ "***",
          coefs[, "p.value"] < 0.01 ~ "**",
          coefs[, "p.value"] < 0.05 ~ "*",
          coefs[, "p.value"] < 0.1 ~ ".",
          TRUE ~ ""
        ),
        AIC = AIC(model)
      )
      
      # Append to the summary data frame
      coef_summary <- rbind(coef_summary, model_data)
    }
  }
  
  return(coef_summary)
}

###
#### Function to perform model averaging ----
###

perform_model_averaging <- function(model_results, threshold = 2) {
  # Get models with deltaAIC < threshold
  model_subset <- model_results$comparison %>%
    filter(deltaAIC < threshold)
  
  # Extract model names and weights
  model_names <- model_subset$Model
  model_formulas <- model_subset$Formula
  weights <- model_subset$weight / sum(model_subset$weight)  # Renormalize weights
  
  # Get all models
  models <- model_results$models[model_names]
  
  # Get all unique parameters across models
  all_params <- unique(unlist(lapply(models, function(m) names(coef(m)))))
  
  # Create matrix to store coefficients
  coef_matrix <- matrix(0, nrow = length(all_params), ncol = length(models))
  rownames(coef_matrix) <- all_params
  colnames(coef_matrix) <- model_names
  
  # Fill in coefficient matrix
  for (i in seq_along(models)) {
    model_name <- model_names[i]
    model_coefs <- coef(models[[model_name]])
    for (param in names(model_coefs)) {
      coef_matrix[param, model_name] <- model_coefs[param]
    }
  }
  
  # Calculate weighted average coefficients
  avg_coefs <- coef_matrix %*% weights
  
  # Create a data frame with results
  avg_coef_df <- data.frame(
    Parameter = rownames(coef_matrix),
    Estimate = avg_coefs,
    stringsAsFactors = FALSE
  )
  
  # Add formula information
  avg_coef_df$AvgFormula <- paste("Model average of:", paste(model_formulas, collapse = "; "))
  
  # Add individual model coefficients, formulas, and weights for reference
  for (i in seq_along(model_names)) {
    avg_coef_df[[paste0("Model", i)]] <- coef_matrix[, i]
    avg_coef_df[[paste0("Formula", i)]] <- model_formulas[i]
    avg_coef_df[[paste0("Weight", i)]] <- weights[i]
  }
  
  return(avg_coef_df)
}


###
#### Function to calculate marginal effects ----
###
calculate_marginal_effects <- function(model, data, pred_var) {
  # Only proceed if the predictor is in the model
  if (!(pred_var %in% names(coef(model)))) {
    return(NA)
  }
  
  # Create two versions of the data
  data_pred_0 <- data
  data_pred_1 <- data
  
  # Set predictor to 0 in first dataset
  data_pred_0[[pred_var]] <- 0
  
  # Set predictor to 1 in second dataset
  data_pred_1[[pred_var]] <- 1
  
  # Get linear predictors for each
  X0 <- model.matrix(formula(model), data_pred_0)
  X1 <- model.matrix(formula(model), data_pred_1)
  
  # Get predicted probabilities
  probs_0 <- plogis(X0 %*% coef(model))
  probs_1 <- plogis(X1 %*% coef(model))
  
  # Calculate average marginal effect
  marg_effect <- mean(probs_1 - probs_0)
  
  return(marg_effect)
}

###
#### Function for predicted probability differences by scenario ----
###

create_probability_scenarios <- function() {
  # Create an empty data frame
  scenarios <- NULL
  
  for (key in names(all_top_models)) {
    parts <- strsplit(key, "_")[[1]]
    resp_var <- parts[1]
    terr_var <- paste(parts[-1], collapse = "_")
    
    # Determine predictor based on response
    pred_var <- ifelse(resp_var == "FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
    
    # For each top model
    for (model_name in names(all_top_models[[key]])) {
      model <- all_top_models[[key]][[model_name]]
      
      # Get mean mass value
      mean_mass <- mean(all_model_results[[key]]$data$logMass_AVONET, na.rm = TRUE)
      
      # Create scenarios with generic column names
      temp_data <- expand.grid(
        logMass_AVONET = mean_mass,
        predictor_value = c(0, 1),
        terr_value = c(0, 1)
      )
      
      # Extract coefficients
      coefs <- coef(model)
      
      # Calculate linear predictor
      temp_data$linear_pred <- coefs["(Intercept)"]
      
      # Add main effects if present
      if (pred_var %in% names(coefs)) {
        temp_data$linear_pred <- temp_data$linear_pred + coefs[pred_var] * temp_data$predictor_value
      }
      
      if (terr_var %in% names(coefs)) {
        temp_data$linear_pred <- temp_data$linear_pred + coefs[terr_var] * temp_data$terr_value
      }
      
      if ("logMass_AVONET" %in% names(coefs)) {
        temp_data$linear_pred <- temp_data$linear_pred + coefs["logMass_AVONET"] * temp_data$logMass_AVONET
      }
      
      # Add interaction if present
      interaction_term <- paste0(pred_var, ":", terr_var)
      if (interaction_term %in% names(coefs)) {
        temp_data$linear_pred <- temp_data$linear_pred + 
          coefs[interaction_term] * temp_data$predictor_value * temp_data$terr_value
      }
      
      # Calculate predicted probabilities
      temp_data$pred_prob <- plogis(temp_data$linear_pred)
      
      # Create standardized output dataframe with consistent column names
      result_data <- data.frame(
        Response = resp_var,
        Territoriality = terr_var,
        Model = model_name,
        Formula = paste(deparse(formula(model)), collapse = " "),
        AIC = AIC(model),
        logMass_AVONET = temp_data$logMass_AVONET,
        Predictor = pred_var,
        Predictor_Value = temp_data$predictor_value,
        Territoriality_Variable = terr_var,
        Territoriality_Value = temp_data$terr_value,
        Predicted_Probability = temp_data$pred_prob,
        stringsAsFactors = FALSE
      )
      
      # Create scenario descriptions
      result_data$Scenario <- paste(
        ifelse(result_data$Predictor_Value == 0, paste0(pred_var, ": Absent"), paste0(pred_var, ": Present")),
        ifelse(result_data$Territoriality_Value == 0, paste0(terr_var, ": 0"), paste0(terr_var, ": 1")),
        sep = ", "
      )
      
      # Initialize scenarios if it's NULL
      if (is.null(scenarios)) {
        scenarios <- result_data
      } else {
        # Append to results, now with matching column names
        scenarios <- rbind(scenarios, result_data)
      }
    }
  }
  
  return(scenarios)
}

###
#### Main Analysis ----
###

# Define the territoriality variables to test
territoriality_vars <- c(
  "TerritorialityWeakVsStrong", 
  "TerritorialityWeakVsStrongHighConf", 
  "Territory_12vs3", 
  "TerritorialityPermissiveExclusive", 
  "TerritorialityPermissiveExclusiveHighConf", 
  "TerritorialityPermissiveColonialCoopVsExclusive"
)

# Define the response variables
response_vars <- c("FemaleSong_Agg01", "HighConfidence_Coop")

# Create lists to store results
all_model_results <- list()
all_top_models <- list()
all_plots <- list()
all_coefficients <- list()
all_model_averages <- list()

# Run all analyses
for (resp_var in response_vars) {
  for (terr_var in territoriality_vars) {
    # Skip if the territoriality variable doesn't exist
    if (!(terr_var %in% names(dfIn_clean))) {
      message(paste("Skipping", terr_var, "- variable not found in data"))
      next
    }
    
    # Create a unique key for this analysis
    analysis_key <- paste(resp_var, terr_var, sep = "_")
    
    # Run the initial models
    message(paste("\nRunning models for response:", resp_var, "and territoriality:", terr_var))
    all_model_results[[analysis_key]] <- run_model_set(resp_var, terr_var, dfIn_clean, tree_clean)
    
    # Run the top models with bootstrapping
    message(paste("Running top models with bootstrapping for:", analysis_key))
    all_top_models[[analysis_key]] <- run_top_models_with_boot(
      all_model_results[[analysis_key]], 
      n_boot = 100
    )
    
    # Create interaction plots for relevant top models
    for (model_name in names(all_top_models[[analysis_key]])) {
      plot_key <- paste(analysis_key, model_name, sep = "_")
      message(paste("Attempting to create plot for", plot_key))
      
      tryCatch({
        all_plots[[plot_key]] <- create_effect_plot(
          all_top_models[[analysis_key]][[model_name]], 
          all_model_results[[analysis_key]]$data, 
          resp_var, 
          terr_var
        )
        
        if (!is.null(all_plots[[plot_key]])) {
          message(paste("Successfully created plot for", plot_key))
        } else {
          message(paste("No plot was generated for", plot_key))
        }
      }, error = function(e) {
        message(paste("Error creating plot for", plot_key, ":", e$message))
        message("Skipping this plot and continuing with others")
        all_plots[[plot_key]] <- NULL
      })
    }
    
    # Create coefficient summary
    all_coefficients[[analysis_key]] <- create_coefficient_summary(
      all_top_models[[analysis_key]],
      resp_var,
      terr_var
    )
    
    # Perform model averaging
    all_model_averages[[analysis_key]] <- perform_model_averaging(
      all_model_results[[analysis_key]]
    )
  }
}

###
#### Create summary tables and visualizations ----
###

##### Combine all model comparisons into one table ----
all_comparisons <- data.frame()
for (key in names(all_model_results)) {
  # Extract response and territoriality variables from the key
  parts <- strsplit(key, "_")[[1]]
  resp_var <- parts[1]
  terr_var <- paste(parts[-1], collapse = "_")  # Rejoin in case territoriality has underscores
  
  # Add this information to the comparison
  comparison <- all_model_results[[key]]$comparison
  comparison$Response <- resp_var
  comparison$Territoriality <- terr_var
  
  # Append to the combined table
  all_comparisons <- rbind(all_comparisons, comparison)
}

# Sort by response, territoriality, and AIC
all_comparisons <- all_comparisons %>%
  arrange(Response, Territoriality, AIC)
# Save the results
write.csv(all_comparisons, "all_model_comparisons.csv", row.names = FALSE)


##### Calculate odds ratios and confidence intervals and summaries of coefficients (formerly "all_coefficient_summaries.csv") ----
all_coef_summary <- data.frame()

for (key in names(all_coefficients)) {
  coef_data <- all_coefficients[[key]]
  
  # Skip if empty
  if (nrow(coef_data) == 0) next
  
  # Calculate odds ratios
  coef_data$OddsRatio <- exp(coef_data$Estimate)
  coef_data$OR_LowerCI <- exp(coef_data$Estimate - 1.96 * coef_data$StdErr)
  coef_data$OR_UpperCI <- exp(coef_data$Estimate + 1.96 * coef_data$StdErr)
  
  # Add to combined data
  all_coef_summary <- rbind(all_coef_summary, coef_data)
}
# Write to CSV
write.csv(all_coef_summary, "coefficient_summaries_effect_sizes.csv", row.names = FALSE)


all_coef_summary %>% group_by(Response, Model) %>% summarize(NumberOfTerrVars = length(unique(Territoriality)))
all_coef_summary %>% group_by(Response, Parameter) %>% summarize(NumberOfTerrVars = length(unique(Territoriality)), median_p = median(p_value), nSig05 = sum(p_value < 0.05), nSig1 = sum(p_value < 0.1), nPresent = n(), minCoef = min(Estimate), maxCoef = max(Estimate)) %>% print(n=35)
all_coef_summary %>% group_by(Response, Model, Parameter) %>% summarize(NumberOfTerrVars = length(unique(Territoriality)), median_p = median(p_value), nSig05 = sum(p_value < 0.05), nSig1 = sum(p_value < 0.1), nPresent = n(), minCoef = min(Estimate), maxCoef = max(Estimate)) %>% print(n=95)

# Filter for parameters of interest
or_plot_data <- coef_data %>%
  filter(Parameter %in% c("HighConfidence_Coop", "TerritorialityPermissiveExclusiveHighConf"))

# Create odds ratio plot
ggplot(or_plot_data, aes(x = Model, y = OddsRatio, color = Parameter)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = OR_LowerCI, ymax = OR_UpperCI), 
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  theme_bw() +
  labs(title = "Odds Ratios Across Models",
       x = "Model", y = "Odds Ratio (log scale)")


##### Summarize effect sizes by territoriality variable ----
effect_size_summary <- all_coef_summary %>%
  filter(Parameter %in% c("HighConfidence_Coop", "FemaleSong_Agg01", 
                          "TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf",
                          "Territory_12vs3", "TerritorialityPermissiveExclusive",
                          "TerritorialityPermissiveExclusiveHighConf",
                          "TerritorialityPermissiveColonialCoopVsExclusive")) %>%
  group_by(Response, Territoriality, Parameter) %>%
  summarize(
    n_models = n(),
    mean_OR = mean(OddsRatio),
    min_OR = min(OddsRatio),
    max_OR = max(OddsRatio),
    sig_models = sum(p_value < 0.05),
    pct_significant = sig_models / n_models * 100
  ) %>%
  arrange(Response, Territoriality, Parameter)

# Write summary
write.csv(effect_size_summary, "effect_size_summary.csv", row.names = FALSE)


##### Calculate marginal effects using helper function ----
marginal_effects <- data.frame()

for (key in names(all_top_models)) {
  parts <- strsplit(key, "_")[[1]]
  resp_var <- parts[1]
  terr_var <- paste(parts[-1], collapse = "_")
  
  # Determine the predictor variable
  pred_var <- ifelse(resp_var == "FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
  
  # Apply to each model
  for (model_name in names(all_top_models[[key]])) {
    model <- all_top_models[[key]][[model_name]]
    data <- all_model_results[[key]]$data
    
    # Calculate marginal effects for both main predictors
    marg_effect_pred <- calculate_marginal_effects(model, data, pred_var)
    marg_effect_terr <- calculate_marginal_effects(model, data, terr_var)
    
    # Get the formula
    formula_str <- deparse(formula(model))
    
    # Create a row for this model
    model_row <- data.frame(
      Response = resp_var,
      Territoriality = terr_var,
      Model = model_name,
      Formula = formula_str,
      Predictor = pred_var,
      ME_Predictor = marg_effect_pred,
      ME_Territoriality = marg_effect_terr,
      AIC = AIC(model)
    )
    
    # Add to results
    marginal_effects <- rbind(marginal_effects, model_row)
  }
}
# Write to CSV
write.csv(marginal_effects, "marginal_effects.csv", row.names = FALSE)


##### Predicted probability differences ----
# Generate and save scenario-based probabilities
predicted_probabilities <- create_probability_scenarios() # doesn't currently seem correct - all "Scenarios" within each model have the same Predicted_Probability
#write.csv(predicted_probabilities, "predicted_probabilities.csv", row.names = FALSE)


##### Revised code for combining model averages ----
all_avg_summary <- data.frame()

for (key in names(all_model_averages)) {
  # Extract response and territoriality variables from the key
  parts <- strsplit(key, "_")[[1]]
  resp_var <- parts[1]
  terr_var <- paste(parts[-1], collapse = "_")
  
  # Add this information to the averages
  averages <- all_model_averages[[key]]
  
  # Only keep essential columns to ensure consistent structure
  essential_cols <- c("Parameter", "Estimate", "AvgFormula")
  if (all(essential_cols %in% colnames(averages))) {
    # Keep only essential columns plus response and territoriality
    averages_slim <- averages[, essential_cols]
    averages_slim$Response <- resp_var
    averages_slim$Territoriality <- terr_var
    
    # Store the number of models used in averaging
    model_columns <- grep("^Model[0-9]+$", colnames(averages))
    averages_slim$ModelsAveraged <- length(model_columns)
    
    # Add maximum weight information
    weight_columns <- grep("^Weight[0-9]+$", colnames(averages))
    if (length(weight_columns) > 0) {
      averages_slim$MaxWeight <- max(unlist(averages[, weight_columns]))
    } else {
      averages_slim$MaxWeight <- NA
    }
    
    # Extract formula for highest weighted model
    if (length(weight_columns) > 0) {
      max_weight_idx <- which.max(unlist(averages[1, weight_columns]))
      formula_col <- paste0("Formula", max_weight_idx)
      if (formula_col %in% colnames(averages)) {
        averages_slim$TopModelFormula <- averages[1, formula_col]
      }
    }
    
    # Append to the combined table
    all_avg_summary <- rbind(all_avg_summary, averages_slim)
  } else {
    warning(paste("Skipping", key, "- missing essential columns"))
  }
}

# Save the results - now with formulas included
write.csv(all_avg_summary, "all_model_averages.csv", row.names = FALSE)

##### Also save individual model averaging results with full details ----
for (key in names(all_model_averages)) {
  # Add response and territoriality info
  parts <- strsplit(key, "_")[[1]]
  resp_var <- parts[1]
  terr_var <- paste(parts[-1], collapse = "_")
  
  averages <- all_model_averages[[key]]
  averages$Response <- resp_var
  averages$Territoriality <- terr_var
  
  # Create a clean filename
  clean_key <- gsub(":", "_", key)
  filename <- paste0("model_avg_", clean_key, ".csv")
  
  # Save individual result
  write.csv(averages, filename, row.names = FALSE)
}

##### Create a combined PDF with all interaction plots ----
pdf("all_interaction_plots2.pdf", width = 10, height = 8)
for (plot_key in names(all_plots)) {
  if (!is.null(all_plots[[plot_key]])) {
    print(all_plots[[plot_key]])
  }
}
dev.off()

##### Create summary of best models for each analysis ----
best_models_summary <- data.frame()
for (key in names(all_model_results)) {
  # Extract response and territoriality variables
  parts <- strsplit(key, "_")[[1]]
  resp_var <- parts[1]
  terr_var <- paste(parts[-1], collapse = "_")
  
  # Get the best model
  best_model <- all_model_results[[key]]$comparison[1, ]
  best_model$Response <- resp_var
  best_model$Territoriality <- terr_var
  
  # Append to summary
  best_models_summary <- rbind(best_models_summary, best_model)
}

# Sort by response and AIC
best_models_summary <- best_models_summary %>%
  arrange(Response, AIC)

# Save this summary
write.csv(best_models_summary, "best_models_summary.csv", row.names = FALSE)


##### Create heat map plot of coefficients ----
library(reshape2)

# Create matrix of coefficients
coef_matrix <- coef_data %>%
  select(Model, Parameter, Estimate) %>%
  dcast(Parameter ~ Model, value.var = "Estimate")

# Create significance matrix
sig_matrix <- coef_data %>%
  select(Model, Parameter, p_value) %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ ""
  )) %>%
  dcast(Parameter ~ Model, value.var = "Significance")

# Plot heatmap
ggplot(melt(coef_matrix, id.vars = "Parameter"), 
       aes(x = variable, y = Parameter, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Coefficient") +
  geom_text(data = melt(sig_matrix, id.vars = "Parameter"),
            aes(x = variable, y = Parameter, label = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coefficients Across Models",
       x = "Model", y = "Parameter")



# Print a summary of findings
cat("\n\n=== SUMMARY OF FINDINGS ===\n\n")
cat("Total number of analyses run:", length(all_model_results), "\n")
cat("Total number of top models with bootstrapping:", sum(sapply(all_top_models, length)), "\n")
cat("Total number of interaction plots created:", sum(!sapply(all_plots, is.null)), "\n\n")

# Print best models for each response variable
cat("Best models for predicting Female Song:\n")
print(best_models_summary %>% filter(Response == "FemaleSong_Agg01") %>% 
        select(Territoriality, Model, Formula, AIC, deltaAIC, weight, n_species) %>%
        arrange(AIC))

cat("\nBest models for predicting Cooperative Breeding:\n")
print(best_models_summary %>% filter(Response == "HighConfidence_Coop") %>% 
        select(Territoriality, Model, Formula, AIC, deltaAIC, weight, n_species) %>%
        arrange(AIC))

