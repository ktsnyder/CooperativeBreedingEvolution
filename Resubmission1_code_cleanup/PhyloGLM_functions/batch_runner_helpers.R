# Helper functions for PhyloGLM Batch Runner

#' Bootstrap top models
#' 
#' @param top_models List of top model objects
#' @param data Prepared data
#' @param tree Prepared tree
#' @param n_boot Number of bootstrap replicates
#' @return List with bootstrap results
bootstrap_top_models <- function(top_models, data, tree, n_boot = 100) {
  
  boot_results <- list()
  
  for (model_name in names(top_models)) {
    model <- top_models[[model_name]]
    
    cat(paste("    Bootstrapping", model_name, "..."))
    
    # Refit with bootstrap
    boot_fit <- tryCatch({
      phyloglm(
        formula = formula(model),
        data = data,
        phy = tree,
        method = "logistic_MPLE",
        btol = 30,
        log.alpha.bound = 4,
        boot = n_boot
      )
    }, error = function(e) {
      cat(" ERROR\n")
      NULL
    })
    
    if (!is.null(boot_fit)) {
      cat(" Done\n")
      boot_results[[model_name]] <- boot_fit
    }
  }
  
  return(boot_results)
}

#' Extract coefficients from all models
#' 
#' @param model_results Output from fit_all_models
#' @param boot_results Optional bootstrap results
#' @return Data frame with all coefficients
extract_all_coefficients <- function(model_results, boot_results = NULL) {
  
  models <- model_results$models
  all_coefs <- data.frame()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Get coefficients
    coef_summary <- summary(model)$coefficients
    
    # Convert to data frame
    coef_df <- data.frame(
      Model = model_name,
      AIC = AIC(model),  # Add AIC for each model
      Parameter = rownames(coef_summary),
      Estimate = coef_summary[, "Estimate"],
      StdErr = coef_summary[, "StdErr"],
      z_value = coef_summary[, "z.value"],
      p_value = coef_summary[, "p.value"],
      CI_lower = NA,  # Initialize CI columns
      CI_upper = NA,
      stringsAsFactors = FALSE
    )
    
    # Add bootstrap CIs if available
    if (!is.null(boot_results) && model_name %in% names(boot_results)) {
      boot_model <- boot_results[[model_name]]
      if (!is.null(boot_model$bootconfint95)) {
        ci <- boot_model$bootconfint95
        # Match by parameter name
        for (i in 1:nrow(coef_df)) {
          param_name <- coef_df$Parameter[i]
          if (param_name %in% rownames(ci)) {
            coef_df$CI_lower[i] <- ci[param_name, 1]
            coef_df$CI_upper[i] <- ci[param_name, 2]
          }
        }
      }
    }
    
    all_coefs <- rbind(all_coefs, coef_df)
  }
  
  # Add significance stars
  all_coefs$Significance <- ifelse(all_coefs$p_value < 0.001, "***",
                                   ifelse(all_coefs$p_value < 0.01, "**",
                                         ifelse(all_coefs$p_value < 0.05, "*",
                                               ifelse(all_coefs$p_value < 0.1, ".", ""))))
  
  return(all_coefs)
}

#' Calculate effect sizes
#' 
#' @param coefficients Data frame of coefficients
#' @return Data frame with effect sizes
calculate_effect_sizes <- function(coefficients) {
  
  # For logistic regression, convert to odds ratios
  effects <- coefficients
  effects$OddsRatio <- exp(effects$Estimate)
  
  # Always calculate CIs - use bootstrap if available, otherwise normal approximation
  if ("CI_lower" %in% colnames(effects) && !all(is.na(effects$CI_lower))) {
    # Use bootstrap CIs where available
    effects$OR_CI_lower <- ifelse(!is.na(effects$CI_lower),
                                  exp(effects$CI_lower),
                                  exp(effects$Estimate - 1.96 * effects$StdErr))
    effects$OR_CI_upper <- ifelse(!is.na(effects$CI_upper),
                                  exp(effects$CI_upper),
                                  exp(effects$Estimate + 1.96 * effects$StdErr))
  } else {
    # Use normal approximation for all
    effects$OR_CI_lower <- exp(effects$Estimate - 1.96 * effects$StdErr)
    effects$OR_CI_upper <- exp(effects$Estimate + 1.96 * effects$StdErr)
  }
  
  # Add CI width for diagnostics
  effects$CI_width <- effects$OR_CI_upper - effects$OR_CI_lower
  
  # # Exclude intercept for effect sizes
  # effects <- effects[effects$Parameter != "(Intercept)", ]
  
  return(effects)
}

#' Create best model effects summary
#' 
#' @param effects Data frame of effect sizes for all models
#' @param comparison Model comparison results
#' @param boot_results Bootstrap results for empirical p-values
#' @return Data frame with best model effects only
create_best_model_effects <- function(effects, comparison, boot_results = NULL) {
  
  # Get the best model (lowest AIC)
  best_model_name <- comparison$comparison$Model[1]
  best_effects <- effects[effects$Model == best_model_name, ]
  
  if (nrow(best_effects) == 0) {
    return(data.frame())  # Return empty if no effects
  }
  
  # Create output with requested columns
  best_summary <- data.frame(
    Model = best_effects$Model,
    AIC = best_effects$AIC,
    Parameter = best_effects$Parameter,
    OddsRatio = best_effects$OddsRatio,
    OR_CI_Lower = best_effects$OR_CI_lower,
    OR_CI_Upper = best_effects$OR_CI_upper,
    stringsAsFactors = FALSE
  )
  
  # Calculate empirical p-values from bootstrap if available
  if (!is.null(boot_results) && best_model_name %in% names(boot_results)) {
    boot_model <- boot_results[[best_model_name]]

    if (!is.null(boot_model$bootstrap)) {
      # Extract bootstrap coefficient estimates
      boot_coefs <- boot_model$bootstrap
      
      # Calculate empirical p-values for each parameter
      best_summary$p_value <- NA
      
      for (i in 1:nrow(best_summary)) {
        param_name <- best_summary$Parameter[i]
        
        if (param_name %in% colnames(boot_coefs)) {
          # Get bootstrap estimates for this parameter
          boot_estimates <- boot_coefs[, param_name ]
          # Convert to odds ratios
          boot_ors <- exp(boot_estimates)
          
          # Count how many bootstrap ORs are on the opposite side of 1.0 from the point estimate
          point_or <- best_summary$OddsRatio[i]
          
          if (point_or > 1.0) {
            # If point estimate OR > 1, count bootstrap ORs <= 1
            n_opposite <- sum(boot_ors <= 1.0, na.rm = TRUE)
          } else {
            # If point estimate OR < 1, count bootstrap ORs >= 1  
            n_opposite <- sum(boot_ors >= 1.0, na.rm = TRUE)
          }
          
          # Empirical p-value = 2 * (fraction on opposite side)
          best_summary$p_value[i] <- 2 * (n_opposite / length(boot_ors))
        } else {
          print(param_name)
          print("issue calculating effect p-values")
        }
      }
    }
  }
  
  return(best_summary)
}


# Removed model_average_effects()

# Removed create_analysis_plots()

#' Integrate results across batch analyses
#' 
#' @param all_results List of all analysis results
#' @param output_dir Directory for integrated outputs
#' @return List of integrated summaries
integrate_batch_results <- function(all_results, output_dir) {
  
  # Extract successful analyses
  successful <- all_results[sapply(all_results, function(r) r$success)]
  
  if (length(successful) == 0) {
    warning("No successful analyses to integrate")
    return(NULL)
  }
  
  # 1. Combine all model comparisons
  all_comparisons <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    if (!is.null(result$comparison) && !is.null(result$comparison$comparison) && 
        nrow(result$comparison$comparison) > 0) {
      comp <- result$comparison$comparison
      comp$Analysis <- name
      all_comparisons <- rbind(all_comparisons, comp)
    }
  }
  
  write.csv(all_comparisons, 
            file.path(output_dir, "all_model_comparisons.csv"),
            row.names = FALSE)
  
  # 2. Combine all coefficients
  all_coefs <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    if (!is.null(result$coefficients) && nrow(result$coefficients) > 0) {
      coefs <- result$coefficients
      coefs$Analysis <- name
      all_coefs <- rbind(all_coefs, coefs)
    }
  }
  
  write.csv(all_coefs,
            file.path(output_dir, "all_coefficients.csv"),
            row.names = FALSE)
  
  # 3. Create effect size comparison
  effect_comparison <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    if (!is.null(result$effects)) {
      # Get effects from best model
      best_model <- result$comparison$comparison$Model[1]
      best_effects <- result$effects[result$effects$Model == best_model, ]
      
      # Only add if there are effects
      if (nrow(best_effects) > 0) {
        best_effects$Analysis <- name
        effect_comparison <- rbind(effect_comparison, best_effects)
      } else {
        cat("Warning: No effects found for best model", best_model, "in analysis", name, "\n")
      }
    }
  }
  
  write.csv(effect_comparison,
            file.path(output_dir, "effect_comparison.csv"),
            row.names = FALSE)
  
  # 4. Create summary statistics
  summary_stats <- data.frame(
    Analysis = names(successful),
    N_species = sapply(successful, function(r) r$prepared_data$n_species),
    N_models = sapply(successful, function(r) length(r$models$models)),
    Best_model = sapply(successful, function(r) r$comparison$comparison$Model[1]),
    Best_AIC = sapply(successful, function(r) r$comparison$comparison$AIC[1]),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_stats,
            file.path(output_dir, "analysis_summary_stats.csv"),
            row.names = FALSE)
  
  return(list(
    comparisons = all_comparisons,
    coefficients = all_coefs,
    effects = effect_comparison,
    summary = summary_stats
  ))
}