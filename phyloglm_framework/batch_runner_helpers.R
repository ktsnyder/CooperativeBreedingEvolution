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
  
  # Exclude intercept for effect sizes
  effects <- effects[effects$Parameter != "(Intercept)", ]
  
  return(effects)
}

#' Model averaging for effect estimates
#' 
#' @param comparison Model comparison results
#' @param coefficients All coefficients
#' @param threshold deltaAIC threshold for inclusion
#' @return Data frame with model-averaged effects
model_average_effects <- function(comparison, coefficients, threshold = 2) {
  
  # Get models within threshold
  top_models <- comparison$comparison[comparison$comparison$deltaAIC < threshold, ]
  top_model_names <- top_models$Model
  weights <- top_models$weight / sum(top_models$weight)  # Renormalize
  
  # Get unique parameters across top models
  top_coefs <- coefficients[coefficients$Model %in% top_model_names, ]
  unique_params <- unique(top_coefs$Parameter)
  unique_params <- unique_params[unique_params != "(Intercept)"]
  
  # Average each parameter
  averaged <- data.frame()
  
  for (param in unique_params) {
    # Get all estimates for this parameter
    param_coefs <- top_coefs[top_coefs$Parameter == param, ]
    
    # Match with model weights
    param_weights <- numeric(nrow(param_coefs))
    for (i in 1:nrow(param_coefs)) {
      model_idx <- which(top_models$Model == param_coefs$Model[i])
      param_weights[i] <- weights[model_idx]
    }
    
    # If parameter not in all models, adjust weights
    if (nrow(param_coefs) < length(top_model_names)) {
      # Models without this parameter contribute 0
      total_weight <- sum(param_weights)
    } else {
      total_weight <- 1
    }
    
    # Calculate weighted average
    avg_estimate <- sum(param_coefs$Estimate * param_weights) / total_weight
    
    # Weighted SE (Burnham & Anderson approximation)
    avg_se <- sqrt(sum(param_weights * (param_coefs$StdErr^2 + 
                                        (param_coefs$Estimate - avg_estimate)^2)) / total_weight)
    
    averaged <- rbind(averaged, data.frame(
      Parameter = param,
      Estimate = avg_estimate,
      StdErr = avg_se,
      CI_lower = avg_estimate - 1.96 * avg_se,
      CI_upper = avg_estimate + 1.96 * avg_se,
      N_models = nrow(param_coefs),
      Total_weight = total_weight,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add odds ratios
  averaged$OddsRatio <- exp(averaged$Estimate)
  averaged$OR_CI_lower <- exp(averaged$CI_lower)
  averaged$OR_CI_upper <- exp(averaged$CI_upper)
  
  return(averaged)
}

#' Create analysis plots
#' 
#' @param model_results Model fitting results
#' @param comparison Model comparison results
#' @param coefficients Coefficient data frame
#' @param prepared_data Prepared data object
#' @param output_dir Directory for plots
create_analysis_plots <- function(model_results, comparison, coefficients,
                                 prepared_data, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # 1. Model comparison plot
  library(ggplot2)
  
  comp_plot <- ggplot(comparison$comparison, 
                      aes(x = reorder(Model, -deltaAIC), y = deltaAIC)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(x = "Model", y = "ΔAIC", 
         title = "Model Comparison",
         subtitle = "Models within 2 AIC units are considered equivalent") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "model_comparison.png"), 
         comp_plot, width = 8, height = 6)
  
  # 2. Coefficient plot for best model
  best_model_name <- comparison$comparison$Model[1]
  best_coefs <- coefficients[coefficients$Model == best_model_name & 
                            coefficients$Parameter != "(Intercept)", ]
  
  if (nrow(best_coefs) > 0) {
    coef_plot <- ggplot(best_coefs, aes(x = Parameter, y = Estimate)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Estimate - 1.96 * StdErr,
                       ymax = Estimate + 1.96 * StdErr),
                   width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      coord_flip() +
      labs(x = "", y = "Coefficient Estimate",
           title = paste("Coefficients:", best_model_name)) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "best_model_coefficients.png"),
           coef_plot, width = 8, height = 6)
  }
  
  # 3. Effect size plot (odds ratios)
  effects <- calculate_effect_sizes(coefficients)
  top_effects <- effects[effects$Model %in% 
                        comparison$comparison$Model[comparison$comparison$deltaAIC < 2], ]
  
  if (nrow(top_effects) > 0) {
    or_plot <- ggplot(top_effects, aes(x = Parameter, y = OddsRatio, color = Model)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = OR_CI_lower, ymax = OR_CI_upper),
                   position = position_dodge(width = 0.5),
                   width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      scale_y_log10() +
      coord_flip() +
      labs(x = "", y = "Odds Ratio (95% CI)",
           title = "Effect Sizes from Top Models") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "effect_sizes.png"),
           or_plot, width = 8, height = 8)
  }
}

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