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