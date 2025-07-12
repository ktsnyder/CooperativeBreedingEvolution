# PhyloGLM Summary Plots - Forest plot and interaction visualization
# Based on the reference figure design

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggtext)

#' Create forest plot showing bidirectional effects
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param predictor_var Main predictor variable name
#' @param response_var Response variable name
#' @param predictor_label Display label for predictor
#' @param response_label Display label for response
#' @param deltaAIC_threshold Threshold for including models
#' @return ggplot object
create_bidirectional_forest_plot <- function(batch_results,
                                           predictor_var = "HighConfidence_Coop",
                                           response_var = "FemaleSong_Agg01",
                                           predictor_label = "Cooperative Breeding",
                                           response_label = "Female Song",
                                           deltaAIC_threshold = 2) {
  
  # Extract all coefficients and model comparisons
  all_coefs <- data.frame()
  all_comparisons <- data.frame()
  
  for (name in names(batch_results$results)) {
    # Skip Territory_12vs3 analyses
    if (grepl("Terr3", name)) next
    
    result <- batch_results$results[[name]]
    
    # Only include successful analyses with territoriality
    if (result$success && 
        "TerritorialityWeakVsStrong" %in% result$config$predictors) {
      # Get coefficients
      coefs <- result$coefficients
      coefs$Analysis <- name
      all_coefs <- rbind(all_coefs, coefs)
      
      # Get model comparisons
      comps <- result$comparison$comparison
      comps$Analysis <- name
      all_comparisons <- rbind(all_comparisons, comps)
    }
  }
  
  # Identify analyses with our variables of interest
  # Forward direction: predictor -> response
  forward_analyses <- all_coefs %>%
    filter(Parameter == predictor_var) %>%
    pull(Analysis) %>%
    unique()
  
  # Reverse direction: response -> predictor  
  reverse_analyses <- all_coefs %>%
    filter(Parameter == response_var) %>%
    pull(Analysis) %>%
    unique()
  
  # Get effects for competitive models only
  bidirectional_effects <- rbind(
    # Forward effects
    all_coefs %>%
      filter(Analysis %in% forward_analyses,
             Parameter == predictor_var) %>%
      inner_join(all_comparisons %>% select(Analysis, Model, deltaAIC),
                 by = c("Analysis", "Model")) %>%
      filter(deltaAIC < deltaAIC_threshold) %>%
      mutate(Direction = paste0(predictor_label, " → ", response_label)),
    
    # Reverse effects
    all_coefs %>%
      filter(Analysis %in% reverse_analyses,
             Parameter == response_var) %>%
      inner_join(all_comparisons %>% select(Analysis, Model, deltaAIC),
                 by = c("Analysis", "Model")) %>%
      filter(deltaAIC < deltaAIC_threshold) %>%
      mutate(Direction = paste0(response_label, " → ", predictor_label))
  )
  
  # Process effects - only keep best models from analyses with territoriality
  bidirectional_effects <- bidirectional_effects %>%
    mutate(
      Is_Best = deltaAIC == 0,
      # Extract control variable from analysis name
      Control = case_when(
        grepl("Mass", Analysis) ~ "Mass",
        grepl("WingDim", Analysis) ~ "Wing Dimorphism",
        grepl("PlumDim", Analysis) ~ "Plumage Dimorphism",
        TRUE ~ "Other"
      ),
      PlotLabel = paste0(Direction, "\n(", Control, ")")
    ) %>%
    # Only keep BEST model per analysis
    group_by(Analysis) %>%
    filter(deltaAIC == min(deltaAIC)) %>%
    ungroup()
  
  # Calculate odds ratios with proper CIs
  bidirectional_effects <- bidirectional_effects %>%
    mutate(
      OddsRatio = exp(Estimate),
      # Use bootstrap CIs if available, otherwise use normal approximation
      OR_LowerCI = ifelse(!is.na(CI_lower), 
                          exp(CI_lower), 
                          exp(Estimate - 1.96 * StdErr)),
      OR_UpperCI = ifelse(!is.na(CI_upper), 
                          exp(CI_upper), 
                          exp(Estimate + 1.96 * StdErr))
    )
  
  # Create ordered factor for y-axis
  plot_order <- unique(bidirectional_effects$PlotLabel)
  bidirectional_effects$PlotLabel <- factor(bidirectional_effects$PlotLabel,
                                            levels = rev(plot_order))
  
  # Define colors
  terr_colors <- c("TerritorialityWeakVsStrong" = "#E64B35",
                   "Territory_12vs3" = "#4DBBD5",
                   "None" = "#999999")
  
  # Create forest plot
  p <- ggplot(bidirectional_effects,
              aes(y = reorder(PlotLabel, desc(as.numeric(factor(Direction)))), x = OddsRatio)) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = OR_LowerCI, xmax = OR_UpperCI),
                   height = 0.2,
                   linewidth = 0.8,
                   color = "#E64B35") +
    geom_point(aes(shape = p_value < 0.05),
               size = 3.5,
               color = "#E64B35") +
    # Add model labels with better positioning
    geom_text(aes(label = Model,
                  x = OR_UpperCI * 1.1),
              hjust = 0, size = 2.5, color = "gray40") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
                  limits = c(0.2, 20)) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
                       name = "p < 0.05") +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3), guide = "none") +
    labs(x = "Odds Ratio (95% CI)",
         y = "",
         title = paste0("A. Bidirectional Effects: ", predictor_label, " ⇄ ", response_label),
         subtitle = "Best models with territoriality moderation shown") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

#' Create interaction plot showing best model predictions
#' 
#' @param result Single analysis result containing the best interaction model
#' @param predictor_var Predictor variable name
#' @param terr_var Territoriality variable name
#' @param response_var Response variable name
#' @param predictor_label Display label for predictor
#' @param response_label Display label for response
#' @return ggplot object
create_interaction_plot <- function(result,
                                   predictor_var = "HighConfidence_Coop",
                                   terr_var = "TerritorialityWeakVsStrong",
                                   response_var = "FemaleSong_Agg01",
                                   predictor_label = "Cooperative Breeding",
                                   response_label = "Female Song") {
  
  # Find best interaction model
  best_model_name <- result$comparison$comparison$Model[1]
  
  # Check if it's actually an interaction model
  if (!grepl("x", best_model_name)) {
    # Find best interaction model
    interaction_models <- result$comparison$comparison %>%
      filter(grepl("x", Model)) %>%
      arrange(AIC)
    
    if (nrow(interaction_models) > 0) {
      best_model_name <- interaction_models$Model[1]
    } else {
      stop("No interaction model found")
    }
  }
  
  # Get the model
  model <- result$models$models[[best_model_name]]
  
  # Create prediction data
  pred_data <- expand.grid(
    predictor = c(0, 1),
    territory = c(0, 1)
  )
  names(pred_data) <- c(predictor_var, terr_var)
  
  # Add control variables at their means
  data <- result$prepared_data$data
  if ("logMass_AVONET" %in% names(data)) {
    pred_data$logMass_AVONET <- mean(data$logMass_AVONET, na.rm = TRUE)
  }
  
  # Get predictions - phyloglm doesn't have a predict method, so calculate manually
  # Extract coefficients
  coefs <- coef(model)
  
  # Create model matrix
  formula_obj <- formula(model)
  X <- model.matrix(formula_obj, data = pred_data)
  
  # Calculate linear predictor and convert to probabilities
  pred_data$predicted <- plogis(X %*% coefs)
  
  # Calculate confidence intervals using bootstrap if available
  if (!is.null(result$bootstrap) && best_model_name %in% names(result$bootstrap)) {
    boot_model <- result$bootstrap[[best_model_name]]
    
    if (!is.null(boot_model$boot)) {
      # Use bootstrap predictions
      n_boot <- nrow(boot_model$boot)
      boot_preds <- matrix(NA, nrow = nrow(pred_data), ncol = n_boot)
      
      for (i in 1:n_boot) {
        boot_coefs <- boot_model$boot[i, ]
        names(boot_coefs) <- names(coef(model))
        
        # Calculate linear predictor
        X <- model.matrix(formula(model), data = pred_data)
        boot_preds[, i] <- plogis(X %*% boot_coefs)
      }
      
      pred_data$lower_ci <- apply(boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE)
      pred_data$upper_ci <- apply(boot_preds, 1, quantile, probs = 0.975, na.rm = TRUE)
    }
  }
  
  # If no bootstrap CIs, use delta method approximation
  if (!("lower_ci" %in% names(pred_data))) {
    # Get variance-covariance matrix
    vcov_mat <- vcov(model)
    X <- model.matrix(formula(model), data = pred_data)
    
    # Calculate standard errors of predictions
    pred_se <- sqrt(diag(X %*% vcov_mat %*% t(X)))
    
    # Transform to probability scale
    logit_pred <- qlogis(pred_data$predicted)
    pred_data$lower_ci <- plogis(logit_pred - 1.96 * pred_se)
    pred_data$upper_ci <- plogis(logit_pred + 1.96 * pred_se)
  }
  
  # Prepare for plotting
  pred_data[[terr_var]] <- factor(pred_data[[terr_var]],
                                  levels = c(0, 1),
                                  labels = c("Weak/No Territory", "Strong Territory"))
  
  pred_data[[predictor_var]] <- factor(pred_data[[predictor_var]],
                                       levels = c(0, 1),
                                       labels = c("Absent", "Present"))
  
  # Extract interaction statistics
  coefs <- summary(model)$coefficients
  interaction_term <- paste0(predictor_var, ":", terr_var)
  
  if (interaction_term %in% rownames(coefs)) {
    int_beta <- coefs[interaction_term, "Estimate"]
    int_p <- coefs[interaction_term, "p.value"]
    interaction_text <- paste0("Interaction: β = ", round(int_beta, 2), 
                              ", p = ", format.pval(int_p, digits = 3))
  } else {
    interaction_text <- "No interaction term"
  }
  
  # Create interaction plot
  p <- ggplot(pred_data, aes_string(x = predictor_var, y = "predicted",
                                    color = terr_var, group = terr_var)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = .data[[terr_var]]),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = c("Weak/No Territory" = "#4DBBD5",
                                  "Strong Territory" = "#E64B35")) +
    scale_fill_manual(values = c("Weak/No Territory" = "#4DBBD5",
                                 "Strong Territory" = "#E64B35")) +
    labs(x = predictor_label,
         y = paste0("P(", response_label, ")"),
         color = "Territory",
         fill = "Territory",
         title = paste0("B. Best Interaction Model: ", best_model_name),
         subtitle = paste0("Formula: ", response_label, " ~ ", predictor_label, 
                          " * Territoriality + controls\n", interaction_text)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  return(p)
}

#' Create combined summary figure with bidirectional interaction plots
#' 
#' @param batch_results Results from batch analysis
#' @param forward_analysis Name of analysis for forward direction interaction
#' @param reverse_analysis Name of analysis for reverse direction interaction
#' @param ... Additional arguments passed to plotting functions
#' @return Combined patchwork plot
create_phyloglm_summary_figure <- function(batch_results,
                                          forward_analysis = NULL,
                                          reverse_analysis = NULL,
                                          predictor_var = "HighConfidence_Coop",
                                          response_var = "FemaleSong_Agg01",
                                          predictor_label = "Cooperative Breeding",
                                          response_label = "Female Song",
                                          output_file = NULL) {
  
  # Create forest plot
  panel_a <- create_bidirectional_forest_plot(
    batch_results,
    predictor_var = predictor_var,
    response_var = response_var,
    predictor_label = predictor_label,
    response_label = response_label
  )
  
  # Find forward analysis if not specified (predictor -> response)
  if (is.null(forward_analysis)) {
    for (name in names(batch_results$results)) {
      result <- batch_results$results[[name]]
      if (result$success && 
          result$config$response == response_var &&
          predictor_var %in% result$config$predictors &&
          "TerritorialityWeakVsStrong" %in% result$config$predictors &&
          !grepl("Terr3", name)) {  # Exclude Territory_12vs3
        forward_analysis <- name
        cat("Found forward analysis:", name, "\n")
        break
      }
    }
    if (is.null(forward_analysis)) {
      cat("Warning: No forward analysis found with", response_var, "~", predictor_var, "+ TerritorialityWeakVsStrong\n")
    }
  }
  
  # Find reverse analysis if not specified (response -> predictor)
  if (is.null(reverse_analysis)) {
    for (name in names(batch_results$results)) {
      result <- batch_results$results[[name]]
      if (result$success && 
          result$config$response == predictor_var &&
          response_var %in% result$config$predictors &&
          "TerritorialityWeakVsStrong" %in% result$config$predictors &&
          !grepl("Terr3", name)) {  # Exclude Territory_12vs3
        reverse_analysis <- name
        cat("Found reverse analysis:", name, "\n")
        break
      }
    }
    if (is.null(reverse_analysis)) {
      cat("Warning: No reverse analysis found with", predictor_var, "~", response_var, "+ TerritorialityWeakVsStrong\n")
    }
  }
  
  # Create forward interaction plot
  panel_b <- NULL
  if (!is.null(forward_analysis) && forward_analysis %in% names(batch_results$results)) {
    result <- batch_results$results[[forward_analysis]]
    
    terr_var <- if ("TerritorialityWeakVsStrong" %in% result$config$predictors) {
      "TerritorialityWeakVsStrong"
    } else {
      NULL
    }
    
    if (!is.null(terr_var)) {
      panel_b <- create_interaction_plot(
        result,
        predictor_var = predictor_var,
        terr_var = terr_var,
        response_var = response_var,
        predictor_label = predictor_label,
        response_label = response_label
      )
      # Update title
      panel_b <- panel_b + 
        labs(title = paste0("B. Best Interaction Model: ", predictor_label, " → ", response_label))
    }
  }
  
  # Create reverse interaction plot
  panel_c <- NULL
  if (!is.null(reverse_analysis) && reverse_analysis %in% names(batch_results$results)) {
    result <- batch_results$results[[reverse_analysis]]
    
    terr_var <- if ("TerritorialityWeakVsStrong" %in% result$config$predictors) {
      "TerritorialityWeakVsStrong"
    } else {
      NULL
    }
    
    if (!is.null(terr_var)) {
      panel_c <- create_interaction_plot(
        result,
        predictor_var = response_var,  # Swap roles
        terr_var = terr_var,
        response_var = predictor_var,  # Swap roles
        predictor_label = response_label,  # Swap labels
        response_label = predictor_label   # Swap labels
      )
      # Update title
      panel_c <- panel_c + 
        labs(title = paste0("C. Best Interaction Model: ", response_label, " → ", predictor_label))
    }
  }
  
  # Combine panels based on what's available
  if (!is.null(panel_b) && !is.null(panel_c)) {
    # All three panels
    combined <- panel_a / (panel_b | panel_c)
    combined <- combined + 
      plot_layout(heights = c(1, 1))
  } else if (!is.null(panel_b)) {
    # Just forward direction
    combined <- panel_a | panel_b
  } else if (!is.null(panel_c)) {
    # Just reverse direction
    combined <- panel_a | panel_c
  } else {
    # Forest plot only
    combined <- panel_a
  }
  
  combined <- combined + 
    plot_annotation(
      title = "Phylogenetic GLM Analysis: Cooperative Breeding and Female Song",
      subtitle = "Examining bidirectional relationships and territorial moderation",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )
  
  # Save if requested
  if (!is.null(output_file)) {
    height <- if (!is.null(panel_b) && !is.null(panel_c)) 10 else 6
    ggsave(output_file, combined, width = 14, height = height, dpi = 300)
  }
  
  return(combined)
}