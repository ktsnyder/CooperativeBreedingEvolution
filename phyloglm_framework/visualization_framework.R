# Multi-Analysis Visualization Framework for PhyloGLM
# Creates comparative plots across multiple analyses

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

#' Create a forest plot comparing effects across analyses
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param parameter_name Name of parameter to compare
#' @param analysis_labels Optional named vector for custom analysis labels
#' @return ggplot object
create_comparative_forest_plot <- function(batch_results, 
                                         parameter_name,
                                         analysis_labels = NULL) {
  
  # Extract successful analyses
  successful <- batch_results$results[sapply(batch_results$results, function(r) r$success)]
  
  # Collect effect sizes for the parameter
  effect_data <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    
    # Get best model effects
    if (!is.null(result$effects)) {
      best_model <- result$comparison$comparison$Model[1]
      param_effects <- result$effects[
        result$effects$Model == best_model & 
        result$effects$Parameter == parameter_name, 
      ]
      
      if (nrow(param_effects) > 0) {
        param_effects$Analysis <- name
        effect_data <- rbind(effect_data, param_effects[1, ])
      }
    }
  }
  
  if (nrow(effect_data) == 0) {
    stop(paste("No effects found for parameter:", parameter_name))
  }
  
  # Apply custom labels if provided
  if (!is.null(analysis_labels)) {
    effect_data$Analysis_Label <- ifelse(
      effect_data$Analysis %in% names(analysis_labels),
      analysis_labels[effect_data$Analysis],
      effect_data$Analysis
    )
  } else {
    effect_data$Analysis_Label <- effect_data$Analysis
  }
  
  # Create forest plot
  p <- ggplot(effect_data, aes(x = OddsRatio, y = reorder(Analysis_Label, OddsRatio))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), 
                   height = 0.2, size = 0.5) +
    geom_point(aes(color = p_value < 0.05), size = 3) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "gray50"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "") +
    labs(x = "Odds Ratio (95% CI)",
         y = "",
         title = paste("Comparative Effects:", parameter_name),
         subtitle = if (!is.null(batch_results$bootstrap_n)) {
           paste("Effects from best models in each analysis (", batch_results$bootstrap_n, " bootstrap iterations)", sep = "")
         } else {
           "Effects from best models in each analysis"
         }) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
  
  return(p)
}

#' Create a heatmap of model selection across analyses
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param top_n Number of top models to show per analysis
#' @return ggplot object
create_model_selection_heatmap <- function(batch_results, top_n = 5) {
  
  successful <- batch_results$results[sapply(batch_results$results, function(r) r$success)]
  
  # Collect model rankings
  model_data <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    
    if (!is.null(result$comparison)) {
      comp <- result$comparison$comparison
      comp <- comp[1:min(top_n, nrow(comp)), ]
      comp$Analysis <- name
      comp$Rank <- 1:nrow(comp)
      model_data <- rbind(model_data, comp[, c("Analysis", "Model", "deltaAIC", "Rank")])
    }
  }
  
  # Create heatmap
  p <- ggplot(model_data, aes(x = Analysis, y = factor(Rank), fill = deltaAIC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Model), size = 2, color = "black") +  # Reduced from 3 to 2
    scale_fill_viridis(option = "D", direction = -1,
                       name = "ΔAIC",
                       limits = c(0, max(model_data$deltaAIC))) +
    scale_y_discrete(labels = paste("Rank", 1:top_n)) +
    labs(x = "", y = "",
         title = "Model Selection Across Analyses",
         subtitle = paste("Top", top_n, "models per analysis")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10))
  
  return(p)
}

#' Create effect size comparison matrix
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param parameters Vector of parameter names to include
#' @return ggplot object
create_effect_matrix_plot <- function(batch_results, 
                                     parameters = NULL) {
  
  successful <- batch_results$results[sapply(batch_results$results, function(r) r$success)]
  
  # Collect all effects
  all_effects <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    
    if (!is.null(result$effects)) {
      best_model <- result$comparison$comparison$Model[1]
      effects <- result$effects[result$effects$Model == best_model, ]
      effects$Analysis <- name
      all_effects <- rbind(all_effects, effects)
    }
  }
  
  # Filter parameters if specified
  if (!is.null(parameters)) {
    all_effects <- all_effects[all_effects$Parameter %in% parameters, ]
  }
  
  # Create matrix plot
  p <- ggplot(all_effects, aes(x = Analysis, y = Parameter)) +
    geom_tile(aes(fill = log2(OddsRatio)), color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", OddsRatio),
                  color = p_value < 0.05), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         name = "log2(OR)") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60"),
                       guide = "none") +
    labs(x = "", y = "",
         title = "Effect Size Matrix",
         subtitle = "Odds ratios from best models (bold = p < 0.05)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
  return(p)
}

#' Create sample size comparison plot
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @return ggplot object
create_sample_size_plot <- function(batch_results) {
  
  successful <- batch_results$results[sapply(batch_results$results, function(r) r$success)]
  
  # Collect sample sizes
  sample_data <- data.frame(
    Analysis = names(successful),
    N_species = sapply(successful, function(r) r$prepared_data$n_species),
    N_models = sapply(successful, function(r) length(r$models$models)),
    N_converged = sapply(successful, function(r) {
      sum(sapply(r$models$convergence, function(c) c$converged))
    })
  )
  
  # Create long format for plotting
  sample_long <- pivot_longer(sample_data, 
                             cols = c(N_species, N_models, N_converged),
                             names_to = "Metric",
                             values_to = "Count")
  
  # Clean metric names
  sample_long$Metric <- factor(sample_long$Metric,
                               levels = c("N_species", "N_models", "N_converged"),
                               labels = c("Species", "Models Fitted", "Models Converged"))
  
  # Create plot
  p <- ggplot(sample_long, aes(x = Analysis, y = Count, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_viridis_d(option = "C") +
    labs(x = "", y = "Count",
         title = "Sample Sizes and Model Fitting Success",
         fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  return(p)
}

#' Create interaction visualization for a specific analysis
#' 
#' @param analysis_result Single analysis result
#' @param predictor1 First predictor variable name
#' @param predictor2 Second predictor variable name (for interaction)
#' @param n_points Number of points for continuous predictors
#' @return ggplot object
create_interaction_plot <- function(analysis_result, 
                                   predictor1, 
                                   predictor2,
                                   n_points = 100) {
  
  # Get best interaction model
  models <- analysis_result$models$models
  interaction_models <- names(models)[grepl(paste0(predictor1, ".*:", predictor2), 
                                           names(models)) |
                                     grepl(paste0(predictor2, ".*:", predictor1), 
                                           names(models))]
  
  if (length(interaction_models) == 0) {
    stop("No interaction model found for specified predictors")
  }
  
  # Use best interaction model
  model_name <- interaction_models[1]
  model <- models[[model_name]]
  
  # Determine variable types
  var_info <- analysis_result$prepared_data$variable_info
  type1 <- var_info[[predictor1]]$type
  type2 <- var_info[[predictor2]]$type
  
  # Create prediction data based on types
  data <- analysis_result$prepared_data$data
  
  if (type1 == "binary" && type2 == "binary") {
    # 2x2 interaction
    pred_data <- expand.grid(
      var1 = c(0, 1),
      var2 = c(0, 1)
    )
    names(pred_data) <- c(predictor1, predictor2)
    
    # Add control variables at mean
    for (var in analysis_result$config$controls) {
      if (is.numeric(data[[var]])) {
        pred_data[[var]] <- mean(data[[var]], na.rm = TRUE)
      }
    }
    
    # Predict
    pred_data$predicted <- predict(model, newdata = pred_data, type = "response")
    
    # Create plot
    pred_data[[predictor2]] <- factor(pred_data[[predictor2]], 
                                      labels = c("Absent", "Present"))
    
    p <- ggplot(pred_data, aes_string(x = predictor1, y = "predicted", 
                                      color = predictor2, group = predictor2)) +
      geom_line(size = 1.5) +
      geom_point(size = 3) +
      scale_x_continuous(breaks = c(0, 1), labels = c("Absent", "Present")) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = var_info[[predictor1]]$label,
           y = paste("P(", var_info[[analysis_result$config$response]]$label, ")", sep = ""),
           color = var_info[[predictor2]]$label,
           title = "Interaction Effect",
           subtitle = paste("Model:", model_name)) +
      theme_minimal()
    
  } else if (type1 %in% c("continuous", "continuous_special") && type2 == "binary") {
    # Continuous x Binary interaction
    cont_range <- range(data[[predictor1]], na.rm = TRUE)
    
    pred_data <- expand.grid(
      var1 = seq(cont_range[1], cont_range[2], length.out = n_points),
      var2 = c(0, 1)
    )
    names(pred_data) <- c(predictor1, predictor2)
    
    # Add controls
    for (var in analysis_result$config$controls) {
      if (var != predictor1 && is.numeric(data[[var]])) {
        pred_data[[var]] <- mean(data[[var]], na.rm = TRUE)
      }
    }
    
    # Predict
    pred_data$predicted <- predict(model, newdata = pred_data, type = "response")
    
    # Create plot
    pred_data[[predictor2]] <- factor(pred_data[[predictor2]], 
                                      labels = c("Absent", "Present"))
    
    p <- ggplot(pred_data, aes_string(x = predictor1, y = "predicted", 
                                      color = predictor2, fill = predictor2)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = predicted - 0.1, ymax = predicted + 0.1), 
                  alpha = 0.2, color = NA) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = var_info[[predictor1]]$label,
           y = paste("P(", var_info[[analysis_result$config$response]]$label, ")", sep = ""),
           color = var_info[[predictor2]]$label,
           fill = var_info[[predictor2]]$label,
           title = "Interaction Effect") +
      theme_minimal()
    
  } else {
    stop("Interaction plot not implemented for this variable type combination")
  }
  
  return(p)
}

#' Create a comprehensive summary figure for batch results
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param main_parameter Main parameter of interest
#' @param output_file Path to save the figure
#' @return Combined plot object
create_batch_summary_figure <- function(batch_results,
                                       main_parameter,
                                       output_file = NULL) {
  
  # Create individual plots
  p1 <- create_comparative_forest_plot(batch_results, main_parameter)
  p2 <- create_model_selection_heatmap(batch_results, top_n = 3)
  p3 <- create_sample_size_plot(batch_results)
  
  # Try to create effect matrix
  p4 <- tryCatch({
    create_effect_matrix_plot(batch_results)
  }, error = function(e) {
    # Fallback plot if not enough data
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "Insufficient data for effect matrix",
               size = 5) +
      theme_void()
  })
  
  # Combine using patchwork
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "PhyloGLM Batch Analysis Summary",
      subtitle = paste("Comparing", length(batch_results$results), "analyses",
                      if (!is.null(batch_results$bootstrap_n)) 
                        paste("(", batch_results$bootstrap_n, " bootstrap iterations)", sep = "") 
                      else ""),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 14, height = 10, dpi = 300)
  }
  
  return(combined)
}

#' Create model averaging comparison plot
#' 
#' @param batch_results Results from run_phyloglm_batch
#' @param parameter_name Parameter to compare
#' @return ggplot object
create_model_averaging_plot <- function(batch_results, parameter_name) {
  
  successful <- batch_results$results[sapply(batch_results$results, function(r) r$success)]
  
  # Collect averaged effects
  avg_data <- data.frame()
  
  for (name in names(successful)) {
    result <- successful[[name]]
    
    if (!is.null(result$averaged)) {
      param_avg <- result$averaged[result$averaged$Parameter == parameter_name, ]
      
      if (nrow(param_avg) > 0) {
        param_avg$Analysis <- name
        avg_data <- rbind(avg_data, param_avg[1, ])
      }
    }
  }
  
  if (nrow(avg_data) == 0) {
    stop("No model-averaged effects found for parameter:", parameter_name)
  }
  
  # Create plot
  p <- ggplot(avg_data, aes(x = OddsRatio, y = reorder(Analysis, OddsRatio))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), 
                   height = 0.2, size = 0.5) +
    geom_point(aes(size = N_models), color = "darkgreen") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    scale_size_continuous(range = c(2, 6), name = "Models averaged") +
    labs(x = "Model-averaged Odds Ratio (95% CI)",
         y = "",
         title = paste("Model-Averaged Effects:", parameter_name),
         subtitle = "Effects weighted by model support") +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}