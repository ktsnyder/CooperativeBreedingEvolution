# Clean version of PhyloGLM Summary Plots
# Simplified to work with the actual data structure

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Create clean forest plot for bidirectional effects
create_clean_forest_plot <- function(batch_results,
                                    predictor_var = "HighConfidence_Coop",
                                    response_var = "FemaleSong_Agg01",
                                    predictor_label = "Cooperative Breeding",
                                    response_label = "Female Song") {
  
  # Collect relevant effects
  plot_data <- data.frame()
  
  # Process each analysis
  for (name in names(batch_results$results)) {
    # Skip unwanted analyses
    if (grepl("Terr3", name)) next
    if (grepl("Migration|Region", name)) next
    
    result <- batch_results$results[[name]]
    if (!result$success) next
    
    # Check if this analysis has our variables
    config <- result$config
    has_forward <- config$response == response_var && predictor_var %in% config$predictors
    has_reverse <- config$response == predictor_var && response_var %in% config$predictors
    
    if (!has_forward && !has_reverse) next
    
    # Must have territoriality
    if (!"TerritorialityWeakVsStrong" %in% config$predictors) next
    
    # Get best model
    best_model <- result$comparison$comparison[1, ]
    
    # Get coefficient for the relevant parameter
    param_name <- if (has_forward) predictor_var else response_var
    coef_data <- result$coefficients %>%
      filter(Model == best_model$Model, 
             Parameter == param_name)
    
    if (nrow(coef_data) == 0) next
    
    # Extract control variable
    control <- if ("logMass_AVONET" %in% config$controls) {
      "Body Mass"
    } else if ("PercentAbsLogWingDimorphism" %in% config$controls) {
      "Wing Dimorphism"  
    } else if ("logMaleFemalePlumageDiffAbs" %in% config$controls) {
      "Plumage Dimorphism"
    } else {
      "Other"
    }
    
    # Create plot entry
    entry <- data.frame(
      Analysis = name,
      Direction = if (has_forward) {
        paste0(predictor_label, " → ", response_label)
      } else {
        paste0(response_label, " → ", predictor_label)
      },
      Control = control,
      Model = best_model$Model,
      Estimate = coef_data$Estimate[1],
      StdErr = coef_data$StdErr[1],
      p_value = coef_data$p_value[1],
      OddsRatio = exp(coef_data$Estimate[1]),
      CI_lower = if ("CI_lower" %in% names(coef_data) && length(coef_data$CI_lower) > 0 && !is.na(coef_data$CI_lower[1])) {
        coef_data$CI_lower[1]
      } else {
        coef_data$Estimate[1] - 1.96 * coef_data$StdErr[1]
      },
      CI_upper = if ("CI_upper" %in% names(coef_data) && length(coef_data$CI_upper) > 0 && !is.na(coef_data$CI_upper[1])) {
        coef_data$CI_upper[1]
      } else {
        coef_data$Estimate[1] + 1.96 * coef_data$StdErr[1]
      }
    )
    
    plot_data <- rbind(plot_data, entry)
  }
  
  # Calculate OR CIs
  plot_data$OR_CI_lower <- exp(plot_data$CI_lower)
  plot_data$OR_CI_upper <- exp(plot_data$CI_upper)
  
  # Create plot labels
  plot_data$PlotLabel <- paste0(plot_data$Direction, "\n(", plot_data$Control, ")")
  
  # Order by direction then control
  plot_data <- plot_data %>%
    arrange(Direction, Control)
  
  # Create forest plot
  p <- ggplot(plot_data, aes(x = OddsRatio, y = factor(PlotLabel, levels = rev(unique(PlotLabel))))) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper),
                   height = 0.2, linewidth = 0.8, color = "#E64B35") +
    geom_point(aes(shape = p_value < 0.05), size = 4, color = "#E64B35") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16)) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("Significant", "Not significant"),
                       name = "") +
    labs(x = "Odds Ratio (95% CI)",
         y = "",
         title = paste0("A. Bidirectional Effects: ", predictor_label, " ⇄ ", response_label),
         subtitle = "Best models with territorial moderation") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

#' Create interaction plot (fixed for phyloglm)
create_clean_interaction_plot <- function(result,
                                         predictor_var,
                                         terr_var = "TerritorialityWeakVsStrong",
                                         response_var,
                                         predictor_label,
                                         response_label) {
  
  # Find best interaction model
  interaction_models <- result$comparison$comparison %>%
    filter(grepl("x", Model)) %>%
    arrange(AIC)
  
  if (nrow(interaction_models) == 0) {
    stop("No interaction model found")
  }
  
  best_model_name <- interaction_models$Model[1]
  model <- result$models$models[[best_model_name]]
  
  # Create prediction grid
  pred_data <- expand.grid(
    pred = c(0, 1),
    terr = c(0, 1)
  )
  names(pred_data) <- c(predictor_var, terr_var)
  
  # Add control variable
  data <- result$prepared_data$data
  control_var <- result$config$controls[1]
  if (!is.null(control_var) && control_var %in% names(data)) {
    pred_data[[control_var]] <- mean(data[[control_var]], na.rm = TRUE)
  }
  
  # Manual prediction for phyloglm
  coefs <- coef(model)
  X <- model.matrix(formula(model), data = pred_data)
  pred_data$predicted <- plogis(X %*% coefs)
  
  # Simple confidence intervals
  vcov_mat <- vcov(model)
  pred_se <- sqrt(diag(X %*% vcov_mat %*% t(X)))
  pred_data$lower <- plogis(qlogis(pred_data$predicted) - 1.96 * pred_se)
  pred_data$upper <- plogis(qlogis(pred_data$predicted) + 1.96 * pred_se)
  
  # Format for plotting
  pred_data[[terr_var]] <- factor(pred_data[[terr_var]],
                                  levels = c(0, 1),
                                  labels = c("Weak/No Territory", "Strong Territory"))
  pred_data[[predictor_var]] <- factor(pred_data[[predictor_var]],
                                       levels = c(0, 1),
                                       labels = c("Absent", "Present"))
  
  # Get interaction statistics
  coef_summary <- summary(model)$coefficients
  int_term <- paste0(predictor_var, ":", terr_var)
  if (int_term %in% rownames(coef_summary)) {
    int_p <- coef_summary[int_term, "p.value"]
    int_text <- paste0("Interaction p = ", format.pval(int_p, digits = 3))
  } else {
    int_text <- "No interaction term"
  }
  
  # Create plot
  p <- ggplot(pred_data, aes_string(x = predictor_var, y = "predicted",
                                    color = terr_var, group = terr_var)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[terr_var]]),
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
         title = paste0("Best Model: ", best_model_name),
         subtitle = int_text) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right"
    )
  
  return(p)
}

#' Create the complete summary figure
create_complete_summary_figure <- function(batch_results,
                                          predictor_var = "HighConfidence_Coop",
                                          response_var = "FemaleSong_Agg01",
                                          predictor_label = "Cooperative Breeding",
                                          response_label = "Female Song",
                                          output_file = NULL) {
  
  # Panel A: Forest plot
  panel_a <- create_clean_forest_plot(
    batch_results,
    predictor_var = predictor_var,
    response_var = response_var,
    predictor_label = predictor_label,
    response_label = response_label
  )
  
  # Panel B: Forward interaction (CB -> FS)
  forward_analysis <- "FS_CB_TerrWS_Mass"  # Use the main analysis
  panel_b <- NULL
  
  if (forward_analysis %in% names(batch_results$results)) {
    result <- batch_results$results[[forward_analysis]]
    if (result$success) {
      panel_b <- create_clean_interaction_plot(
        result,
        predictor_var = predictor_var,
        terr_var = "TerritorialityWeakVsStrong",
        response_var = response_var,
        predictor_label = predictor_label,
        response_label = response_label
      )
      panel_b <- panel_b + 
        labs(title = paste0("B. ", predictor_label, " → ", response_label))
    }
  }
  
  # Panel C: Reverse interaction (FS -> CB)
  reverse_analysis <- "CB_FS_TerrWS_Mass"
  panel_c <- NULL
  
  if (reverse_analysis %in% names(batch_results$results)) {
    result <- batch_results$results[[reverse_analysis]]
    if (result$success) {
      panel_c <- create_clean_interaction_plot(
        result,
        predictor_var = response_var,  # Swap
        terr_var = "TerritorialityWeakVsStrong",
        response_var = predictor_var,  # Swap
        predictor_label = response_label,  # Swap
        response_label = predictor_label   # Swap
      )
      panel_c <- panel_c + 
        labs(title = paste0("C. ", response_label, " → ", predictor_label))
    }
  }
  
  # Combine panels
  if (!is.null(panel_b) && !is.null(panel_c)) {
    combined <- panel_a / (panel_b | panel_c) + 
      plot_layout(heights = c(1.2, 1))
  } else if (!is.null(panel_b)) {
    combined <- panel_a | panel_b
  } else {
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
  
  if (!is.null(output_file)) {
    height <- if (!is.null(panel_b) && !is.null(panel_c)) 10 else 6
    ggsave(output_file, combined, width = 14, height = height, dpi = 300)
  }
  
  return(combined)
}