# Visualization Functions for Stepwise Model Expansion Results

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Create forest plot for expanded model showing all predictors
#' @param expansion_result Result from stepwise_expand_model
#' @param title Plot title
#' @return ggplot object
plot_expanded_model_forest <- function(expansion_result, 
                                     title = "Expanded Model Effects") {
  
  # Check if we have effects to plot
  if (is.null(expansion_result$effects) || nrow(expansion_result$effects) == 0) {
    # Calculate effects if not present
    effects <- calculate_effect_sizes(
      expansion_result$final_model,
      model_name = "Expanded_Model"
    )
  } else {
    effects <- expansion_result$effects
  }
  
  # Remove intercept
  effects <- effects[effects$Variable != "(Intercept)", ]
  
  # Clean variable names for display
  effects$Variable_Clean <- effects$Variable
  effects$Variable_Clean <- gsub("HighConfidence_Coop", "Cooperative Breeding", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("FemaleSong_Agg01", "Female Song", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("logMass_AVONET", "Body Mass (log)", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("Territory_transformed", "Territory (1-3)", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("Migration_AVONET_transformed", "Migration (1-3)", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("Centroid.Latitude_AVONET_transformed", "Latitude (abs)", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("GeographicRegion_Jetz", "Geographic Region", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("Griesser2017FamilialLiving", "Familial Living", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("logMaleFemalePlumageDiffAbs", "Plumage Dimorphism", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("PercentAbsLogWingDimorphism", "Wing Dimorphism", effects$Variable_Clean)
  effects$Variable_Clean <- gsub(":", " × ", effects$Variable_Clean)
  effects$Variable_Clean <- gsub("Strong", " (Strong)", effects$Variable_Clean)
  
  # Determine if we have bootstrap CIs
  has_ci <- !all(is.na(effects$CI_lower)) && !all(is.na(effects$CI_upper))
  
  if (!has_ci && !all(is.na(effects$Std.Error))) {
    # Calculate CIs from standard errors
    effects$CI_lower <- effects$Estimate - 1.96 * effects$Std.Error
    effects$CI_upper <- effects$Estimate + 1.96 * effects$Std.Error
  }
  
  # Calculate odds ratios
  effects$OR <- exp(effects$Estimate)
  effects$OR_CI_lower <- exp(effects$CI_lower)
  effects$OR_CI_upper <- exp(effects$CI_upper)
  
  # Mark significance
  effects$Significant <- effects$P.Value < 0.05
  
  # Categorize predictors
  effects$Category <- case_when(
    grepl("Cooperative Breeding|Female Song", effects$Variable_Clean) & 
      !grepl("×", effects$Variable_Clean) ~ "Main Effects",
    grepl("Territory|Migration|Region", effects$Variable_Clean) & 
      !grepl("×", effects$Variable_Clean) ~ "Environmental",
    grepl("Mass|Dimorphism|Latitude", effects$Variable_Clean) & 
      !grepl("×", effects$Variable_Clean) ~ "Morphological/Geographic",
    grepl("Familial", effects$Variable_Clean) ~ "Social",
    grepl("×", effects$Variable_Clean) ~ "Interactions",
    TRUE ~ "Other"
  )
  
  # Order by category and effect size
  effects <- effects %>%
    arrange(Category, desc(abs(Estimate)))
  
  # Create plot
  p <- ggplot(effects, aes(x = OR, y = reorder(Variable_Clean, seq_along(Variable_Clean)))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper),
                   height = 0.2, color = "gray30") +
    geom_point(aes(color = Category, shape = Significant), size = 3) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                      labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05")) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title,
      subtitle = paste("Final model with", nrow(effects), "predictors"),
      x = "Odds Ratio (95% CI)",
      y = "",
      shape = "Significance",
      color = "Category"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  return(p)
}

#' Create plot showing model improvement at each step
#' @param expansion_result Result from stepwise_expand_model
#' @return ggplot object
plot_stepwise_improvement <- function(expansion_result) {
  
  # Create data for plotting
  if (length(expansion_result$expansion_history) == 0) {
    # No improvements made
    return(
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = "No predictors improved the model", 
                size = 6) +
        theme_void()
    )
  }
  
  # Build step data
  steps <- data.frame(
    Step = 0,
    Predictor = "Base Model",
    AIC = expansion_result$base_aic,
    Cumulative_Improvement = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(expansion_result$expansion_history)) {
    step_info <- expansion_result$expansion_history[[i]]
    steps <- rbind(steps, data.frame(
      Step = i,
      Predictor = step_info$added,
      AIC = step_info$aic,
      Cumulative_Improvement = expansion_result$base_aic - step_info$aic,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create improvement plot
  p1 <- ggplot(steps, aes(x = Step, y = AIC)) +
    geom_line(color = "darkblue", size = 1) +
    geom_point(size = 3, color = "darkblue") +
    geom_text(aes(label = Predictor), 
             hjust = -0.1, vjust = 0.5, angle = 45, size = 3) +
    scale_x_continuous(breaks = steps$Step) +
    labs(
      title = "Model Improvement by Step",
      x = "Step",
      y = "AIC"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10)
    )
  
  # Create bar plot of improvements
  improvements <- expansion_result$improvement_tests
  improvements <- improvements %>%
    arrange(desc(aic_improvement)) %>%
    mutate(Predictor = factor(predictor, levels = predictor))
  
  p2 <- ggplot(improvements, aes(x = Predictor, y = aic_improvement, fill = kept)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                     labels = c("FALSE" = "Rejected", "TRUE" = "Kept")) +
    coord_flip() +
    labs(
      title = "AIC Improvement by Predictor",
      subtitle = "Red line = threshold for inclusion",
      x = "",
      y = "AIC Improvement",
      fill = "Decision"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
  
  # Combine plots
  combined <- p1 / p2 + plot_layout(heights = c(1, 1.5))
  
  return(combined)
}

#' Create comprehensive summary figure for stepwise results
#' @param expansion_results List of expansion results (e.g., for FS and CB)
#' @param output_file File path to save figure
#' @return Combined plot
create_stepwise_summary_figure <- function(expansion_results, 
                                         output_file = NULL) {
  
  plots <- list()
  
  # Create plots for each analysis
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    
    # Determine direction from name
    if (grepl("^FS_", name)) {
      direction <- "Female Song → Cooperative Breeding"
    } else if (grepl("^CB_", name)) {
      direction <- "Cooperative Breeding → Female Song"
    } else {
      direction <- name
    }
    
    # Forest plot
    forest <- plot_expanded_model_forest(
      result,
      title = paste(direction, "- Expanded Model")
    )
    
    # Improvement plot
    improvement <- plot_stepwise_improvement(result)
    
    plots[[name]] <- list(forest = forest, improvement = improvement)
  }
  
  # Combine all plots
  if (length(plots) == 2) {
    # Two analyses - arrange in 2x2 grid
    combined <- (plots[[1]]$forest | plots[[2]]$forest) / 
                (plots[[1]]$improvement | plots[[2]]$improvement) +
      plot_layout(heights = c(1, 1.5)) +
      plot_annotation(
        title = "Stepwise Model Expansion Results",
        subtitle = "Iterative addition of predictors to best base models",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12)
        )
      )
  } else {
    # Single analysis
    combined <- plots[[1]]$forest / plots[[1]]$improvement +
      plot_layout(heights = c(1, 1.5))
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    # Add extension if not present
    if (!grepl("\\.(png|pdf)$", output_file)) {
      output_file <- paste0(output_file, ".png")
    }
    
    ggsave(output_file, combined, width = 14, height = 12, dpi = 300)
    
    # Also save as PDF
    pdf_file <- gsub("\\.png$", ".pdf", output_file)
    if (pdf_file != output_file) {
      ggsave(pdf_file, combined, width = 14, height = 12, device = "pdf")
    }
  }
  
  return(combined)
}

#' Create comparison table between base and expanded models
#' @param expansion_results List of expansion results
#' @param output_file Optional CSV file to save table
#' @return Data frame with comparison
create_model_comparison_table <- function(expansion_results, 
                                        output_file = NULL) {
  
  comparison <- data.frame()
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    
    # Base model info
    base_info <- data.frame(
      Analysis = name,
      Model = "Base",
      Formula = deparse(result$base_formula),
      N_Predictors = length(all.vars(result$base_formula)) - 1,
      AIC = result$base_aic,
      stringsAsFactors = FALSE
    )
    
    # Expanded model info
    if (result$n_predictors_added > 0) {
      expanded_info <- data.frame(
        Analysis = name,
        Model = "Expanded",
        Formula = deparse(result$final_formula),
        N_Predictors = length(all.vars(result$final_formula)) - 1,
        AIC = result$final_aic,
        stringsAsFactors = FALSE
      )
      
      # Add improvement info
      improvement_info <- data.frame(
        Analysis = name,
        Model = "Improvement",
        Formula = paste("Added:", paste(sapply(result$expansion_history, 
                                              function(x) x$added), 
                                      collapse = ", ")),
        N_Predictors = result$n_predictors_added,
        AIC = -result$total_improvement,  # Negative to show improvement
        stringsAsFactors = FALSE
      )
      
      comparison <- rbind(comparison, base_info, expanded_info, improvement_info)
    } else {
      comparison <- rbind(comparison, base_info)
    }
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    write.csv(comparison, output_file, row.names = FALSE)
  }
  
  return(comparison)
}

#' Plot predictor importance across multiple expanded models
#' @param expansion_results List of expansion results
#' @return ggplot object
plot_predictor_importance <- function(expansion_results) {
  
  # Collect all improvement tests
  all_improvements <- data.frame()
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    improvements <- result$improvement_tests
    improvements$Analysis <- name
    all_improvements <- rbind(all_improvements, improvements)
  }
  
  # Summarize by predictor
  predictor_summary <- all_improvements %>%
    group_by(predictor, description, type) %>%
    summarise(
      mean_improvement = mean(aic_improvement, na.rm = TRUE),
      max_improvement = max(aic_improvement, na.rm = TRUE),
      times_kept = sum(kept, na.rm = TRUE),
      times_tested = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_improvement))
  
  # Create plot
  p <- ggplot(predictor_summary, 
             aes(x = reorder(predictor, mean_improvement), 
                 y = mean_improvement)) +
    geom_bar(stat = "identity", aes(fill = type)) +
    geom_text(aes(label = paste0(times_kept, "/", times_tested)), 
             hjust = -0.2, size = 3) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Predictor Importance Across Analyses",
      subtitle = "Mean AIC improvement when added to base models",
      x = "",
      y = "Mean AIC Improvement",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}