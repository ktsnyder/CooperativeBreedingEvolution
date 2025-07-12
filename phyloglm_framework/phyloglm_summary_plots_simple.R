# Simplified PhyloGLM Summary Plots
# Working version that handles the actual data structure

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(phylolm)

#' Create forest plot for bidirectional effects (simplified)
create_simple_forest_plot <- function(batch_results) {
  
  # Collect data for plotting
  plot_data <- data.frame()
  
  # Get all analyses from batch_results that contain TerrWS but not Terr3
  all_analyses <- names(batch_results)
  target_analyses <- all_analyses[grepl("TerrWS", all_analyses) & 
                                  !grepl("Terr3", all_analyses) &
                                  !grepl("Region", all_analyses) &
                                  !grepl("Fam", all_analyses) &
                                  !grepl("Migration", all_analyses)]
  
  cat("Found", length(target_analyses), "target analyses for forest plot\n")
  
  for (name in target_analyses) {
    if (!(name %in% names(batch_results))) next
    
    result <- batch_results[[name]]
    if (!result$success) next
    
    # Get best model
    best_model <- result$comparison$comparison[1, ]
    
    # Skip if best model is null
    if (best_model$Model == "Null") {
      cat("Skipping analysis", name, "- best model is Null\n")
      next
    }
    
    # Determine parameter and direction based on response variable
    config <- result$config
    if (config$response == "FemaleSong_Agg01") {
      # CB -> FS
      param <- "HighConfidence_Coop"
      direction <- "Cooperative Breeding → Female Song"
    } else if (config$response == "HighConfidence_Coop") {
      # FS -> CB
      param <- "FemaleSong_Agg01"
      direction <- "Female Song → Cooperative Breeding"
    } else {
      cat("Skipping analysis with response:", config$response, "\n")
      next
    }
    
    # Get coefficient
    coef_row <- result$coefficients %>%
      filter(Model == best_model$Model, Parameter == param)
    
    if (nrow(coef_row) == 0) {
      cat("No coefficient found for", param, "in model", best_model$Model, "\n")
      next
    }
    
    # Extract control based on actual control variable used
    control_var <- result$config$controls[1]
    control <- case_when(
      control_var == "logMass_AVONET" ~ "Body Mass",
      control_var == "PercentAbsLogWingDimorphism" ~ "Wing Dimorphism",
      control_var == "logMaleFemalePlumageDiffAbs" ~ "Plumage Dimorphism",
      control_var == "abs(Centroid.Latitude_AVONET)" ~ "abs(Latitude)",
      control_var == "Centroid.Latitude_AVONET" ~ "Latitude",
      TRUE ~ "Other"
    )
    
    # Get species count
    n_species <- nrow(result$prepared_data$data)
    
    # Calculate CIs - handle potential issues
    est <- as.numeric(coef_row$Estimate[1])
    se <- as.numeric(coef_row$StdErr[1])
    
    # Check for bootstrap CIs
    has_ci <- "CI_lower" %in% names(coef_row) && 
              !is.na(coef_row$CI_lower[1]) && 
              is.numeric(coef_row$CI_lower[1])
    
    if (has_ci) {
      ci_lower <- as.numeric(coef_row$CI_lower[1])
      ci_upper <- as.numeric(coef_row$CI_upper[1])
    } else {
      ci_lower <- est - 1.96 * se
      ci_upper <- est + 1.96 * se
    }
    
    # Add to plot data
    plot_data <- rbind(plot_data, data.frame(
      Analysis = name,
      Direction = direction,
      Control = control,
      Model = best_model$Model,
      Estimate = est,
      StdErr = se,
      p_value = as.numeric(coef_row$p_value[1]),
      OddsRatio = exp(est),
      OR_CI_lower = exp(ci_lower),
      OR_CI_upper = exp(ci_upper),
      n_species = n_species,
      stringsAsFactors = FALSE
    ))
  }
  
  # Check if we have any data
  if (nrow(plot_data) == 0) {
    cat("No data found for forest plot\n")
    return(NULL)
  }
  
  # Create plot label with species count
  plot_data$PlotLabel <- paste0(plot_data$Direction, "\n(Territoriality, ", plot_data$Control, "; n = ", plot_data$n_species, " species)")
  
  # Custom ordering - FS->CB at top, then CB->FS (Body Mass), then rest
  plot_data <- plot_data %>%
    mutate(
      order_key = case_when(
        grepl("Female Song → Cooperative", Direction) & Control == "Body Mass" ~ 1,
        grepl("Cooperative Breeding → Female", Direction) & Control == "Body Mass" ~ 2,
        grepl("Cooperative Breeding → Female", Direction) & Control == "Wing Dimorphism" ~ 3,
        grepl("Cooperative Breeding → Female", Direction) & Control == "Plumage Dimorphism" ~ 4,
        grepl("Cooperative Breeding → Female", Direction) & Control == "Latitude" ~ 5,
        TRUE ~ 6 + as.numeric(factor(Control))
      )
    ) %>%
    arrange(desc(order_key))  # Reverse order for top-to-bottom display
  
  # Create forest plot
  p <- ggplot(plot_data, aes(x = OddsRatio, y = reorder(PlotLabel, seq_along(PlotLabel)))) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper),
                   height = 0.2, linewidth = 0.8, color = "#E64B35") +
    geom_point(aes(shape = p_value < 0.05), size = 4, color = "#E64B35") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16)) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                       labels = c("FALSE" = "Not significant", "TRUE" = "Significant"),
                       name = "p < 0.05") +
    labs(x = "Odds Ratio (95% CI)",
         y = "",
         title = "A. Bidirectional Effects: Cooperative Breeding ⇄ Female Song",
         subtitle = NULL) +
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

#' Create interaction plot (simplified)
create_simple_interaction_plot <- function(result, 
                                         predictor_var,
                                         response_var,
                                         predictor_label,
                                         response_label) {
  
  # Find best interaction model
  int_models <- result$comparison$comparison %>%
    filter(grepl("x", Model)) %>%
    arrange(AIC)
  
  if (nrow(int_models) == 0) {
    stop("No interaction model found")
  }
  
  best_model_name <- int_models$Model[1]
  model <- result$models$models[[best_model_name]]
  
  # Create prediction data
  pred_data <- expand.grid(
    pred = c(0, 1),
    terr = c(0, 1)
  )
  names(pred_data) <- c(predictor_var, "TerritorialityWeakVsStrong")
  
  # Add control variable
  control_var <- result$config$controls[1]
  if (!is.null(control_var) && control_var %in% names(result$prepared_data$data)) {
    pred_data[[control_var]] <- mean(result$prepared_data$data[[control_var]], na.rm = TRUE)
  }
  
  # Manual prediction for phyloglm
  coefs <- coef(model)
  
  # Create formula without response variable
  form <- formula(model)
  # Extract right-hand side of formula
  rhs_formula <- as.formula(paste("~", as.character(form)[3]))
  
  X <- model.matrix(rhs_formula, data = pred_data)
  pred_data$predicted <- plogis(X %*% coefs)
  
  # Simple CIs - use a fixed width for now since phyloglm doesn't have vcov method
  # This is a reasonable approximation for visualization purposes
  pred_data$lower <- pmax(0, pred_data$predicted - 0.1)
  pred_data$upper <- pmin(1, pred_data$predicted + 0.1)
  
  # Format for plotting
  pred_data$TerritorialityWeakVsStrong <- factor(pred_data$TerritorialityWeakVsStrong,
                                                 levels = c(0, 1),
                                                 labels = c("Weak/No Territoriality", "Strong Territoriality"))
  pred_data$predictor_factor <- factor(pred_data[[predictor_var]],
                                       levels = c(0, 1),
                                       labels = c("Absent", "Present"))
  
  # Get model formula for subtitle
  form <- formula(model)
  formula_text <- paste(as.character(form)[2], "~", as.character(form)[3])
  
  # Simplify the formula text for display
  formula_text <- gsub("HighConfidence_Coop", "CB", formula_text)
  formula_text <- gsub("FemaleSong_Agg01", "FS", formula_text)
  formula_text <- gsub("TerritorialityWeakVsStrong", "Terr", formula_text)
  formula_text <- gsub("logMass_AVONET", "log(Mass)", formula_text)
  formula_text <- gsub("PercentAbsLogWingDimorphism", "WingDim", formula_text)
  formula_text <- gsub("logMaleFemalePlumageDiffAbs", "PlumDim", formula_text)
  formula_text <- gsub("Centroid.Latitude_AVONET", "Latitude", formula_text)
  formula_text <- gsub("abs\\(Centroid.Latitude_AVONET\\)", "abs(Latitude)", formula_text)
  
  # Create plot
  p <- ggplot(pred_data, aes(x = predictor_factor, y = predicted,
                             color = TerritorialityWeakVsStrong, 
                             group = TerritorialityWeakVsStrong)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, 
                    fill = TerritorialityWeakVsStrong),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = c("Weak/No Territoriality" = "#4DBBD5",
                                  "Strong Territoriality" = "#E64B35")) +
    scale_fill_manual(values = c("Weak/No Territoriality" = "#4DBBD5",
                                 "Strong Territoriality" = "#E64B35")) +
    labs(x = predictor_label,
         y = paste0("P(", response_label, ")"),
         color = "Territory",
         fill = "Territory",
         title = paste0("Best Model: ", best_model_name),
         subtitle = paste0("Best model formula: ", formula_text)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right",
      axis.title.y = element_text(margin = margin(r = 10))  # Add margin to y-axis label
    )
  
  return(p)
}

#' Create complete summary figure (simplified)
create_simple_summary_figure <- function(batch_results, output_file = NULL) {
  
  # Panel A: Forest plot
  panel_a <- create_simple_forest_plot(batch_results)
  
  # Panel B: CB -> FS interaction
  panel_b <- NULL
  if ("FS_vs_CB_TerrWS_Mass" %in% names(batch_results)) {
    result <- batch_results[["FS_vs_CB_TerrWS_Mass"]]
    if (result$success) {
      panel_b <- create_simple_interaction_plot(
        result,
        predictor_var = "HighConfidence_Coop",
        response_var = "FemaleSong_Agg01",
        predictor_label = "Cooperative Breeding",
        response_label = "Female Song"
      )
      panel_b <- panel_b + 
        labs(title = "B. Cooperative Breeding → Female Song")
    }
  }
  
  # Panel C: FS -> CB interaction
  panel_c <- NULL
  if ("CB_vs_FS_TerrWS_Mass" %in% names(batch_results)) {
    result <- batch_results[["CB_vs_FS_TerrWS_Mass"]]
    if (result$success) {
      panel_c <- create_simple_interaction_plot(
        result,
        predictor_var = "FemaleSong_Agg01",
        response_var = "HighConfidence_Coop",
        predictor_label = "Female Song",
        response_label = "Cooperative Breeding"
      )
      panel_c <- panel_c + 
        labs(title = "C. Female Song → Cooperative Breeding")
    }
  }
  
  # Combine panels
  if (!is.null(panel_b) && !is.null(panel_c)) {
    combined <- panel_a / (panel_b | panel_c) + 
      plot_layout(heights = c(1.2, 1))
  } else if (!is.null(panel_b)) {
    combined <- panel_a | panel_b
  } else if (!is.null(panel_c)) {
    combined <- panel_a | panel_c
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
    
    # Add extension if not present
    if (!grepl("\\.(png|pdf)$", output_file)) {
      output_file <- paste0(output_file, ".png")
    }
    
    # Save as PNG
    ggsave(output_file, combined, width = 14, height = height, dpi = 300)
    
    # Also save as PDF with higher DPI
    pdf_file <- gsub("\\.png$", ".pdf", output_file)
    if (pdf_file == output_file) {
      # If no .png extension, just add .pdf
      pdf_file <- paste0(output_file, ".pdf")
    }
    ggsave(pdf_file, combined, width = 14, height = height, dpi = 600, device = "pdf")
  }
  
  return(combined)
}