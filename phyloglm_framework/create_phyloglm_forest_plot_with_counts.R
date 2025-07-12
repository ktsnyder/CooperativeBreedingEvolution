# Forest plot function for phyloglm combined results with species counts
library(ggplot2)
library(dplyr)

#' Create a forest plot from combined phyloglm results
#' 
#' @param all_results Combined results from phyloglm analyses (e.g., from combined_CB_response/all_results.rds)
#' @param predictor_name Name of the predictor to plot (e.g., "FemaleSong_Agg01")
#' @param predictor_label Label for the predictor in the plot (e.g., "Female Song")
#' @param response_label Label for the response in the plot (e.g., "Cooperative Breeding")
#' @param reference_value Optional reference value to show as vertical line
#' @param output_file Path to save the plot
#' 
create_phyloglm_forest_plot_with_counts <- function(all_results,
                                                   predictor_name = "FemaleSong_Agg01",
                                                   predictor_label = "Female Song",
                                                   response_label = "Cooperative Breeding",
                                                   reference_value = NULL,
                                                   output_file = NULL) {
  
  # Extract effects for the specified predictor from each analysis
  effect_data <- data.frame()
  
  for (analysis_name in names(all_results)) {
    result <- all_results[[analysis_name]]
    
    # Skip failed analyses
    if (!result$success || is.null(result$effects)) next
    
    # Get effects for the predictor from the best model
    best_model <- result$comparison$comparison$Model[1]
    predictor_effects <- result$effects %>%
      filter(Model == best_model, 
             Parameter == predictor_name)
    
    if (nrow(predictor_effects) > 0) {
      # Extract sample size from prepared_data
      n_species <- result$prepared_data$n_species
      
      # Extract control/subset information from analysis name
      subset_info <- case_when(
        grepl("Jackknife", analysis_name) ~ "Jackknife by species",
        grepl("Holarctic.*Non.*coop", analysis_name) ~ "Holarctic Non-cooperative",
        grepl("Tropical.*Coop", analysis_name) ~ "Tropical Cooperative",
        grepl("Global.*Coop", analysis_name) ~ "Global Cooperative",
        grepl("Strong.*Terr", analysis_name) ~ "Strong Territoriality",
        grepl("Year.*round.*Terr", analysis_name) ~ "Year-round Territoriality",
        grepl("Wing.*Dim", analysis_name) ~ "Wing Dimorphism",
        grepl("Plum.*Dim", analysis_name) ~ "Plumage Dichromatism",
        grepl("Mass", analysis_name) ~ "Controlling for Mass",
        grepl("absLat", analysis_name) ~ "Controlling for Latitude",
        grepl("Region", analysis_name) ~ "Controlling for Region",
        grepl("Migration", analysis_name) ~ "Controlling for Migration",
        grepl("Fam", analysis_name) ~ "Controlling for Family",
        TRUE ~ analysis_name
      )
      
      # Calculate CIs if not available
      ci_lower <- predictor_effects$CI_lower[1]
      ci_upper <- predictor_effects$CI_upper[1]
      
      if (is.na(ci_lower) || is.na(ci_upper)) {
        # Calculate from standard error
        ci_lower <- predictor_effects$Estimate[1] - 1.96 * predictor_effects$StdErr[1]
        ci_upper <- predictor_effects$Estimate[1] + 1.96 * predictor_effects$StdErr[1]
      }
      
      effect_data <- rbind(effect_data, data.frame(
        Analysis = analysis_name,
        Subset = subset_info,
        Estimate = predictor_effects$Estimate[1],
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        OddsRatio = predictor_effects$OddsRatio[1],
        OR_CI_Lower = predictor_effects$OR_CI_lower[1],
        OR_CI_Upper = predictor_effects$OR_CI_upper[1],
        p_value = predictor_effects$p_value[1],
        Significance = predictor_effects$Significance[1],
        N_Species = n_species,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add species counts to labels
  effect_data$Subset_Label <- paste0(effect_data$Subset, "\n(n = ", effect_data$N_Species, " species)")
  
  # Define the desired order (from bottom to top in the plot)
  subset_order <- c("Jackknife by species", "Holarctic Non-cooperative", "Tropical Cooperative", 
                    "Global Cooperative", "Strong Territoriality", "Year-round Territoriality", 
                    "Wing Dimorphism", "Plumage Dichromatism", "Controlling for Mass",
                    "Controlling for Latitude", "Controlling for Region", 
                    "Controlling for Migration", "Controlling for Family")
  
  # Create ordered factor
  effect_data$Subset_Label <- factor(effect_data$Subset_Label,
                                    levels = effect_data$Subset_Label[match(subset_order, effect_data$Subset)])
  
  # Remove NA levels
  effect_data <- effect_data[!is.na(effect_data$Subset_Label), ]
  
  # Create the forest plot
  p <- ggplot(effect_data, aes(x = Estimate, y = Subset_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                   height = 0.2, linewidth = 1) +
    geom_point(aes(color = p_value < 0.05), size = 4) +
    geom_text(aes(label = sprintf("%.3f", Estimate)), 
              vjust = -1.2, size = 3.5) +
    scale_color_manual(values = c("TRUE" = "#1B4F72", "FALSE" = "gray50"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "Significance") +
    labs(
      title = paste0(predictor_label, " → ", response_label, " Path Coefficients"),
      subtitle = "Coefficient estimates and 95% CI from best models",
      x = "Coefficient Estimate",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 20, 10, 10),
      legend.position = "bottom"
    )
  
  # Add reference line if provided
  if (!is.null(reference_value)) {
    p <- p + geom_vline(xintercept = reference_value, 
                        linetype = "dotted", 
                        color = "#2E86AB", 
                        linewidth = 1.2)
  }
  
  # Save if output file specified
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    message("Plot saved to: ", output_file)
  }
  
  return(p)
}

# Function to create forest plots for all predictors in the results
create_all_phyloglm_forest_plots <- function(results_file,
                                           output_dir = "Outputs/Figures/phyloglm_forest_plots") {
  
  # Load results
  all_results <- readRDS(results_file)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get response variable name from first successful analysis
  first_success <- NULL
  for (r in all_results) {
    if (r$success) {
      first_success <- r
      break
    }
  }
  
  if (is.null(first_success)) {
    stop("No successful analyses found in results")
  }
  
  response_var <- first_success$config$response
  response_label <- case_when(
    response_var == "HighConfidence_Coop" ~ "Cooperative Breeding",
    response_var == "FemaleSong_Agg01" ~ "Female Song",
    TRUE ~ response_var
  )
  
  # Get unique predictors across all analyses
  all_predictors <- unique(unlist(lapply(all_results, function(r) {
    if (r$success) r$config$predictors
  })))
  
  # Create a plot for each predictor
  plots <- list()
  
  for (pred in all_predictors) {
    pred_label <- case_when(
      pred == "FemaleSong_Agg01" ~ "Female Song",
      pred == "HighConfidence_Coop" ~ "Cooperative Breeding",
      pred == "TerritorialityWeakVsStrong" ~ "Territoriality",
      pred == "logMass_AVONET" ~ "Body Mass (log)",
      TRUE ~ pred
    )
    
    # Skip if this is the response variable
    if (pred == response_var) next
    
    cat("Creating forest plot for:", pred_label, "→", response_label, "\n")
    
    tryCatch({
      p <- create_phyloglm_forest_plot_with_counts(
        all_results = all_results,
        predictor_name = pred,
        predictor_label = pred_label,
        response_label = response_label,
        output_file = file.path(output_dir, paste0("forest_plot_", 
                                                   gsub(" ", "_", pred_label), "_to_",
                                                   gsub(" ", "_", response_label), ".png"))
      )
      plots[[pred]] <- p
    }, error = function(e) {
      cat("Error creating plot for", pred, ":", e$message, "\n")
    })
  }
  
  return(plots)
}