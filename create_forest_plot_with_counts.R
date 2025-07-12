# Required libraries
library(ggplot2)
library(dplyr)

# Note: This function requires the extract_phylopath_results() function
# from create_phylopath_bias_robustness_figure.R
# Source that file first or ensure the function is available

#' Create an improved forest plot with species counts for phylopath downsampling
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param rate_to_plot Character string specifying which rate to plot. Options include:
#'   "COOP->FS" (default), "FS->COOP", "TERR->FS", "TERR->COOP", "MASS->FS", "MASS->COOP", "MASS->TERR"
#' @param reference_value Reference coefficient value from full dataset. If NULL (default), 
#'   will be calculated from phylopath_full_dataset_result.rds if available
#' @param full_dataset_n Number of species in the full dataset. If NULL (default),
#'   will be extracted from phylopath_full_dataset_result.rds if available
#' @param output_file Path for saving the figure. If NULL (default), will be automatically
#'   generated based on rate_to_plot (e.g., "Outputs/Figures/forest_plot_COOP_FS_with_counts.pdf")


create_forest_plot_with_counts <- function(bias_results_list = NULL,
                                           rate_to_plot = "COOP->FS",
                                           reference_value = NULL,
                                           full_dataset_n = NULL,
                                           output_file = NULL) {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results(results_dir = "Outputs")
  }
  
  # Load full dataset reference value if not provided
  if (is.null(reference_value)) {
    if (file.exists("phylopath_full_dataset_result.rds")) {
      full_result <- readRDS("phylopath_full_dataset_result.rds")
      
      # Calculate conditional average for the selected rate
      if (!is.null(full_result$result)) {
        # Get conditional average using phylopath
        avg_result <- phylopath::average(full_result$result, avg_method = "conditional")
        
        # Extract the coefficient based on the selected rate
        rate_mapping <- list(
          "COOP->FS" = c(from = "HighConfidence_Coop", to = "FemaleSong_Agg01"),
          "FS->COOP" = c(from = "FemaleSong_Agg01", to = "HighConfidence_Coop"),
          "TERR->FS" = c(from = "TerritorialityWeakVsStrong", to = "FemaleSong_Agg01"),
          "TERR->COOP" = c(from = "TerritorialityWeakVsStrong", to = "HighConfidence_Coop"),
          "MASS->FS" = c(from = "logMass_AVONET", to = "FemaleSong_Agg01"),
          "MASS->COOP" = c(from = "logMass_AVONET", to = "HighConfidence_Coop"),
          "MASS->TERR" = c(from = "logMass_AVONET", to = "TerritorialityWeakVsStrong")
        )
        
        if (rate_to_plot %in% names(rate_mapping)) {
          selected_rate <- rate_mapping[[rate_to_plot]]
          reference_value <- avg_result$coef[selected_rate["from"], selected_rate["to"]]
          
          # Also get the sample size if not provided
          if (is.null(full_dataset_n) && !is.null(full_result$n_species)) {
            full_dataset_n <- full_result$n_species
          }
        }
      }
    } else {
      warning("phylopath_full_dataset_result.rds not found. Using default reference value.")
      reference_value <- 0.556  # Default fallback
    }
  }
  
  # Set default full_dataset_n if still NULL
  if (is.null(full_dataset_n)) {
    full_dataset_n <- 875  # Default fallback
  }
  
  # Generate output file path if not provided
  if (is.null(output_file)) {
    # Convert rate_to_plot to safe filename format
    safe_rate_name <- gsub("->", "_", rate_to_plot)
    output_file <- paste0("Outputs/Figures/forest_plot_", safe_rate_name, "_with_counts.pdf")
  }
  
  # Extract coefficients and species counts
  coef_summary <- data.frame(
    Bias_Correction = character(),
    Mean_Coefficient = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    N_Iterations = integer(),
    N_Species = integer(),
    stringsAsFactors = FALSE
  )
  
  # Define mapping for rate patterns
  rate_patterns <- list(
    "COOP->FS" = c(from = "Coop", to = "FemaleSong"),
    "FS->COOP" = c(from = "FemaleSong", to = "Coop"),
    "TERR->FS" = c(from = "Terr", to = "FemaleSong"),
    "TERR->COOP" = c(from = "Terr", to = "Coop"),
    "MASS->FS" = c(from = "Mass", to = "FemaleSong"),
    "MASS->COOP" = c(from = "Mass", to = "Coop"),
    "MASS->TERR" = c(from = "Mass", to = "Terr")
  )
  
  # Validate rate_to_plot
  if (!rate_to_plot %in% names(rate_patterns)) {
    stop("Invalid rate_to_plot. Must be one of: ", paste(names(rate_patterns), collapse = ", "))
  }
  
  # Get the pattern for the selected rate
  selected_pattern <- rate_patterns[[rate_to_plot]]
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      # Find coefficient column for selected rate
      # Look for columns that match the pattern: from.*to.*to.*est
      pattern_str <- paste0(selected_pattern["from"], ".*to.*", selected_pattern["to"], ".*est")
      coef_columns <- grep(pattern_str, names(result$detailed_models), value = TRUE, ignore.case = TRUE)
      
      if (length(coef_columns) > 0) {
        coef_col <- coef_columns[1]
        
        # Extract coefficients - calculate conditional average for each seed
        coef_data <- result$detailed_models %>%
          filter(delta_CICc < 2) %>%
          group_by(seed) %>%
          summarise(
            coef = {
              # Get data for this seed (all models with delta_CICc < 2)
              seed_data <- pick(everything())
              
              # Calculate CICc weights for ALL models with delta_CICc < 2
              all_weights <- exp(-0.5 * seed_data$delta_CICc)
              all_weights <- all_weights / sum(all_weights)
              
              # Now filter to rows with non-NA coefficients
              valid_rows <- !is.na(seed_data[[coef_col]])
              
              if (sum(valid_rows) > 0) {
                # Get weights for models that have the coefficient
                valid_weights <- all_weights[valid_rows]
                # Renormalize these weights to sum to 1
                valid_weights <- valid_weights / sum(valid_weights)
                # Calculate conditional average
                sum(seed_data[[coef_col]][valid_rows] * valid_weights)
              } else {
                NA_real_
              }
            },
            n_species = first(nSpecies),
            .groups = "drop"
          ) %>%
          filter(!is.na(coef))
        
        if (nrow(coef_data) > 0) {
          coef_summary <- rbind(coef_summary, data.frame(
            Bias_Correction = bias_name,
            Mean_Coefficient = mean(coef_data$coef),
            CI_Lower = quantile(coef_data$coef, 0.025),
            CI_Upper = quantile(coef_data$coef, 0.975),
            N_Iterations = length(unique(coef_data$seed)),
            N_Species = round(mean(coef_data$n_species, na.rm = TRUE)),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Add species counts to the bias correction labels
  coef_summary$Bias_Label <- paste0(coef_summary$Bias_Correction, "\n(n ~ ", coef_summary$N_Species, " species)")
  
  # Define the desired order (from bottom to top in the plot)
  bias_order <- c("Jackknife by species", "Holarctic Non-cooperative", "Tropical Cooperative", "Global Cooperative",
                  "Strong Territoriality", "Year-round Territoriality", 
                  "Wing Dimorphism", "Plumage Dichromatism")
  
  # Create ordered factor for plotting
  coef_summary$Bias_Label <- factor(coef_summary$Bias_Label,
                                   levels = coef_summary$Bias_Label[match(bias_order, coef_summary$Bias_Correction)])
  
  # Create the forest plot
  forest_plot <- ggplot(coef_summary, 
                        aes(x = Mean_Coefficient, 
                            y = Bias_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = reference_value, linetype = "dotted", 
               color = "#2E86AB", linewidth = 1.2) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                   height = 0.2, linewidth = 1) +
    geom_point(size = 4, color = "#1B4F72") +
    geom_text(aes(label = sprintf("%.3f", Mean_Coefficient)), 
              vjust = -1.2, size = 3.5) +
    scale_x_continuous(
      limits = c(
        min(0, min(coef_summary$CI_Lower) - 0.1),
        max(0.8, max(coef_summary$CI_Upper) + 0.1)
      ),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    labs(
      title = paste0(gsub("->", " → ", rate_to_plot), " Path Coefficients"),
      subtitle = "Mean and 95% CI of conditional averages from 500 downsampling iterations",
      x = "Path Coefficient",
      y = NULL
    ) +
    annotate("text", x = reference_value + 0.02, y = 0.5, 
             label = paste0("Full dataset (n = ", full_dataset_n, " species)"), 
             color = "#2E86AB", size = 3, hjust = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 20, 10, 10)
    )
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, forest_plot, width = 8, height = 6, dpi = 600)
    ggsave(gsub(".pdf", ".png", output_file), forest_plot, 
           width = 8, height = 6, dpi = 300)
  }
  
  return(forest_plot)
}