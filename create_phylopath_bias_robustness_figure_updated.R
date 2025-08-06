# create_bias_robustness_figure.R
# Creates a multi-panel figure showing robustness of CB-FS association to various biases
# 7/8/2025 - search for detailed_models csv outputs recursively in extract_phylopath_results - HOWEVER, create_bias_robustness_figure() is not the most up-to-date. In fact, this script might be mostly obsolete since the functions I actually use have been moved to their own script files (create_enhanced_DAG.R, create_bias_model_consistency_heatmap.R). Yes, it seems it might be obsolete except for extract_phylopath_results().

library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(tidyr)

#' Extract results from phylopath downsampling outputs
#'
#' @param results_dir Directory containing phylopath downsampling results
#' @param trait_set Character string specifying which trait combination to extract
#'   Default is "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
#' @param n_iterations Optional filter for number of iterations (e.g., 500, 10)
#' @param verbose Logical - whether to print debug information about file searches
#' @return List of data frames with processed results for each bias correction
extract_phylopath_results <- function(results_dir = "Outputs", 
                                      trait_set = NULL,
                                      n_iterations = NULL,
                                      verbose = FALSE) {
  
  # Set default trait set for backward compatibility
  if (is.null(trait_set)) {
    trait_set <- "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
  }
  
  n_iterations_input = n_iterations # because we need to have n_iterations be NULL for JackknifeSpecies but then reset n_iterations back to the input for all the other ones
  
  # Define the bias corrections we're looking for
  bias_types <- list(
    "Full Dataset" = NULL,  # We'll need to handle this separately
    "Holarctic Non-cooperative" = "Remove83HolarcticNoncoop",
    "Tropical Cooperative" = "Remove24TropicalCoop", 
    "Global Cooperative" = "Remove15GlobalCoop",
    "Strong Territoriality" = "Remove266StrongTerr",
    "Year-round Territoriality" = "Remove155Terr3",
    "Jackknife by species" = "JackknifeSpecies",
    "Wing Dimorphism" = "Remove71HighWingDimorphism",
    "Plumage Dichromatism" = "Remove167HighPlumageDimorphism"
  )
  
  results_list <- list()
  
  # Process each bias correction
  for (bias_name in names(bias_types)) {
    pattern <- bias_types[[bias_name]]
    
    n_iterations = n_iterations_input # resets n_iterations to the input value if JackknifeSpecies was the previous bias_type
    
    if (is.null(pattern)) {
      # Handle full dataset case - look for non-downsampled results
      # For now, we'll skip this and handle it separately
      next
    }
    
    # Determine search directories based on bias type
    if (pattern == "JackknifeSpecies") {
      # JackknifeSpecies results are in PhylopathJackknife subdirectories
      search_dirs <- c(
        file.path(results_dir, "PhylopathJackknife", paste0(trait_set, " models")),
        # Also check without "models" suffix in case directory structure varies
        file.path(results_dir, "PhylopathJackknife", trait_set),
        # Check with underscores instead of spaces
        file.path(results_dir, "PhylopathJackknife", paste0(gsub(" ", "_", trait_set), " models")),
        file.path(results_dir, "PhylopathJackknife", gsub(" ", "_", trait_set))
      )
      n_iterations = NULL
    } else {
      # Regular downsampling results are in PhylopathDownsampled
      search_dirs <- c(
        file.path(results_dir, "PhylopathDownsampled", paste0(trait_set, " models")),
        # Also check without "models" suffix
        file.path(results_dir, "PhylopathDownsampled", trait_set),
        # Check with underscores instead of spaces
        file.path(results_dir, "PhylopathDownsampled", paste0(gsub(" ", "_", trait_set), " models")),
        file.path(results_dir, "PhylopathDownsampled", gsub(" ", "_", trait_set))
      )
    }
    
    # Find the detailed models file
    detailed_file <- NULL
    if (verbose) {
      cat("\nSearching for", bias_name, "(", pattern, ")\n")
    }
    
    for (search_dir in search_dirs) {
      if (verbose) {
        cat("  Checking directory:", search_dir, "\n")
      }
      
      if (dir.exists(search_dir)) {
        # Create pattern for finding files
        file_pattern <- paste0("detailed_models_", pattern, ".*\\.csv$")
        
        # Add n_iterations filter if specified
        if (!is.null(n_iterations)) {
          # Try multiple patterns for iteration specification
          file_patterns <- c(
            paste0("detailed_models_", pattern, "_", n_iterations, "_.*\\.csv$"),
            paste0("detailed_models_", pattern, "_n", n_iterations, "_.*\\.csv$")
          )
        } else {
          file_patterns <- file_pattern
        }
        
        # Try each pattern
        for (fp in file_patterns) {
          if (verbose) {
            cat("    Looking for pattern:", fp, "\n")
          }
          
          found_files <- list.files(search_dir, 
                                    pattern = fp,
                                    full.names = TRUE, recursive = FALSE)
          
          if (length(found_files) > 0) {
            detailed_file <- found_files
            if (verbose) {
              cat("    Found", length(found_files), "file(s)\n")
            }
            break
          }
        }
        
        if (!is.null(detailed_file) && length(detailed_file) > 0) {
          break
        }
      } else {
        if (verbose) {
          cat("    Directory does not exist\n")
        }
      }
    }
    
    if (!is.null(detailed_file) && length(detailed_file) > 0) {
      # Read the most recent file if multiple exist
      detailed_file <- detailed_file[length(detailed_file)]
      
      # Read and process the data
      detailed_models <- read.csv(detailed_file, stringsAsFactors = FALSE)
      
      # Store the results
      results_list[[bias_name]] <- list(
        detailed_models = detailed_models,
        filename = basename(detailed_file),
        trait_set = trait_set
      )
    }
  }
  
  # Also search for dimorphism results in subdirectories within the trait set directory
  dimorphism_dirs <- c("WingDimorphism", "PlumageDimorphism")
  
  for (dim_dir in dimorphism_dirs) {
    dim_path <- file.path(results_dir, "PhylopathDownsampled", paste0(trait_set, " models"), dim_dir)
    
    if (verbose) {
      cat("\nChecking for", dim_dir, "in:", dim_path, "\n")
    }
    
    if (dir.exists(dim_path)) {
      # Find all detailed_models files in subdirectory
      dim_files <- list.files(dim_path, 
                             pattern = "detailed_models.*\\.csv$",
                             full.names = TRUE,
                             recursive = FALSE)
      
      if (verbose && length(dim_files) > 0) {
        cat("  Found", length(dim_files), "file(s)\n")
      }
      
      # Filter files based on n_iterations if specified
      if (!is.null(n_iterations)) {
        dim_files <- dim_files[grepl(paste0("_", n_iterations, "_|_n", n_iterations, "_"), dim_files)]
      }
      
      for (dim_file in dim_files) {
        # Extract bias type from filename more cleanly
        filename <- basename(dim_file)
        # For files like detailed_models_Remove71HighWingDimorphism_n500_2025-06-12.csv
        if (grepl("Remove[0-9]+High", filename)) {
          # Extract just the dimorphism type with updated label
          if (dim_dir == "WingDimorphism") {
            bias_label <- "Wing Dimorphism"
          } else if (dim_dir == "PlumageDimorphism") {
            bias_label <- "Plumage Dichromatism"
          } else {
            bias_label <- gsub("(.+)Dimorphism", "\\1 Dimorphism", dim_dir)
          }
        } else {
          bias_label <- dim_dir
        }
        
        # Read and process the data
        detailed_models <- read.csv(dim_file, stringsAsFactors = FALSE)
        
        # Store the results
        results_list[[bias_label]] <- list(
          detailed_models = detailed_models,
          filename = basename(dim_file),
          trait_set = trait_set
        )
      }
    }
  }
  
  return(results_list)
}

#' Create a comprehensive bias robustness figure
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param output_file Path for saving the figure
create_bias_robustness_figure <- function(bias_results_list = NULL, 
                                          output_file = "Outputs/PhylopathFigures/bias_robustness_figure.pdf") {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Panel A: Summary of CB→FS path coefficients across bias corrections
  # Extract CB to FS coefficients from each bias correction
  coef_summary <- data.frame(
    Bias_Correction = character(),
    Mean_Coefficient = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    N_Iterations = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process each bias correction result
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    # Extract CB→FS coefficients from detailed models
    if (!is.null(result$detailed_models)) {
      # Find the relevant column - it should contain both Coop and FemaleSong
      # Try both orders: Coop to FemaleSong and FemaleSong to Coop
      cb_fs_columns <- c(
        grep("Coop.*to.*FemaleSong.*est", names(result$detailed_models), value = TRUE),
        grep("FemaleSong.*to.*Coop.*est", names(result$detailed_models), value = TRUE)
      )
      # Keep only CB to FS direction
      cb_fs_columns <- cb_fs_columns[grepl("Coop.*to.*FemaleSong", cb_fs_columns)]
      
      if (length(cb_fs_columns) > 0) {
        # Use the first matching column
        cb_fs_col <- cb_fs_columns[1]
        
        # Extract coefficients for each iteration
        cb_fs_coef <- result$detailed_models %>%
          filter(!is.na(.data[[cb_fs_col]])) %>%
          group_by(seed) %>%
          summarise(
            coef = mean(.data[[cb_fs_col]], na.rm = TRUE),
            .groups = "drop"
          )
        
        # Calculate summary statistics
        if (nrow(cb_fs_coef) > 0) {
          coef_summary <- rbind(coef_summary, data.frame(
            Bias_Correction = bias_name,
            Mean_Coefficient = mean(cb_fs_coef$coef),
            CI_Lower = quantile(cb_fs_coef$coef, 0.025),
            CI_Upper = quantile(cb_fs_coef$coef, 0.975),
            N_Iterations = length(unique(cb_fs_coef$seed)),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Create forest plot
  # Order the bias corrections for display
  coef_summary$Bias_Correction <- factor(coef_summary$Bias_Correction, 
                                         levels = rev(coef_summary$Bias_Correction))
  
  # Add reference value for full dataset (you can update this with actual value)
  full_dataset_estimate <- 0.56  # Update this with your actual full dataset estimate
  
  panel_a <- ggplot(coef_summary, aes(x = Mean_Coefficient, y = Bias_Correction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Add shaded region for full dataset estimate
    geom_rect(aes(xmin = full_dataset_estimate - 0.02, xmax = full_dataset_estimate + 0.02,
                  ymin = -Inf, ymax = Inf), 
              fill = "#2E86AB", alpha = 0.1, inherit.aes = FALSE) +
    geom_vline(xintercept = full_dataset_estimate, linetype = "dotted", 
               color = "#2E86AB", linewidth = 1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), 
                   height = 0.2, color = "black", linewidth = 0.8) +
    geom_point(aes(color = ifelse(CI_Lower > 0, "Significant", "Non-significant")), 
               size = 4) +
    scale_color_manual(values = c("Significant" = "#2E86AB", "Non-significant" = "#A7A7A7"),
                       guide = "none") +
    geom_text(aes(label = paste0(sprintf("%.3f", Mean_Coefficient), 
                                 ifelse(CI_Lower > 0, "*", ""),
                                 " (n=", N_Iterations, ")")), 
              vjust = -1, size = 3) +
    scale_x_continuous(limits = c(-0.1, 0.8), 
                       breaks = seq(0, 0.8, 0.2)) +
    labs(
      title = "A. Cooperative Breeding → Female Song Path Coefficient",
      subtitle = "Mean and 95% CI across downsampling iterations",
      x = "Path Coefficient",
      y = NULL
    ) +
    annotate("text", x = full_dataset_estimate + 0.03, y = 0.5, 
             label = "Full dataset\nestimate", 
             color = "#2E86AB", size = 3, hjust = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      panel.grid.major.y = element_blank()
    )
  
  # Extract model frequencies for bar plots
  # First, let's create a function to get top models from each bias correction
  get_top_models <- function(detailed_models, n_top = 15) {
    if (is.null(detailed_models) || nrow(detailed_models) == 0) return(NULL)
    
    top_models <- detailed_models %>%
      filter(delta_CICc < 2) %>%
      count(model) %>%
      arrange(desc(n)) %>%
      head(n_top) %>%
      mutate(
        percentage = n / length(unique(detailed_models$seed)) * 100,
        label = paste0(n, " (", round(percentage, 1), "%)")
      )
    
    return(top_models)
  }
  
  # Get top models for first two bias corrections for panel B and C
  bias_names <- names(bias_results_list)
  
  # Panel B: Top models for first bias correction (or geographic)
  if (length(bias_names) >= 1) {
    first_bias <- bias_names[1]
    top_models_1 <- get_top_models(bias_results_list[[first_bias]]$detailed_models)
    
    if (!is.null(top_models_1)) {
      panel_b <- ggplot(top_models_1, aes(x = n, y = reorder(model, n))) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = label), hjust = -0.1, size = 3) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(
          title = paste0("B. Top Models - ", first_bias),
          subtitle = paste0("Models appearing in Δ CICc < 2 set across ", 
                           bias_results_list[[first_bias]]$n_iterations, " iterations"),
          x = "Frequency in Top Model Set",
          y = NULL
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9),
          panel.grid.major.y = element_blank()
        )
    }
  }
  
  # Panel C: Top models for second bias correction (or territoriality)
  if (length(bias_names) >= 2) {
    second_bias <- bias_names[2]
    top_models_2 <- get_top_models(bias_results_list[[second_bias]]$detailed_models)
    
    if (!is.null(top_models_2)) {
      panel_c <- ggplot(top_models_2, aes(x = n, y = reorder(model, n))) +
        geom_bar(stat = "identity", fill = "#A7C957", alpha = 0.8) +
        geom_text(aes(label = label), hjust = -0.1, size = 3) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(
          title = paste0("C. Top Models - ", second_bias),
          subtitle = paste0("Models appearing in Δ CICc < 2 set across ", 
                           bias_results_list[[second_bias]]$n_iterations, " iterations"),
          x = "Frequency in Top Model Set",
          y = NULL
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9),
          panel.grid.major.y = element_blank()
        )
    }
  }
  
  # Combine panels with better layout
  # Handle cases where we might not have all panels
  if (exists("panel_b") && exists("panel_c")) {
    combined_plot <- (panel_a | panel_b) / panel_c +
      plot_layout(widths = c(1.2, 1), heights = c(1, 1)) +
      plot_annotation(
        title = "Robustness of Cooperative Breeding-Female Song Association to Sampling Biases",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
  } else if (exists("panel_b")) {
    combined_plot <- panel_a | panel_b +
      plot_layout(widths = c(1.2, 1)) +
      plot_annotation(
        title = "Robustness of Cooperative Breeding-Female Song Association to Sampling Biases",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
  } else {
    combined_plot <- panel_a +
      plot_annotation(
        title = "Robustness of Cooperative Breeding-Female Song Association to Sampling Biases",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
  }
  
  # Save figure with larger dimensions
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ggsave(output_file, combined_plot, width = 12, height = 9, dpi = 600)
  ggsave(gsub(".pdf", ".png", output_file), combined_plot, 
         width = 12, height = 9, dpi = 600)
  
  return(combined_plot)
}

# Alternative: Create a simpler version focusing just on the key coefficient
create_simple_robustness_figure <- function(coef_data, output_file) {
  
  # Create a clean coefficient comparison plot
  coef_plot <- ggplot(coef_data, aes(x = reorder(Correction, -Coefficient), 
                                     y = Coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0.56, linetype = "dotted", color = "#2E86AB",
               size = 1) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), 
                  width = 0.2, size = 0.8) +
    geom_point(size = 5, color = "#2E86AB") +
    geom_text(aes(label = sprintf("%.2f", Coefficient)), 
              vjust = -1.5, size = 4) +
    annotate("text", x = 4.5, y = 0.58, label = "Original estimate", 
             color = "#2E86AB", size = 3, hjust = 1) +
    scale_y_continuous(limits = c(-0.1, 0.8)) +
    labs(
      title = "Cooperative Breeding → Female Song Path Coefficient",
      subtitle = "Robustness across different bias corrections (500 iterations each)",
      x = "Bias Correction Applied",
      y = "Path Coefficient (95% CI)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      panel.grid.major.x = element_blank()
    )
  
  ggsave(output_file, coef_plot, width = 8, height = 6, dpi = 300)
  return(coef_plot)
}

#' Create a comprehensive model frequency comparison figure
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param output_file Path for saving the figure
#' @param n_top Number of top models to show per panel
create_model_frequency_figure <- function(bias_results_list = NULL,
                                        output_file = "Outputs/Figures/model_frequency_comparison.pdf",
                                        n_top = 10) {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Function to get top models
  get_top_models_df <- function(detailed_models, bias_name, n_top = 10) {
    if (is.null(detailed_models) || nrow(detailed_models) == 0) return(NULL)
    
    top_models <- detailed_models %>%
      filter(delta_CICc < 2) %>%
      count(model) %>%
      arrange(desc(n)) %>%
      head(n_top) %>%
      mutate(
        percentage = n / length(unique(detailed_models$seed)) * 100,
        label = paste0(n, " (", round(percentage, 1), "%)"),
        bias_correction = bias_name
      )
    
    return(top_models)
  }
  
  # Combine all top models
  all_top_models <- bind_rows(
    lapply(names(bias_results_list), function(bias_name) {
      get_top_models_df(bias_results_list[[bias_name]]$detailed_models, bias_name, n_top)
    })
  )
  
  # Create faceted plot
  if (!is.null(all_top_models) && nrow(all_top_models) > 0) {
    model_plot <- ggplot(all_top_models, aes(x = n, y = reorder(model, n))) +
      geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
      geom_text(aes(label = label), hjust = -0.1, size = 2.5) +
      facet_wrap(~ bias_correction, scales = "free", ncol = 2) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        title = "Top Models Across Different Bias Corrections",
        subtitle = "Models appearing in Δ CICc < 2 set across downsampling iterations",
        x = "Frequency in Top Model Set",
        y = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 8),
        panel.grid.major.y = element_blank()
      )
    
    # Save figure
    # Create output directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, model_plot, width = 14, height = 10, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), model_plot, 
           width = 14, height = 10, dpi = 300)
    
    return(model_plot)
  }
}

#' Create a violin plot showing coefficient distributions
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param output_file Path for saving the figure
create_coefficient_violin_plot <- function(bias_results_list = NULL,
                                         output_file = "Outputs/Figures/coefficient_distributions.pdf") {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Extract all coefficients into a long format data frame
  all_coefs <- data.frame()
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      # Find CB→FS coefficient column
      cb_fs_columns <- grep("Coop.*FemaleSong.*est", names(result$detailed_models), value = TRUE)
      
      if (length(cb_fs_columns) > 0) {
        cb_fs_col <- cb_fs_columns[1]
        
        # Extract coefficients for each model in each iteration
        coef_data <- result$detailed_models %>%
          filter(!is.na(.data[[cb_fs_col]])) %>%
          select(seed, model, coef = all_of(cb_fs_col)) %>%
          mutate(bias_correction = bias_name)
        
        all_coefs <- bind_rows(all_coefs, coef_data)
      }
    }
  }
  
  if (nrow(all_coefs) > 0) {
    # Create violin plot with overlaid box plot
    violin_plot <- ggplot(all_coefs, aes(x = bias_correction, y = coef)) +
      geom_violin(fill = "#2E86AB", alpha = 0.3, trim = FALSE) +
      geom_boxplot(width = 0.2, outlier.alpha = 0.3, alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0.56, linetype = "dotted", color = "#2E86AB", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
      scale_y_continuous(limits = c(-0.2, 1.0)) +
      labs(
        title = "Distribution of Cooperative Breeding → Female Song Path Coefficients",
        subtitle = "Across all models and iterations for each bias correction",
        x = "Bias Correction",
        y = "Path Coefficient"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11),
        panel.grid.major.x = element_blank()
      ) +
      annotate("text", x = 1, y = 0.58, label = "Full dataset estimate", 
               color = "#2E86AB", size = 3, hjust = 0)
    
    # Save figure
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, violin_plot, width = 10, height = 6, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), violin_plot, 
           width = 10, height = 6, dpi = 300)
    
    return(violin_plot)
  }
}

#' Create a ridge plot showing coefficient distributions by model type
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param output_file Path for saving the figure
create_coefficient_ridge_plot <- function(bias_results_list = NULL,
                                        output_file = "Outputs/Figures/coefficient_ridge_plot.pdf") {
  
  library(ggridges)
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Extract coefficients grouped by model structure
  model_coefs <- data.frame()
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      # Find CB→FS coefficient column
      cb_fs_columns <- grep("Coop.*FemaleSong.*est", names(result$detailed_models), value = TRUE)
      
      if (length(cb_fs_columns) > 0) {
        cb_fs_col <- cb_fs_columns[1]
        
        # Extract coefficients and simplify model names
        coef_data <- result$detailed_models %>%
          filter(!is.na(.data[[cb_fs_col]]), delta_CICc < 2) %>%
          select(seed, model, coef = all_of(cb_fs_col)) %>%
          mutate(
            bias_correction = bias_name,
            # Simplify model names
            has_CB_FS = grepl("COOP→FS", model),
            has_TERR_FS = grepl("TERR→FS", model),
            has_TERR_CB = grepl("TERR→COOP", model),
            has_MASS = grepl("MASS→", model),
            model_type = case_when(
              has_CB_FS & has_TERR_FS & has_TERR_CB ~ "CB→FS + T→FS + T→CB",
              has_CB_FS & has_TERR_FS & !has_TERR_CB ~ "CB→FS + T→FS",
              has_CB_FS & !has_TERR_FS & has_TERR_CB ~ "CB→FS + T→CB",
              has_CB_FS & !has_TERR_FS & !has_TERR_CB ~ "CB→FS only",
              !has_CB_FS ~ "No CB→FS",
              TRUE ~ "Other"
            )
          )
        
        model_coefs <- bind_rows(model_coefs, coef_data)
      }
    }
  }
  
  if (nrow(model_coefs) > 0) {
    # Create ridge plot
    ridge_plot <- ggplot(model_coefs, aes(x = coef, y = model_type, fill = model_type)) +
      geom_density_ridges(alpha = 0.8, scale = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0.56, linetype = "dotted", color = "#2E86AB", linewidth = 1) +
      scale_fill_manual(values = c(
        "CB→FS + T→FS + T→CB" = "#2E86AB",
        "CB→FS + T→FS" = "#5BA3C9",
        "CB→FS + T→CB" = "#88C0E6",
        "CB→FS only" = "#B5DDFF",
        "No CB→FS" = "#E63946",
        "Other" = "#A7A7A7"
      ), guide = "none") +
      facet_wrap(~ bias_correction, ncol = 2) +
      labs(
        title = "Path Coefficient Distributions by Model Structure",
        subtitle = "Showing CB→FS coefficients across different model types for each bias correction",
        x = "CB → FS Path Coefficient",
        y = "Model Structure"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)
      )
    
    # Save figure
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, ridge_plot, width = 12, height = 10, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), ridge_plot, 
           width = 12, height = 10, dpi = 300)
    
    return(ridge_plot)
  }
}

#' Create a heatmap of coefficients across models and bias corrections
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param output_file Path for saving the figure
create_coefficient_heatmap <- function(bias_results_list = NULL,
                                     output_file = "Outputs/Figures/coefficient_heatmap.pdf") {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Extract mean coefficients for top models
  heatmap_data <- data.frame()
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      # Get top 10 models
      top_models <- result$detailed_models %>%
        filter(delta_CICc < 2) %>%
        count(model) %>%
        arrange(desc(n)) %>%
        head(10) %>%
        pull(model)
      
      # Extract mean coefficients for these models
      cb_fs_col <- grep("Coop.*FemaleSong.*est", names(result$detailed_models), value = TRUE)[1]
      
      if (!is.na(cb_fs_col)) {
        model_means <- result$detailed_models %>%
          filter(model %in% top_models) %>%
          group_by(model) %>%
          summarise(
            mean_coef = mean(.data[[cb_fs_col]], na.rm = TRUE),
            n_iterations = n(),
            .groups = "drop"
          ) %>%
          mutate(
            bias_correction = bias_name,
            # Shorten model names
            model_short = gsub("_MASS→FS", "", model),
            model_short = gsub("_MASS→COOP", "", model_short),
            model_short = gsub("_MASS→TERR", "", model_short),
            model_short = gsub("→", "→", model_short)
          )
        
        heatmap_data <- bind_rows(heatmap_data, model_means)
      }
    }
  }
  
  if (nrow(heatmap_data) > 0) {
    # Create heatmap
    heatmap_plot <- ggplot(heatmap_data, 
                          aes(x = bias_correction, y = model_short, fill = mean_coef)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.3f", mean_coef)), size = 3) +
      scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2E86AB",
                          midpoint = 0, limits = c(-0.2, 0.8),
                          name = "Mean\nCoefficient") +
      labs(
        title = "Mean CB→FS Coefficients Across Top Models and Bias Corrections",
        subtitle = "Values shown for models with Δ CICc < 2",
        x = "Bias Correction",
        y = "Model Structure"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 11),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    # Save figure
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, heatmap_plot, width = 12, height = 8, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), heatmap_plot, 
           width = 12, height = 8, dpi = 300)
    
    return(heatmap_plot)
  }
}

#' Create a comprehensive summary figure with multiple panels
#'
#' @param bias_results_list List containing results from different bias corrections  
#' @param output_file Path for saving the figure
create_comprehensive_robustness_figure <- function(bias_results_list = NULL,
                                                 output_file = "Outputs/Figures/comprehensive_robustness.pdf") {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Create individual plots
  # 1. Forest plot
  coef_summary <- data.frame(
    Bias_Correction = character(),
    Mean_Coefficient = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    N_Iterations = integer(),
    stringsAsFactors = FALSE
  )
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      cb_fs_columns <- c(
        grep("Coop.*to.*FemaleSong.*est", names(result$detailed_models), value = TRUE),
        grep("FemaleSong.*to.*Coop.*est", names(result$detailed_models), value = TRUE)
      )
      cb_fs_columns <- cb_fs_columns[grepl("Coop.*to.*FemaleSong", cb_fs_columns)]
      
      if (length(cb_fs_columns) > 0) {
        cb_fs_col <- cb_fs_columns[1]
        cb_fs_coef <- result$detailed_models %>%
          filter(!is.na(.data[[cb_fs_col]])) %>%
          group_by(seed) %>%
          summarise(coef = mean(.data[[cb_fs_col]], na.rm = TRUE), .groups = "drop")
        
        if (nrow(cb_fs_coef) > 0) {
          coef_summary <- rbind(coef_summary, data.frame(
            Bias_Correction = bias_name,
            Mean_Coefficient = mean(cb_fs_coef$coef),
            CI_Lower = quantile(cb_fs_coef$coef, 0.025),
            CI_Upper = quantile(cb_fs_coef$coef, 0.975),
            N_Iterations = length(unique(cb_fs_coef$seed)),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Simplified forest plot
  panel_a <- ggplot(coef_summary, aes(x = Mean_Coefficient, y = reorder(Bias_Correction, Mean_Coefficient))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0.56, linetype = "dotted", color = "#2E86AB", linewidth = 1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 0.8) +
    geom_point(size = 4, color = "#2E86AB") +
    geom_text(aes(label = sprintf("%.3f", Mean_Coefficient)), vjust = -1, size = 3) +
    scale_x_continuous(limits = c(0, 0.8)) +
    labs(title = "A. CB→FS Path Coefficients", x = "Coefficient (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # 2. Violin plot for coefficient distributions
  all_coefs <- data.frame()
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    if (!is.null(result$detailed_models)) {
      cb_fs_columns <- grep("Coop.*to.*FemaleSong.*est", names(result$detailed_models), value = TRUE)
      if (length(cb_fs_columns) > 0) {
        coef_data <- result$detailed_models %>%
          filter(!is.na(.data[[cb_fs_columns[1]]])) %>%
          select(coef = all_of(cb_fs_columns[1])) %>%
          mutate(bias_correction = bias_name)
        all_coefs <- bind_rows(all_coefs, coef_data)
      }
    }
  }
  
  panel_b <- ggplot(all_coefs, aes(x = bias_correction, y = coef)) +
    geom_violin(fill = "#2E86AB", alpha = 0.3) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    geom_hline(yintercept = 0.56, linetype = "dotted", color = "#2E86AB") +
    labs(title = "B. Coefficient Distributions", x = NULL, y = "Path Coefficient") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9)
    )
  
  # 3. Model frequency summary
  model_freq_summary <- data.frame()
  for (bias_name in names(bias_results_list)[1:min(3, length(bias_results_list))]) {
    result <- bias_results_list[[bias_name]]
    if (!is.null(result$detailed_models)) {
      top_model <- result$detailed_models %>%
        filter(delta_CICc < 2) %>%
        count(model) %>%
        arrange(desc(n)) %>%
        head(1) %>%
        mutate(
          bias_correction = bias_name,
          model_simple = gsub(".*_", "", model)
        )
      model_freq_summary <- bind_rows(model_freq_summary, top_model)
    }
  }
  
  panel_c <- ggplot(model_freq_summary, aes(x = bias_correction, y = n, fill = bias_correction)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = model_simple), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c("#2E86AB", "#5BA3C9", "#88C0E6"), guide = "none") +
    labs(title = "C. Top Model Frequency", x = NULL, y = "Frequency") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9)
    )
  
  # Combine panels
  combined_plot <- panel_a / (panel_b | panel_c) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = "Robustness of Cooperative Breeding-Female Song Association",
      subtitle = "Across multiple bias corrections and model structures",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  # Save figure
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ggsave(output_file, combined_plot, width = 10, height = 10, dpi = 300)
  ggsave(gsub(".pdf", ".png", output_file), combined_plot, 
         width = 10, height = 10, dpi = 300)
  
  return(combined_plot)
}

#' Create a DAG plot from phylopath results
#'
#' @param phylopath_result Result from phylopath analysis or path to RDS file
#' @param output_file Path for saving the figure
#' @param title Title for the plot
create_phylopath_dag <- function(phylopath_result = NULL,
                                output_file = "Outputs/Figures/phylopath_dag.pdf",
                                title = "Phylogenetic Path Analysis") {
  
  library(phylopath)
  
  # If result is a path, load it
  if (is.character(phylopath_result)) {
    if (file.exists(phylopath_result)) {
      phylopath_result <- readRDS(phylopath_result)
    } else {
      stop("File not found: ", phylopath_result)
    }
  }
  
  # Extract the phylopath result object if it's wrapped
  if (!is.null(phylopath_result$result)) {
    result <- phylopath_result$result
  } else {
    result <- phylopath_result
  }
  
  # Set up layout positions for 4-variable DAG
  phylopath_map_positions <- data.frame(
    name = c("FemaleSong_Agg01", "HighConfidence_Coop", 
             "TerritorialityWeakVsStrong", "logMass_AVONET"),
    x = c(8, 2, 5, 5),
    y = c(9, 9, 5, 1),
    stringsAsFactors = FALSE
  )
  
  # Create conditional averaged model
  cond_avg <- average(result, cut_off = 2, avg_method = "conditional")
  
  # Create the DAG plot
  dag_plot <- plot(cond_avg,
                  manual_layout = phylopath_map_positions,
                  text_size = 4,
                  box_x = 24, box_y = 18,
                  edge_width = 2) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(30, 30, 30, 30)) +
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, dag_plot, width = 8, height = 8, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), dag_plot, 
           width = 8, height = 8, dpi = 300)
  }
  
  return(dag_plot)
}

#### Moved: Create an enhanced DAG plot with better styling ----
#'
#' Function create_enhanced_dag() is now in create_enhanced_DAG.R



#' Extract species counts from detailed models data
#'
#' @param detailed_models Data frame from phylopath downsampling
#' @return Numeric species count
get_species_count <- function(detailed_models) {
  if (!is.null(detailed_models$nSpecies)) {
    return(unique(detailed_models$nSpecies)[1])
  }
  # If nSpecies column doesn't exist, estimate from data
  # This is approximate - better to have actual count
  return(NA)
}

#' Create an improved coefficient distribution plot
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param paths_to_show Which paths to display (default shows main paths)
#' @param output_file Path for saving the figure  
create_all_paths_coefficient_plot <- function(bias_results_list = NULL,
                                            paths_to_show = c("CB_to_FS", "TERR_to_FS", 
                                                            "TERR_to_CB", "FS_to_CB", "MASS_to_FS", "MASS_to_CB", "MASS_to_TERR"),
                                            output_file = "Outputs/Figures/all_paths_coefficients.pdf") {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Extract coefficients for all paths
  all_path_coefs <- data.frame()
  
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    
    if (!is.null(result$detailed_models)) {
      df <- result$detailed_models
      
      # Find all coefficient columns (ending in _est)
      coef_cols <- grep("_est$", names(df), value = TRUE)
      
      # Extract coefficients for each path
      for (col in coef_cols) {
        # Parse the path name
        path_parts <- strsplit(col, "_to_")[[1]]
        if (length(path_parts) >= 2) {
          from_var <- gsub("_est$", "", path_parts[length(path_parts)])
          to_var <- paste(path_parts[1:(length(path_parts)-1)], collapse = "_to_")
          
          # Simplify variable names
          from_simple <- gsub("HighConfidence_Coop", "CB", from_var)
          from_simple <- gsub("FemaleSong_Agg01", "FS", from_simple)
          from_simple <- gsub("TerritorialityWeakVsStrong", "TERR", from_simple)
          from_simple <- gsub("logMass_AVONET", "Mass", from_simple)
          
          to_simple <- gsub("HighConfidence_Coop", "CB", to_var)
          to_simple <- gsub("FemaleSong_Agg01", "FS", to_simple) 
          to_simple <- gsub("TerritorialityWeakVsStrong", "TERR", to_simple)
          to_simple <- gsub("logMass_AVONET", "Mass", to_simple)
          
          path_name <- paste0(from_simple, " → ", to_simple)
          
          # Extract non-NA values
          coef_vals <- df[[col]][!is.na(df[[col]])]
          
          if (length(coef_vals) > 0) {
            path_df <- data.frame(
              bias_correction = bias_name,
              path = path_name,
              coefficient = coef_vals,
              stringsAsFactors = FALSE
            )
            all_path_coefs <- bind_rows(all_path_coefs, path_df)
          }
        }
      }
    }
  }
  
  # Create the plot
  if (nrow(all_path_coefs) > 0) {
    # Order paths by mean coefficient
    path_order <- all_path_coefs %>%
      group_by(path) %>%
      summarise(mean_coef = mean(coefficient), .groups = "drop") %>%
      arrange(mean_coef) %>%
      pull(path)
    
    all_path_coefs$path <- factor(all_path_coefs$path, levels = path_order)
    
    # Create violin plot
    coef_plot <- ggplot(all_path_coefs, aes(x = coefficient, y = path)) +
      geom_violin(aes(fill = bias_correction), alpha = 0.6, scale = "width") +
      geom_boxplot(aes(group = interaction(path, bias_correction)), 
                   width = 0.1, alpha = 0.8, outlier.size = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~ bias_correction, ncol = 1) +
      labs(
        title = "Distribution of Path Coefficients Across Bias Corrections",
        subtitle = "Based on conditional averaging of models with Δ CICc < 2",
        x = "Path Coefficient",
        y = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9),
        legend.position = "none"
      )
    
    # Save figure
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, coef_plot, width = 10, height = 12, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), coef_plot, 
           width = 10, height = 12, dpi = 300)
    
    return(coef_plot)
  }
}

#' Create model selection bar plot showing CICc values
#'
#' @param detailed_models Data frame from phylopath analysis
#' @param n_models Number of top models to show (default 20)
#' @param output_file Path for saving the figure
create_model_selection_plot <- function(detailed_models,
                                      n_models = 20,
                                      output_file = "Outputs/Figures/model_selection.pdf",
                                      title = "Model Selection Results") {
  
  # Get unique models and their best CICc
  model_summary <- detailed_models %>%
    group_by(model) %>%
    summarise(
      best_CICc = min(CICc),
      best_delta_CICc = min(delta_CICc),
      frequency = n(),
      .groups = "drop"
    ) %>%
    arrange(best_CICc) %>%
    head(n_models)
  
  # Create the plot
  model_plot <- ggplot(model_summary, 
                      aes(x = reorder(model, -best_CICc), 
                          y = best_CICc,
                          fill = best_delta_CICc < 2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", best_delta_CICc)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#A7A7A7"),
                      labels = c("TRUE" = "Δ CICc < 2", "FALSE" = "Δ CICc ≥ 2"),
                      name = "Model Set") +
    labs(
      title = title,
      subtitle = paste("Top", n_models, "models by CICc. Values show Δ CICc"),
      x = NULL,
      y = "CICc"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, model_plot, width = 10, height = 8, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), model_plot, 
           width = 10, height = 8, dpi = 300)
  }
  
  return(model_plot)
}


#### Moved: Create an improved forest plot with species counts ----
#'
#' Function create_forest_plot_with_counts() is now in create_forest_plot_with_counts_updated.R



#' Create final publication-ready robustness figure
#'
#' This creates a figure matching the user's caption with:
#' Panel A: DAG with conditional-averaged coefficients  
#' Panel B: Model selection bar plot
#' Panels C-E: Violin plots for different bias corrections
#'
#' @param output_file Path for saving the figure
create_publication_robustness_figure <- function(output_file = "Outputs/Figures/publication_robustness_figure.pdf") {
  
  # Extract all results
  bias_results <- extract_phylopath_results()
  
  # Panel A: DAG (placeholder for now - needs actual phylopath result)
  # For now, create a text placeholder
  panel_a <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = "Panel A: DAG\n(Requires full dataset\nphylopath result)", 
             size = 6) +
    theme_void() +
    ggtitle("A. Directed Acyclic Graph") +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  # Panel B: Model selection for one bias correction (e.g., Holarctic)
  if ("Geographic (Holarctic)" %in% names(bias_results)) {
    panel_b <- create_model_selection_plot(
      bias_results[["Geographic (Holarctic)"]]$detailed_models,
      n_models = 20,
      output_file = NULL,
      title = "B. Model Selection Results"
    ) + 
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  }
  
  # Panels C-E: Coefficient distributions for specific bias corrections
  # Extract CB→FS coefficients for three key bias corrections
  create_bias_violin <- function(bias_name, panel_letter) {
    if (bias_name %in% names(bias_results)) {
      result <- bias_results[[bias_name]]
      df <- result$detailed_models
      
      # Find CB→FS coefficient column
      cb_fs_col <- grep("Coop.*to.*FemaleSong.*est", names(df), value = TRUE)[1]
      
      if (!is.na(cb_fs_col)) {
        coef_data <- df %>%
          filter(!is.na(.data[[cb_fs_col]])) %>%
          select(coefficient = all_of(cb_fs_col))
        
        n_removed <- case_when(
          grepl("Holarctic", bias_name) ~ 83,
          grepl("Tropical", bias_name) ~ 24,
          grepl("Plumage", bias_name) ~ 167,
          TRUE ~ NA_real_
        )
        
        subtitle_text <- if (!is.na(n_removed)) {
          paste0(n_removed, " species removed per iteration")
        } else {
          "Downsampled dataset"
        }
        
        violin_plot <- ggplot(coef_data, aes(x = "", y = coefficient)) +
          geom_violin(fill = "#2E86AB", alpha = 0.3) +
          geom_boxplot(width = 0.2, outlier.alpha = 0.5) +
          geom_hline(yintercept = 0.56, linetype = "solid", color = "red", linewidth = 1) +
          scale_y_continuous(limits = c(0, 1)) +
          labs(
            title = paste0(panel_letter, ". ", bias_name),
            subtitle = subtitle_text,
            x = NULL,
            y = "CB → FS Coefficient"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 9),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 10)
          )
        
        return(violin_plot)
      }
    }
    return(NULL)
  }
  
  panel_c <- create_bias_violin("Geographic (Holarctic)", "C")
  panel_d <- create_bias_violin("Geographic (Tropical)", "D")
  panel_e <- create_bias_violin("Plumage Dimorphism", "E")
  
  # Combine all panels
  # Top row: DAG and model selection
  # Bottom row: Three violin plots
  top_row <- panel_a | panel_b
  bottom_row <- panel_c | panel_d | panel_e
  
  combined_plot <- top_row / bottom_row +
    plot_layout(heights = c(1, 0.8)) +
    plot_annotation(
      title = "Robustness of the cooperative breeding-female song association to data availability biases",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  # Save figure
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, combined_plot, width = 14, height = 10, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), combined_plot, 
           width = 14, height = 10, dpi = 300)
  }
  
  return(combined_plot)
}

#' Create a plot showing all path coefficients like existing violin-box plots
#'
#' This replicates the style of existing plots in PhylopathDownsampled
#' Shows horizontal violin/box plots for all paths with vertical reference lines
#'
#' @param detailed_models Data frame from a single bias correction
#' @param reference_values Named vector of reference values for each path
#' @param output_file Path for saving the figure
#' @param title Title for the plot
create_all_paths_horizontal_plot <- function(detailed_models,
                                           reference_values = NULL,
                                           output_file = NULL,
                                           title = "Distribution of Path Coefficients Across Iterations") {
  
  # Extract all coefficient columns
  coef_cols <- grep("_est$", names(detailed_models), value = TRUE)
  
  # Create long format data
  path_data <- data.frame()
  
  for (col in coef_cols) {
    # Parse path name
    path_parts <- strsplit(col, "_to_")[[1]]
    if (length(path_parts) >= 2) {
      from_var <- gsub("_est$", "", path_parts[length(path_parts)])
      to_var <- paste(path_parts[1:(length(path_parts)-1)], collapse = "_to_")
      
      # Create readable labels
      from_label <- case_when(
        grepl("HighConfidence_Coop", from_var) ~ "Cooperation",
        grepl("FemaleSong", from_var) ~ "Female Song",
        grepl("TerritorialityWeakVsStrong", from_var) ~ "Strong Territoriality",
        grepl("logMass", from_var) ~ "Body Mass",
        TRUE ~ from_var
      )
      
      to_label <- case_when(
        grepl("HighConfidence_Coop", to_var) ~ "Cooperation",
        grepl("FemaleSong", to_var) ~ "Female Song",
        grepl("TerritorialityWeakVsStrong", to_var) ~ "Strong Territoriality",
        grepl("logMass", to_var) ~ "Body Mass",
        TRUE ~ to_var
      )
      
      path_label <- paste(from_label, "→", to_label)
      
      # Extract values
      values <- detailed_models[[col]][!is.na(detailed_models[[col]])]
      
      if (length(values) > 0) {
        path_df <- data.frame(
          path = path_label,
          coefficient = values,
          stringsAsFactors = FALSE
        )
        path_data <- bind_rows(path_data, path_df)
      }
    }
  }
  
  # Calculate mean coefficients for ordering
  path_order <- path_data %>%
    group_by(path) %>%
    summarise(mean_coef = mean(coefficient), .groups = "drop") %>%
    arrange(mean_coef) %>%
    pull(path)
  
  path_data$path <- factor(path_data$path, levels = path_order)
  
  # Create the plot
  p <- ggplot(path_data, aes(x = coefficient, y = path)) +
    geom_violin(fill = "lightblue", alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      subtitle = paste("Based on", length(unique(detailed_models$seed)), 
                      "downsampling iterations (conditional averaging)"),
      x = "Mean Coefficient",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.grid.major.y = element_blank()
    )
  
  # Add reference lines if provided
  if (!is.null(reference_values)) {
    for (i in seq_along(reference_values)) {
      ref_path <- names(reference_values)[i]
      ref_val <- reference_values[i]
      
      # Find matching path in data
      if (any(grepl(ref_path, levels(path_data$path)))) {
        p <- p + geom_vline(xintercept = ref_val, 
                           linetype = "solid", 
                           color = "red", 
                           linewidth = 1)
      }
    }
    
    # Add legend for reference lines
    p <- p + 
      annotate("text", x = max(path_data$coefficient) * 0.9, 
               y = length(levels(path_data$path)) - 0.5,
               label = "Red lines show full dataset conditional averages",
               color = "red", size = 3, hjust = 1)
  }
  
  # Save if output file specified
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), p, 
           width = 10, height = 8, dpi = 300)
  }
  
  return(p)
}

#### Moved: Create a model consistency heatmap across bias corrections ----
#'
# Moved to create_bias_model_consistency_heatmap_ModelsOnY.R (earlier version closer to what was in here) and create_bias_model_consistency_heatmap.R (more polished version)



#' Add usage examples to the script
#'
#' @examples
#' # Load the script
#' source('create_phylopath_bias_robustness_figure.R')
#' 
#' # Extract results
#' results <- extract_phylopath_results()
#' 
#' # Create publication figure
#' pub_fig <- create_publication_robustness_figure()
#' 
#' # Create individual plots
#' forest <- create_forest_plot_with_counts()
#' 
#' # For a specific bias correction
#' if ("Geographic (Holarctic)" %in% names(results)) {
#'   # Model selection plot
#'   models <- create_model_selection_plot(
#'     results[["Geographic (Holarctic)"]]$detailed_models,
#'     title = "Model Selection - Geographic (Holarctic)"
#'   )
#'   
#'   # All paths plot
#'   all_paths <- create_all_paths_horizontal_plot(
#'     results[["Geographic (Holarctic)"]]$detailed_models,
#'     title = "Path Coefficients - Geographic (Holarctic)"
#'   )
#' }

# Example usage - minimal viable version
# To test the figure with your current data:
# 
# # Option 1: Create the 3-panel figure with forest plot and model frequencies
# fig <- create_bias_robustness_figure()
# 
# # Option 2: Create a comprehensive model frequency comparison across all bias corrections
# model_fig <- create_model_frequency_figure()
# 
# # Option 3: Test extraction and see what data is available
# results <- extract_phylopath_results()
# names(results)  # See which bias corrections were found
# 
# # Option 4: Create just the simple forest plot
# if (length(results) > 0) {
#   # Extract actual coefficients for the simple plot
#   coef_data <- data.frame()
#   for (bias_name in names(results)) {
#     cb_fs_col <- grep("Coop.*FemaleSong.*est", 
#                       names(results[[bias_name]]$detailed_models), 
#                       value = TRUE)[1]
#     if (!is.na(cb_fs_col)) {
#       coefs <- results[[bias_name]]$detailed_models[[cb_fs_col]]
#       coef_data <- rbind(coef_data, data.frame(
#         Correction = bias_name,
#         Coefficient = mean(coefs, na.rm = TRUE),
#         CI_Lower = quantile(coefs, 0.025, na.rm = TRUE),
#         CI_Upper = quantile(coefs, 0.975, na.rm = TRUE)
#       ))
#     }
#   }
#   
#   simple_fig <- create_simple_robustness_figure(
#     coef_data, 
#     "Outputs/Figures/test_robustness_figure.png"
#   )
# }
