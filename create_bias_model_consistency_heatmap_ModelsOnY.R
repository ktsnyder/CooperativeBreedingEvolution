# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Note: This function requires the extract_phylopath_results() function
# from create_phylopath_bias_robustness_figure.R
# Source that file first or ensure the function is available

#' Create a model consistency heatmap across bias corrections
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param top_n Number of top models to consider from any bias correction
#' @param output_file Path for saving the figure
#' @param title Title for the heatmap
#' @param full_dataset_best_models Character vector of best models from full dataset analysis
create_bias_model_consistency_heatmap <- function(bias_results_list = NULL,
                                             top_n = 10,
                                             output_file = "Outputs/Figures/model_consistency_heatmap.pdf",
                                             title = "Model Consistency Across Bias Corrections",
                                             full_dataset_best_models = NULL) {
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results()
  }
  
  # Load full dataset best models and delta CICc values if not provided
  full_dataset_delta_cicc <- NULL
  if (is.null(full_dataset_best_models)) {
    if (file.exists("phylopath_full_dataset_result.rds")) {
      full_result <- readRDS("phylopath_full_dataset_result.rds")
      full_dataset_best_models <- full_result$CICsub2_models
      # Extract delta CICc values for ordering
      if (!is.null(full_result$result)) {
        model_summary <- summary(full_result$result)
        full_dataset_delta_cicc <- model_summary$delta_CICc
        names(full_dataset_delta_cicc) <- rownames(model_summary)
      }
    }
  }
  
  # Initialize data structures
  all_models <- character()
  model_data <- list()
  
  # First pass: collect all unique models that appear in top 10 for any bias correction
  for (bias_name in names(bias_results_list)) {
    result <- bias_results_list[[bias_name]]
    detailed_models <- result$detailed_models
    
    if (!is.null(detailed_models) && nrow(detailed_models) > 0) {
      # Group by model and calculate frequency
      model_freq <- detailed_models %>%
        filter(delta_CICc < 2) %>%  # Only consider models within 2 CICc of best
        group_by(model) %>%
        summarise(
          n_significant = n(),
          n_best = sum(delta_CICc == 0, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(n_significant))
      
      # Get top N models for this bias correction
      top_models <- head(model_freq$model, top_n)
      all_models <- unique(c(all_models, top_models))
      
      # Store the frequency data
      model_data[[bias_name]] <- model_freq
    }
  }
  
  # Also include full dataset best models if they're not already in the list
  if (!is.null(full_dataset_best_models)) {
    all_models <- unique(c(all_models, full_dataset_best_models))
  }
  
  # Create matrix for heatmap
  heatmap_matrix <- matrix(0, 
                           nrow = length(all_models), 
                           ncol = length(model_data),
                           dimnames = list(all_models, names(model_data)))
  
  best_matrix <- matrix(0, 
                        nrow = length(all_models), 
                        ncol = length(model_data),
                        dimnames = list(all_models, names(model_data)))
  
  # Fill in the matrices
  for (bias_name in names(model_data)) {
    freq_data <- model_data[[bias_name]]
    total_iterations <- length(unique(bias_results_list[[bias_name]]$detailed_models$seed))
    
    for (model in all_models) {
      if (model %in% freq_data$model) {
        row_data <- freq_data[freq_data$model == model, ]
        # Number of iterations where model was significant (not percentage yet)
        heatmap_matrix[model, bias_name] <- row_data$n_significant
        # Number of times model was best
        best_matrix[model, bias_name] <- row_data$n_best
      }
    }
  }
  
  # Order models by delta CICc from full dataset analysis (descending - worst first)
  if (!is.null(full_dataset_delta_cicc)) {
    # Get delta CICc values for our models, using large value for models not in full dataset
    model_delta_values <- sapply(rownames(heatmap_matrix), function(m) {
      if (m %in% names(full_dataset_delta_cicc)) {
        full_dataset_delta_cicc[m]
      } else {
        999  # Large value for models not in full dataset analysis
      }
    })
    # Order by delta CICc (descending), then by total significance for ties
    model_order <- order(-model_delta_values, -rowSums(heatmap_matrix))
  } else {
    # Fallback: order by whether they're in full dataset best models, then by total significance
    model_in_full <- rownames(heatmap_matrix) %in% full_dataset_best_models
    model_order <- order(model_in_full, -rowSums(heatmap_matrix))
  }
  
  heatmap_matrix <- heatmap_matrix[model_order, , drop = FALSE]
  best_matrix <- best_matrix[model_order, , drop = FALSE]
  
  # Prepare data for ggplot
  heatmap_long <- as.data.frame(heatmap_matrix) %>%
    mutate(Model = rownames(.)) %>%
    pivot_longer(cols = -Model, names_to = "Bias_Correction", values_to = "N_Significant")
  
  best_long <- as.data.frame(best_matrix) %>%
    mutate(Model = rownames(.)) %>%
    pivot_longer(cols = -Model, names_to = "Bias_Correction", values_to = "N_Best")
  
  # Merge the data
  plot_data <- heatmap_long %>%
    left_join(best_long, by = c("Model", "Bias_Correction")) %>%
    mutate(
      Model = factor(Model, levels = rownames(heatmap_matrix)),
      Bias_Correction = factor(Bias_Correction, levels = colnames(heatmap_matrix)),
      Model_Label = ifelse(Model %in% full_dataset_best_models, 
                           paste0(Model, " ★"), 
                           as.character(Model)),
      # Convert 0 values to NA so they appear white
      N_Significant = ifelse(N_Significant == 0, NA, N_Significant)
    )
  
  # Create the heatmap with improved color gradient
  p <- ggplot(plot_data, aes(x = Bias_Correction, y = Model)) +
    geom_tile(aes(fill = N_Significant), color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(N_Best > 0, N_Best, ""),
                  color = N_Significant > 250), 
              size = 3, fontface = "bold") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_fill_gradient(
      low = "#B8E0F0",  # Slightly darker light blue
      high = "#001233", # Darkest blue possible
      limits = c(1, 500),
      breaks = c(1, 100, 200, 300, 400, 500),
      name = "# Iterations\n(ΔCICc < 2)",
      na.value = "white"  # Cells with 0 iterations will be white
    ) +
    scale_y_discrete(labels = function(x) {
      ifelse(x %in% full_dataset_best_models, paste0(x, " ★"), x)
    }) +
    labs(
      title = title,
      subtitle = paste("Models marked with ★ had ΔCICc < 2 in full dataset analysis"),
      x = "Bias Correction",
      y = "Model",
      caption = "Numbers show count of iterations where model was best. White cells = model never within ΔCICc < 2."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10))
    )
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Calculate appropriate height based on number of models
    plot_height <- max(8, 5 + length(all_models) * 0.3)
    
    ggsave(output_file, p, width = 10, height = plot_height, dpi = 300)
    ggsave(gsub(".pdf", ".png", output_file), p, 
           width = 10, height = plot_height, dpi = 300)
  }
  
  return(p)
}
