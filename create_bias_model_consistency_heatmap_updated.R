# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

# Note: This function requires the extract_phylopath_results() function
# from create_phylopath_bias_robustness_figure.R
# Source that file first or ensure the function is available
# Also source phylopath_helper_functions.R for helper functions
# 
# 7/8/2025 - in extract_phylopath_results(), set results_dir = "Outputs/" which will work with new recursive searching in extract_phylopath_results()
# 7/23/2025 - reordered downsampled groups so Jackknife by species is in the last row of the matrix if it is present

# Source helper functions if available
if (file.exists("phylopath_helper_functions.R")) {
  source("phylopath_helper_functions.R")
}

#' Create a model consistency heatmap across bias corrections
#'
#' @param bias_results_list List containing results from different bias corrections
#' @param top_n Number of top models to consider from any bias correction
#' @param output_file Path for saving the figure
#' @param title Title for the heatmap
#' @param full_dataset_best_models Character vector of best models from full dataset analysis
#' @param trait_set Character string specifying which trait combination to use
create_bias_model_consistency_heatmap <- function(bias_results_list = NULL,
                                             top_n = 10,
                                             output_file = NULL,
                                             title = "Model Consistency Across Bias Corrections",
                                             full_dataset_best_models = NULL,
                                             trait_set = NULL) {
  
  # Set default trait set for backward compatibility
  if (is.null(trait_set)) {
    trait_set <- "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
  }
  
  # Generate output file path if not provided
  if (is.null(output_file)) {
    # Get trait set label for filename
    trait_label <- if (exists("get_trait_set_label")) {
      get_trait_set_label(trait_set)
    } else {
      gsub(" ", "_", trait_set)
    }
    output_file <- file.path("Outputs", "PhylopathFigures", trait_label, "model_consistency_heatmap_byrate.pdf")
  }
  
  # If no results provided, extract from default location
  if (is.null(bias_results_list)) {
    bias_results_list <- extract_phylopath_results(results_dir = "Outputs/", trait_set = trait_set)
  }
  
  # Load full dataset best models and delta CICc values if not provided
  full_dataset_delta_cicc <- NULL
  if (is.null(full_dataset_best_models)) {
    # Use helper function to find the appropriate full dataset result
    if (exists("find_or_create_full_dataset_result")) {
      full_result <- find_or_create_full_dataset_result(trait_set, run_if_missing = FALSE)
    } else {
      # Fallback to old method # commented out because this is just asking for errors
      # if (file.exists("phylopath_full_dataset_result.rds")) {
      #   full_result <- readRDS("phylopath_full_dataset_result.rds")
      # } else {
        full_result <- NULL
      # }
    }
    
    if (!is.null(full_result)) {
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
  
  target <- "Jackknife by species"
  if (target %in% names(model_data)) {
    other_names <- names(model_data)[names(model_data) != target]
    model_data <- model_data[c(target, other_names)]
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
  
  # Create a vector to store total iterations for each bias correction
  total_iterations_vec <- numeric()
  
  # Fill in the matrices
  for (bias_name in names(model_data)) {
    freq_data <- model_data[[bias_name]]
    total_iterations <- length(unique(bias_results_list[[bias_name]]$detailed_models$seed))
    total_iterations_vec[bias_name] <- total_iterations
    
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
    model_order <- order(model_delta_values, rowSums(heatmap_matrix))
  } else {
    # Fallback: order by whether they're in full dataset best models, then by total significance
    model_in_full <- rownames(heatmap_matrix) %in% full_dataset_best_models
    model_order <- order(model_in_full, -rowSums(heatmap_matrix))
  }
  
  heatmap_matrix <- heatmap_matrix[model_order, , drop = FALSE]
  best_matrix <- best_matrix[model_order, , drop = FALSE]
  
  # Convert to proportions
  heatmap_matrix_prop <- heatmap_matrix
  best_matrix_prop <- best_matrix
  
  for (bias_name in names(model_data)) {
    if (bias_name %in% colnames(heatmap_matrix)) {
      heatmap_matrix_prop[, bias_name] <- heatmap_matrix[, bias_name] / total_iterations_vec[bias_name]
      best_matrix_prop[, bias_name] <- best_matrix[, bias_name] / total_iterations_vec[bias_name]
    }
  }
  
  # Prepare data for ggplot
  heatmap_long <- as.data.frame(heatmap_matrix_prop) %>%
    mutate(Model = rownames(.)) %>%
    pivot_longer(cols = -Model, names_to = "Bias_Correction", values_to = "Prop_Significant")
  
  best_long <- as.data.frame(best_matrix_prop) %>%
    mutate(Model = rownames(.)) %>%
    pivot_longer(cols = -Model, names_to = "Bias_Correction", values_to = "Prop_Best")
  
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
      Prop_Significant = ifelse(Prop_Significant == 0, NA, Prop_Significant)
    )
  
  # Create presence/absence matrix for rates
  rates <- c("COOP\u2192FS", "FS\u2192COOP", "TERR\u2192FS", "TERR\u2192COOP", "MASS\u2192FS", "MASS\u2192COOP", "MASS\u2192TERR")
  
  # Get the ordered models from the heatmap matrix
  ordered_models <- rownames(heatmap_matrix)
  
  # Create a matrix to store presence/absence
  rate_matrix <- matrix(FALSE, nrow = length(rates), ncol = length(ordered_models),
                        dimnames = list(rates, ordered_models))
  
  # Fill in the matrix based on model names
  for (i in seq_along(ordered_models)) {
    model <- ordered_models[i]
    rate_matrix["COOP\u2192FS", i] <- grepl("COOP→FS", model)
    rate_matrix["FS\u2192COOP", i] <- grepl("FS→COOP", model)
    rate_matrix["TERR\u2192FS", i] <- grepl("TERR→FS", model)
    rate_matrix["TERR\u2192COOP", i] <- grepl("TERR→COOP", model)
    rate_matrix["MASS\u2192FS", i] <- grepl("MASS→FS", model)
    rate_matrix["MASS\u2192COOP", i] <- grepl("MASS→COOP", model)
    rate_matrix["MASS\u2192TERR", i] <- grepl("MASS→TERR", model)
  }
  
  # Convert to data frame for ggplot
  rate_data <- as.data.frame(rate_matrix) %>%
    mutate(Rate = rownames(.)) %>%
    pivot_longer(cols = -Rate, names_to = "Model", values_to = "Present") %>%
    mutate(
      Model = factor(Model, levels = ordered_models),
      Rate = factor(Rate, levels = rev(rates))  # Reverse to show in correct order
    )
  
  # Create a factor for missing values to use in a dummy aesthetic
  plot_data$is_missing <- factor(ifelse(is.na(plot_data$Prop_Significant), "0", "Present"))
  
  # Create dummy data for all legend items
  legend_data <- data.frame(
    Model = factor(rep(ordered_models[1], 3), levels = ordered_models), 
    Bias_Correction = factor(rep(names(model_data)[1], 3), levels = names(model_data)),
    legend_item = factor(c("0", "Rate present", "Model had ΔCICc < 2"), 
                        levels = c("0", "Rate present", "Model had ΔCICc < 2"))
  )
  
  # Create the main heatmap
  p_heatmap <- ggplot(plot_data, aes(x = Model, y = Bias_Correction)) +
    geom_tile(aes(fill = Prop_Significant), color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(Prop_Best > 0, sprintf("%.2f", Prop_Best), ""),
                  color = Prop_Significant > 0.5), 
              size = 3, fontface = "bold") +
    # Add invisible points for legend items
    geom_point(data = legend_data,
               aes(shape = legend_item), 
               size = 0, 
               na.rm = TRUE) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_shape_manual(
      values = c("0" = 22,  # Square with border
                 "Rate present" = 19,  # Filled circle
                 "Model had ΔCICc < 2" = 8),  # Star
      labels = c("0" = "0",
                 "Rate present" = "Rate present",
                 "Model had ΔCICc < 2" = bquote(paste("Model had ", Delta, "CICc < 2 in full dataset"))),
      name = "",
      guide = guide_legend(
        override.aes = list(
          size = c(5, 3, 3),
          fill = c("white", "black", NA),
          color = c("black", "black", "black")
        ),
        order = 2,
        label.theme = element_text(size = 9)
      )
    ) +
    scale_fill_gradient(
      low = "#E5F5FC",  # light blue
      high = "#001233", # Darkest blue possible
      limits = c(0.001, 1),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      name = bquote(paste("Proportion of Iterations (", Delta, "CICc < 2)")),
      na.value = "white",  # Cells with 0 iterations will be white
      guide = guide_colorbar(
        barwidth = 1.5,
        barheight = 10,
        label.theme = element_text(size = 9),
        title.theme = element_text(size = 10),
        title.position = "top",
        title.hjust = 0.5,
        order = 1
      )
    ) +
    scale_x_discrete(labels = NULL) +
    labs(
      title = title,
      subtitle = paste("Numbers show proportion of iterations where model was best."),
      x = NULL,
      y = "Downsampled Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
    )
  
  # Create data frame for delta CICc labels
  delta_labels <- data.frame(
    Model = factor(ordered_models, levels = ordered_models),
    delta_text = sapply(ordered_models, function(m) {
      if (!is.null(full_dataset_delta_cicc) && m %in% names(full_dataset_delta_cicc)) {
        sprintf("%.2f", full_dataset_delta_cicc[m])
      } else {
        "—"  # Em dash for models not in full dataset
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Create the rate presence/absence plot
  p_rates <- ggplot(rate_data, aes(x = Model, y = Rate)) +
    geom_point(data = filter(rate_data, Present), 
               shape = 19,  # Use circle for all points
               size = 3) +
    scale_x_discrete(labels = NULL) +
    labs(x = "Model (defined by rates included in model)", y = "Rate") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )
  
  # Create delta CICc labels plot with left label
  p_delta <- ggplot(delta_labels, aes(x = Model, y = 1)) +
    geom_text(aes(label = delta_text), size = 3.5, vjust = 0.5) +
    scale_x_discrete(labels = NULL, limits = ordered_models) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    # Add label to the left - will need to adjust coordinates
    annotate("text", x = 0.3, y = 1, 
             label = bquote(paste(Delta, "CICc from full dataset analysis:")), 
             hjust = 1, vjust = 0.5, size = 3.5) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL,
         caption = "") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 0, l = 80),
      plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10))
    )
  
  # Combine the plots using patchwork
  library(patchwork)
  p <- p_heatmap / p_delta / p_rates + 
    plot_layout(heights = c(length(model_data), 0.5, length(rates))) +
    # plot_annotation(
    #   caption = ""
    # ) &
    theme(
      legend.position = "right",
      plot.caption = element_text(size = 9, hjust = 0.5)
    )
  
  # Save the plot
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Calculate appropriate dimensions based on number of models and bias corrections
    plot_width <- max(20, 6 + length(all_models) * 0.3)
    plot_height <- max(6, 4 + length(model_data) * 0.4)
    
    ggsave(output_file, p, width = plot_width, height = plot_height, dpi = 600)
    ggsave(gsub(".pdf", ".png", output_file), p, 
           width = plot_width, height = plot_height, dpi = 300)
  }
  
  return(p)
}
