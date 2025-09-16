# create_phylopath_plots_enhanced_fixed.R
# Fixed version that keeps original model names (no renaming)
# Kate Snyder
# Created: 2025-06-08

library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)

#' Create enhanced plots from detailed_models CSV files
#'
#' @param csv_file Path to detailed_models CSV file
#' @param output_prefix Prefix for output files
#' @param save_png Whether to save as PNG
#' @param save_pdf Whether to save as PDF
create_enhanced_phylopath_plots <- function(csv_file, 
                                           output_prefix = NULL,
                                           save_png = TRUE,
                                           save_pdf = TRUE) {
  
  # Load data
  detailed_models_df <- read.csv(csv_file)
  
  if (is.null(output_prefix)) {
    output_prefix <- gsub("\\.csv$", "", basename(csv_file))
  }
  
  # 1. Model frequency analysis
  model_freq <- detailed_models_df %>%
    group_by(model) %>%
    summarise(
      frequency = n(),
      mean_CICc = mean(CICc),
      mean_delta_CICc = mean(delta_CICc),
      .groups = "drop"
    ) %>%
    arrange(desc(frequency))
  
  # Top models plot - USE ORIGINAL MODEL NAMES
  top_n <- min(15, nrow(model_freq))
  p_models <- ggplot(model_freq[1:top_n,], 
                     aes(x = reorder(model, frequency), y = frequency)) +
    geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
    geom_text(aes(label = paste0(frequency, " (", round(mean_delta_CICc, 2), ")")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Top Models by Frequency",
      subtitle = paste("Based on", nrow(detailed_models_df), "model selections across",
                      length(unique(detailed_models_df$seed)), "iterations"),
      x = NULL,
      y = "Frequency (mean Δ CICc shown)"
    ) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_line(color = "gray90", size = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  # 2. Edge coefficient analysis
  # Extract all edge columns
  edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
  
  edge_data <- data.frame()
  for (col in edge_cols) {
    edge_name <- gsub("_est$", "", col)
    p_col <- gsub("_est$", "_p", col)
    
    if (p_col %in% names(detailed_models_df)) {
      valid_rows <- !is.na(detailed_models_df[[col]])
      
      if (sum(valid_rows) > 0) {
        temp_df <- data.frame(
          edge = edge_name,
          coefficient = detailed_models_df[[col]][valid_rows],
          p_value = detailed_models_df[[p_col]][valid_rows],
          seed = detailed_models_df$seed[valid_rows],
          model = detailed_models_df$model[valid_rows]
        )
        edge_data <- rbind(edge_data, temp_df)
      }
    }
  }
  
  # Clean edge names for display (but not model names)
  edge_data$edge_clean <- gsub("_to_", " → ", edge_data$edge)
  edge_data$edge_clean <- gsub("FemaleSong_Agg01", "Female Song", edge_data$edge_clean)
  edge_data$edge_clean <- gsub("HighConfidence_Coop", "Cooperation", edge_data$edge_clean)
  edge_data$edge_clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", edge_data$edge_clean)
  edge_data$edge_clean <- gsub("Territory_12vs3", "Territory Type", edge_data$edge_clean)
  edge_data$edge_clean <- gsub("logMass_AVONET", "Body Mass", edge_data$edge_clean)
  
  # Edge summary
  edge_summary <- edge_data %>%
    group_by(edge, edge_clean) %>%
    summarise(
      n_occurrences = n(),
      mean_coef = mean(coefficient),
      median_coef = median(coefficient),
      prop_significant = mean(p_value < 0.05),
      prop_positive = mean(coefficient > 0),
      .groups = "drop"
    ) %>%
    filter(n_occurrences > 10) %>%  # Only show edges that appear in multiple models
    arrange(desc(n_occurrences))
  
  # Edge coefficient plot
  p_edges <- ggplot(edge_summary, 
                     aes(x = reorder(edge_clean, n_occurrences), y = mean_coef)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = mean_coef - 0.1, ymax = mean_coef + 0.1),
                  width = 0.2, color = "gray60") +
    geom_point(aes(size = n_occurrences, 
                   color = prop_significant,
                   shape = ifelse(prop_positive > 0.5, "Positive", "Negative")),
               alpha = 0.8) +
    coord_flip() +
    labs(
      title = "Edge Coefficients Across Models",
      subtitle = "Size shows frequency, color shows proportion significant",
      x = NULL,
      y = "Mean Coefficient",
      size = "Frequency",
      color = "Prop. Significant",
      shape = "Direction"
    ) +
    scale_color_gradient2(low = "gray70", mid = "#F77F00", high = "#E63946",
                         midpoint = 0.5, limits = c(0, 1)) +
    scale_shape_manual(values = c("Positive" = 16, "Negative" = 17)) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "right",
      legend.box = "vertical"
    )
  
  # 3. Network visualization of most common model structure
  if (nrow(model_freq) > 0) {
    top_model_name <- model_freq$model[1]
    top_model_data <- detailed_models_df %>%
      filter(model == top_model_name)
    
    # Extract edges for this model
    top_model_edges <- edge_data %>%
      filter(model == top_model_name) %>%
      group_by(edge, edge_clean) %>%
      summarise(
        mean_coef = mean(coefficient),
        mean_p = mean(p_value),
        .groups = "drop"
      )
    
    # Create simple network plot
    if (nrow(top_model_edges) > 0) {
      p_network <- ggplot(top_model_edges, aes(x = 0, y = 0)) +
        geom_text(aes(label = paste0(edge_clean, "\n", 
                                    "β = ", round(mean_coef, 3),
                                    "\np = ", round(mean_p, 3))),
                  size = 3.5) +
        labs(
          title = paste("Most Frequent Model:", model_freq$model[1]),  # Use original model name
          subtitle = paste("Appeared in", model_freq$frequency[1], "iterations")
        ) +
        theme_void() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5)
        )
    } else {
      p_network <- NULL
    }
  } else {
    p_network <- NULL
  }
  
  # 4. Create combined plot
  plots_list <- list(p_models, p_edges)
  if (!is.null(p_network)) {
    plots_list <- c(plots_list, list(p_network))
  }
  
  # Save plots
  if (save_png) {
    # Individual plots
    png(paste0(output_prefix, "_models.png"), width = 10, height = 8, units = "in", res = 300)
    print(p_models)
    dev.off()
    
    png(paste0(output_prefix, "_edges.png"), width = 10, height = 8, units = "in", res = 300)
    print(p_edges)
    dev.off()
    
    # Combined plot
    combined_plot <- plot_grid(plotlist = plots_list[1:2], 
                              ncol = 1, 
                              rel_heights = c(1, 1),
                              labels = c("A", "B"),
                              label_size = 16)
    
    png(paste0(output_prefix, "_combined.png"), width = 12, height = 16, units = "in", res = 300)
    print(combined_plot)
    dev.off()
    
    cat("PNG plots saved:\n")
    cat("  -", paste0(output_prefix, "_models.png"), "\n")
    cat("  -", paste0(output_prefix, "_edges.png"), "\n")
    cat("  -", paste0(output_prefix, "_combined.png"), "\n")
  }
  
  if (save_pdf) {
    pdf(paste0(output_prefix, "_all_plots.pdf"), width = 10, height = 8)
    print(p_models)
    print(p_edges)
    if (!is.null(p_network)) print(p_network)
    dev.off()
    
    cat("PDF saved:", paste0(output_prefix, "_all_plots.pdf"), "\n")
  }
  
  # Return plots and summaries
  return(list(
    plots = list(
      models = p_models,
      edges = p_edges,
      network = p_network
    ),
    summaries = list(
      model_freq = model_freq,
      edge_summary = edge_summary
    )
  ))
}

# Example usage function
process_all_detailed_models <- function(pattern = "detailed_models_.*\\.csv$",
                                       output_dir = "Outputs/PhylopathPlots") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Find all detailed model files
  files <- list.files(pattern = pattern, full.names = TRUE)
  
  cat("Found", length(files), "detailed model files to process\n")
  
  for (file in files) {
    cat("\nProcessing:", basename(file), "\n")
    
    # Extract descriptive name
    base_name <- gsub("detailed_models_", "", basename(file))
    base_name <- gsub("\\.csv$", "", base_name)
    
    output_prefix <- file.path(output_dir, base_name)
    
    # Create plots
    results <- create_enhanced_phylopath_plots(
      csv_file = file,
      output_prefix = output_prefix,
      save_png = TRUE,
      save_pdf = TRUE
    )
  }
  
  cat("\nAll files processed!\n")
}