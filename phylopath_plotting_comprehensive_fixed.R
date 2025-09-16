# phylopath_plotting_comprehensive_fixed.R
# Fixed version with proper file naming, DAG layouts, and colors
# Kate Snyder
# Created: 2025-06-08

library(phylopath)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(viridis)

#' Create all plots for non-downsampled phylopath runs
#'
#' @param phylopath_output Output from run_CB_FS_Terr_phylopath()
#' @param output_prefix Prefix for output files - should include all variables
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_nondownsampled_plots <- function(phylopath_output, 
                                       output_prefix = "phylopath",
                                       output_dir = "Outputs/PhylopathPlots",
                                       save_png = TRUE,
                                       save_pdf = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract components
  result <- phylopath_output$result
  summary_result <- summary(result)
  var_map <- phylopath_output$var_map
  nSpecies <- phylopath_output$nSpecies
  
  # Get variable names for file naming
  var_names <- names(var_map)
  var_string <- paste(unlist(var_map), collapse = " ")
  
  # Create file name base with all variables and species count
  file_base <- paste0("phylopath ", var_string, " ", nSpecies, "species")
  
  # Set up layout positions for DAGs with better spacing
  if (!is.null(var_map)) {
    # Create sensible positions based on variable names
    var_names <- names(var_map)
    n_vars <- length(var_names)
    
    if (n_vars == 4) {
      # Standard 4-variable layout with good spacing
      phylopath_map_positions <- data.frame(
        name = unlist(var_map),
        x = c(2, 8, 5, 5),  # Original positions
        y = c(9, 9, 1, 5),  # Original positions
        stringsAsFactors = FALSE
      )
    } else {
      # Circular layout for other numbers
      angles <- seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars]
      phylopath_map_positions <- data.frame(
        name = unlist(var_map),
        x = 5 + 3 * cos(angles),  # Smaller radius to keep away from edges
        y = 5 + 3 * sin(angles),
        stringsAsFactors = FALSE
      )
    }
  } else {
    phylopath_map_positions <- NULL
  }
  
  # 1. Best model DAG
  cat("Creating best model DAG...\n")
  best_model <- best(result)
  p_best <- plot(best_model, 
                 manual_layout = phylopath_map_positions,
                 text_size = 3,  # Reduced from 4
                 box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Best Model (Lowest CICc)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(30, 30, 30, 30)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 2. Fully-averaged best models DAG
  cat("Creating fully-averaged models DAG...\n")
  full_avg <- average(result, cut_off = 2, avg_method = "full")
  p_full_avg <- plot(full_avg,
                     manual_layout = phylopath_map_positions,
                     text_size = 3,  # Reduced from 4
                     box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Fully-Averaged Models (Δ CICc < 2)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(30, 30, 30, 30)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 3. Conditional-averaged best models DAG
  cat("Creating conditional-averaged models DAG...\n")
  cond_avg <- average(result, cut_off = 2, avg_method = "conditional")
  p_cond_avg <- plot(cond_avg,
                     manual_layout = phylopath_map_positions,
                     text_size = 3,  # Reduced from 4
                     box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Conditional-Averaged Models (Δ CICc < 2)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(30, 30, 30, 30)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 4. Bar plot of CICc values - KEEP ORIGINAL MODEL NAMES
  cat("Creating CICc bar plot...\n")
  # Prepare data - DON'T clean the model names
  model_data <- summary_result %>%
    mutate(
      is_top = delta_CICc < 2
    ) %>%
    arrange(CICc)
  
  # Limit to top 20 models for clarity
  n_models <- min(20, nrow(model_data))
  model_data_plot <- model_data[1:n_models, ]
  
  p_cicbar <- ggplot(model_data_plot, 
                     aes(x = reorder(model, -CICc),  # Use original model names
                         y = CICc,
                         fill = is_top)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = round(delta_CICc, 2)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#A7A7A7"),
                      labels = c("TRUE" = "Δ CICc < 2", "FALSE" = "Δ CICc ≥ 2"),
                      name = "Model Set") +
    labs(
      title = "Model Comparison by CICc",
      subtitle = paste("Showing top", n_models, "models. Values show Δ CICc"),
      x = NULL,
      y = "CICc"
    ) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  # 5. Plot of all tested models
  cat("Creating summary plot of all models...\n")
  p_summary <- plot(summary_result) +
    ggtitle(paste("All Tested Models (n =", nrow(summary_result), ")")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  # Combine DAGs into one figure
  dag_combined <- plot_grid(
    p_best, p_full_avg, p_cond_avg,
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 16
  )
  
  # Save plots with proper file names
  if (save_png) {
    # Individual plots - increased sizes to prevent squishing
    ggsave(file.path(output_dir, paste0(file_base, " best_model.png")), 
           p_best, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " full_avg.png")), 
           p_full_avg, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " cond_avg.png")), 
           p_cond_avg, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " cicbar.png")), 
           p_cicbar, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " summary.png")), 
           p_summary, width = 12, height = 10, dpi = 300)
    
    # Combined DAG figure - increased height
    ggsave(file.path(output_dir, paste0(file_base, " dags_combined.png")), 
           dag_combined, width = 38, height = 10, dpi = 300)
    
    cat("PNG files saved to:", output_dir, "\n")
  }
  
  if (save_pdf) {
    pdf(file.path(output_dir, paste0(file_base, " all_plots.pdf")), 
        width = 15, height = 10)
    print(p_best)
    print(p_full_avg)
    print(p_cond_avg)
    print(p_cicbar)
    print(p_summary)
    dev.off()
    cat("PDF saved to:", file.path(output_dir, paste0(file_base, " all_plots.pdf")), "\n")
  }
  
  return(list(
    best_model = p_best,
    full_avg = p_full_avg,
    cond_avg = p_cond_avg,
    cicbar = p_cicbar,
    summary = p_summary,
    dags_combined = dag_combined
  ))
}

#' Create all plots for downsampled phylopath runs
#'
#' @param downsampling_results Output from run_multiple_phylopath()
#' @param downsampling_info Information about the downsampling
#' @param output_prefix Prefix including all variables and Remove[n][group] format
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_downsampled_plots <- function(downsampling_results,
                                   downsampling_info = NULL,
                                   output_prefix = "phylopath_downsampled",
                                   output_dir = "Outputs/PhylopathPlots",
                                   save_png = TRUE,
                                   save_pdf = TRUE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Species counts and fractions plot
  # if (!is.null(downsampling_info)) {
  #   cat("Creating downsampling information plot...\n")
  #   
  #   # Prepare data for plotting
  #   if ("proportions" %in% names(downsampling_info)) {
  #     prop_data <- downsampling_info$proportions
  #     prop_data$group <- factor(prop_data$territoriality, levels = c("Low", "High"))
  #     
  #     # Create stacked bar plot
  #     prop_data_long <- prop_data %>%
  #       mutate(
  #         without_data = total_species - species_with_data
  #       ) %>%
  #       pivot_longer(cols = c(species_with_data, without_data),
  #                   names_to = "data_status",
  #                   values_to = "count")
  #     
  #     p_downsample_info <- ggplot(prop_data_long, 
  #                                 aes(x = group, y = count, fill = data_status)) +
  #       geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
  #       geom_text(data = prop_data,
  #                aes(x = group, y = total_species + 50, 
  #                    label = paste0("n = ", total_species, "\n",
  #                                  round(proportion * 100, 1), "% with data")),
  #                inherit.aes = FALSE, size = 4) +
  #       scale_fill_manual(values = c("species_with_data" = "#2E86AB", 
  #                                  "without_data" = "#F4A261"),  # Changed from gray
  #                        labels = c("species_with_data" = "With FS & CB data",
  #                                  "without_data" = "Missing data"),
  #                        name = "") +
  #       labs(
  #         title = "Data Availability by Territoriality",
  #         subtitle = paste("Target proportion:", round(downsampling_info$target_proportion, 3),
  #                         "\nRemove", downsampling_info$n_to_remove, "species from",
  #                         downsampling_info$downsample_info$group_to_downsample),
  #         x = "Territoriality",
  #         y = "Number of Species"
  #       ) +
  #       theme_cowplot(12) +
  #       theme(legend.position = "bottom")
  #   } else {
  #     p_downsample_info <- NULL
  #   }
  # } else {
     p_downsample_info <- NULL
  # }
  
  # 2. Process detailed models if available
  detailed_models_file <- list.files(pattern = paste0("detailed_models_.*", 
                                                     gsub(" ", "_", output_prefix), 
                                                     ".*\\.csv$"))[1]
  
  if (!is.na(detailed_models_file)) {
    cat("Loading detailed models from:", detailed_models_file, "\n")
    detailed_models_df <- read.csv(detailed_models_file)
    
    # Extract edge information
    edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Aggregate by seed using the aggregate_by_seed function
    source("run_phylopath_fxns.R")  # Make sure we have the aggregation functions
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format for plotting
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # 3. Violin plots of path coefficients
    cat("Creating violin plots of path coefficients...\n")
    
    # Clean path names for display
    seed_level_conditional$path_clean <- gsub("_to_", " → ", seed_level_conditional$path)
    seed_level_conditional$path_clean <- gsub("FemaleSong_Agg01", "Female Song", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("HighConfidence_Coop", "Cooperation", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("Territory_12vs3", "Territory Type", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMass_AVONET", "Body Mass", seed_level_conditional$path_clean)
    
    # Filter to paths that appear in multiple seeds
    path_counts <- seed_level_conditional %>%
      group_by(path) %>%
      summarise(n = n()) %>%
      filter(n > 10)
    
    seed_level_plot <- seed_level_conditional %>%
      filter(path %in% path_counts$path)
    
    p_violin <- ggplot(seed_level_plot, 
                      aes(x = path_clean, y = mean_coefficient)) +
      geom_violin(fill = "lightblue", alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      coord_flip() +
      labs(
        title = "Distribution of Path Coefficients Across Iterations",
        subtitle = paste("Based on", length(unique(seed_level_plot$seed)), 
                        "downsampling iterations (conditional averaging)"),
        x = NULL,
        y = "Mean Coefficient"
      ) +
      theme_cowplot(12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10)
      )
    
    # 4. Heatmap of mean coefficients
    cat("Creating heatmap of mean coefficients...\n")
    
    # Get the across-seeds summary
    path_summary <- agg_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
    
    # Create from-to matrix for heatmap
    if (nrow(path_summary) > 0) {
      # Create a matrix of coefficients
      unique_from <- unique(path_summary$from)
      unique_to <- unique(path_summary$to)
      
      coef_matrix <- matrix(NA, 
                           nrow = length(unique_from), 
                           ncol = length(unique_to),
                           dimnames = list(unique_from, unique_to))
      
      for (i in 1:nrow(path_summary)) {
        coef_matrix[path_summary$from[i], path_summary$to[i]] <- path_summary$mean_coef_conditional[i]
      }
      
      # Convert to long format for ggplot
      coef_long <- as.data.frame(as.table(coef_matrix))
      names(coef_long) <- c("From", "To", "Coefficient")
      coef_long <- coef_long[!is.na(coef_long$Coefficient), ]
      
      # Clean variable names
      clean_names <- function(x) {
        x <- gsub("FemaleSong_Agg01", "Female\nSong", x)
        x <- gsub("HighConfidence_Coop", "Cooperation", x)
        x <- gsub("TerritorialityWeakVsStrong", "Territoriality", x)
        x <- gsub("Territory_12vs3", "Territory\nType", x)
        x <- gsub("logMass_AVONET", "Body Mass", x)
        return(x)
      }
      
      coef_long$From <- clean_names(coef_long$From)
      coef_long$To <- clean_names(coef_long$To)
      
      # Add significance info if available
      sig_info <- path_summary %>%
        mutate(
          From_clean = clean_names(from),
          To_clean = clean_names(to),
          sig_label = ifelse(mean_prop_sig_conditional > 0.95, "***",
                           ifelse(mean_prop_sig_conditional > 0.8, "**",
                                 ifelse(mean_prop_sig_conditional > 0.5, "*", "")))
        )
      
      coef_long <- coef_long %>%
        left_join(sig_info, by = c("From" = "From_clean", "To" = "To_clean"))
      
      p_heatmap <- ggplot(coef_long, aes(x = From, y = To, fill = Coefficient)) +
        geom_tile(color = "white", size = 0.5) +
        geom_text(aes(label = paste0(round(Coefficient, 3), sig_label)),
                 color = ifelse(abs(coef_long$Coefficient) > 0.3, "white", "black"),
                 size = 4) +
        scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2E86AB",
                           midpoint = 0, name = "Mean\nCoefficient",
                           limits = c(-max(abs(coef_long$Coefficient)), 
                                     max(abs(coef_long$Coefficient)))) +
        labs(
          title = "Mean Path Coefficients Across Iterations",
          subtitle = "Significance: *** >95%, ** >80%, * >50% of iterations",
          x = "From Variable",
          y = "To Variable"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    } else {
      p_heatmap <- NULL
    }
    
    # 5. Model frequency bar plot - KEEP ORIGINAL MODEL NAMES
    cat("Creating model frequency plot...\n")
    if (!is.null(downsampling_results$model_frequencies)) {
      model_freq <- downsampling_results$model_frequencies %>%
        arrange(desc(frequency_in_sub2)) %>%
        head(15) %>%
        mutate(
          prop_sub2 = frequency_in_sub2 / max(downsampling_results$results_df$seed)
        )
      
      p_model_freq <- ggplot(model_freq, 
                            aes(x = reorder(model_name, frequency_in_sub2),  # Use original model names
                                y = frequency_in_sub2)) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                    round(prop_sub2 * 100, 1), "%)")),
                 hjust = -0.1, size = 3) +
        coord_flip() +
        labs(
          title = "Top Models Across Downsampling Iterations",
          subtitle = paste("Models appearing in Δ CICc < 2 set across",
                          length(unique(downsampling_results$results_df$seed)), 
                          "iterations"),
          x = NULL,
          y = "Frequency in Top Model Set"
        ) +
        theme_cowplot(12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 9)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    } else {
      p_model_freq <- NULL
    }
    
  } else {
    cat("No detailed models file found\n")
    p_violin <- NULL
    p_heatmap <- NULL
    p_model_freq <- NULL
  }
  
  # Save all plots with proper file names
  if (save_png) {
    if (!is.null(p_downsample_info)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " downsample_info.png")),
             p_downsample_info, width = 8, height = 6, dpi = 300)
    }
    if (!is.null(p_violin)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " violin.png")),
             p_violin, width = 10, height = 8, dpi = 300)
    }
    if (!is.null(p_heatmap)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " heatmap.png")),
             p_heatmap, width = 8, height = 8, dpi = 300)
    }
    if (!is.null(p_model_freq)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.png")),
             p_model_freq, width = 10, height = 8, dpi = 300)
    }
    
    cat("PNG files saved to:", output_dir, "\n")
  }
  
  return(list(
    downsample_info = p_downsample_info,
    violin = p_violin,
    heatmap = p_heatmap,
    model_freq = p_model_freq
  ))
}

#' Wrapper function to create all phylopath plots
#'
#' @param analysis_type Either "nondownsampled" or "downsampled"
#' @param phylopath_output Output from phylopath analysis
#' @param downsampling_info Optional downsampling information
#' @param output_prefix Prefix for output files
#' @param save_png Save as PNG
#' @param save_pdf Save as PDF
create_all_phylopath_plots <- function(analysis_type = c("nondownsampled", "downsampled"),
                                      phylopath_output,
                                      downsampling_info = NULL,
                                      output_prefix = "phylopath",
                                      output_dir = "Outputs/PhylopathPlots",
                                      save_png = TRUE,
                                      save_pdf = FALSE) {
  
  analysis_type <- match.arg(analysis_type)
  
  if (analysis_type == "nondownsampled") {
    plots <- create_nondownsampled_plots(
      phylopath_output = phylopath_output,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = save_png,
      save_pdf = save_pdf
    )
  } else {
    plots <- create_downsampled_plots(
      downsampling_results = phylopath_output,
      downsampling_info = downsampling_info,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = save_png,
      save_pdf = save_pdf
    )
  }
  
  return(plots)
}
