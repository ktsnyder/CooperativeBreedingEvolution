# create_downsampled_plots_complete.R
# Complete version of create_downsampled_plots_flexible that can replace create_downsampled_plots
# for all downsampling outputs including geographic and territoriality biases
# Kate Snyder
# Created: 2025-06-10

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)

#' Complete flexible version of create_downsampled_plots
#'
#' This version can fully replace create_downsampled_plots() for all types of downsampling:
#' - Geographic bias correction (Holarctic, Tropical, Global)
#' - Territoriality bias correction
#' - Dimorphism bias correction
#'
#' @param downsampling_results Output from run_multiple_phylopath or similar
#' @param detailed_models_input Either a dataframe, file path, or NULL
#' @param downsampling_info Information about the downsampling (for downsample_info plot)
#' @param output_prefix Prefix for output files
#' @param output_dir Output directory
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_downsampled_plots_complete <- function(downsampling_results,
                                             detailed_models_input = NULL,
                                             downsampling_info = NULL,
                                             output_prefix = "phylopath_downsampled",
                                             output_dir = "Outputs/PhylopathPlots",
                                             save_png = TRUE,
                                             save_pdf = TRUE) {
  
  # Ensure required functions are available
  if (!exists("aggregate_by_seed")) {
    stop("Required function aggregate_by_seed not found. Please source run_phylopath_fxns.R")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Downsampling information plot (for geographic/territoriality biases)
  if (!is.null(downsampling_info) && "proportions" %in% names(downsampling_info)) {
    cat("Creating downsampling information plot...\n")
    
    # Prepare data for plotting
    prop_data <- downsampling_info$proportions
    
    # Handle different column names (territoriality vs group)
    if ("territoriality" %in% colnames(prop_data)) {
      prop_data$group <- factor(prop_data$territoriality, levels = c("Low", "High"))
    } else if ("group" %in% colnames(prop_data)) {
      # Already has group column, just ensure it's a factor
      prop_data$group <- factor(prop_data$group)
    }
    
    # Create stacked bar plot
    prop_data_long <- prop_data %>%
      mutate(
        without_data = total_species - species_with_data
      ) %>%
      pivot_longer(cols = c(species_with_data, without_data),
                  names_to = "data_status",
                  values_to = "count")
    
    p_downsample_info <- ggplot(prop_data_long, 
                                aes(x = group, y = count, fill = data_status)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
      geom_text(data = prop_data,
               aes(x = group, y = total_species + 50, 
                   label = paste0("n = ", total_species, "\n",
                                 round(proportion * 100, 1), "% with data")),
               inherit.aes = FALSE, size = 4) +
      scale_fill_manual(values = c("species_with_data" = "#2E86AB", 
                                 "without_data" = "#F4A261"),  # Orange instead of gray
                       labels = c("species_with_data" = "With FS & CB data",
                                 "without_data" = "Missing data"),
                       name = "") +
      labs(
        title = "Data Availability by Group",
        subtitle = if (!is.null(downsampling_info$target_proportion)) {
          paste("Target proportion:", round(downsampling_info$target_proportion, 3),
                "\nRemove", downsampling_info$n_to_remove, "species from",
                downsampling_info$downsample_info$group_to_downsample)
        } else {
          NULL
        },
        x = "Group",
        y = "Number of Species"
      ) +
      theme_cowplot(12) +
      theme(legend.position = "bottom")
  } else {
    p_downsample_info <- NULL
  }
  
  # 2. Handle detailed models input
  detailed_models_df <- NULL
  
  if (is.data.frame(detailed_models_input)) {
    # Case A: Already a dataframe
    cat("Using provided detailed models dataframe\n")
    detailed_models_df <- detailed_models_input
    
  } else if (is.character(detailed_models_input) && file.exists(detailed_models_input)) {
    # Case B: File path provided
    cat("Loading detailed models from:", detailed_models_input, "\n")
    detailed_models_df <- read.csv(detailed_models_input)
    
  } else if (is.null(detailed_models_input)) {
    # Case C: Search for file using pattern
    cat("Searching for detailed models file...\n")
    
    # First try the standard pattern
    detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", 
                                                        gsub(" ", "_", output_prefix), 
                                                        ".*\\.csv$"))
    
    # If no files found, try a more general pattern
    if (length(detailed_models_files) == 0) {
      # Extract the key part of the prefix (e.g., "Remove66HighTerr")
      prefix_parts <- strsplit(output_prefix, " ")[[1]]
      key_pattern <- grep("Remove|High", prefix_parts, value = TRUE)
      if (length(key_pattern) > 0) {
        pattern_str <- paste(key_pattern, collapse = ".*")
        detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", pattern_str, ".*\\.csv$"))
      }
    }
    
    # If still no files, try searching in output directory
    if (length(detailed_models_files) == 0 && output_dir != ".") {
      detailed_models_files <- list.files(output_dir, 
                                         pattern = paste0("detailed_models_.*", 
                                                         gsub(" ", "_", output_prefix), 
                                                         ".*\\.csv$"),
                                         full.names = TRUE)
    }
    
    if (length(detailed_models_files) > 0) {
      detailed_models_file <- detailed_models_files[1]
      cat("Loading detailed models from:", detailed_models_file, "\n")
      detailed_models_df <- read.csv(detailed_models_file)
    } else {
      cat("No detailed models file found\n")
    }
    
  } else if (!is.null(downsampling_results$detailed_models)) {
    # Check if detailed_models is included in the results object
    cat("Using detailed models from results object\n")
    detailed_models_df <- downsampling_results$detailed_models
  }
  
  # 3. Process detailed models if available
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    # Extract edge information
    edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Aggregate by seed
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format for plotting
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # 3a. Create violin plots
    cat("Creating violin plots of path coefficients...\n")
    
    # Clean path names for display
    seed_level_conditional$path_clean <- gsub("_to_", " → ", seed_level_conditional$path)
    seed_level_conditional$path_clean <- gsub("FemaleSong_Agg01", "Female Song", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("HighConfidence_Coop", "Cooperation", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("Territory_12vs3", "Territory Type", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMass_AVONET", "Body Mass", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMaleFemalePlumageDiffAbs", "Plumage Dimorphism", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("PercentAbsLogWingDimorphism", "Wing Dimorphism", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("PercentLogWingDimorphism_AVONET", "Wing Dimorphism", seed_level_conditional$path_clean)
    
    # Filter to paths that appear in multiple seeds
    path_counts <- seed_level_conditional %>%
      group_by(path) %>%
      summarise(n = n()) %>%
      filter(n > 10)
    
    seed_level_plot <- seed_level_conditional %>%
      filter(path %in% path_counts$path)
    
    if (nrow(seed_level_plot) > 0) {
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
    } else {
      cat("Not enough data for violin plots\n")
      p_violin <- NULL
    }
    
    # 3b. Create heatmap
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
        x <- gsub("logMaleFemalePlumageDiffAbs", "Plumage\nDimorphism", x)
        x <- gsub("PercentAbsLogWingDimorphism", "Wing\nDimorphism", x)
        x <- gsub("PercentLogWingDimorphism_AVONET", "Wing\nDimorphism", x)
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
        geom_text(aes(label = paste0(round(Coefficient, 3), 
                                    ifelse(!is.na(sig_label), sig_label, ""))),
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
    
    # 3c. Model frequency bar plot
    cat("Creating model frequency plot...\n")
    if (!is.null(downsampling_results$model_frequencies) && nrow(downsampling_results$model_frequencies) > 0) {
      model_freq <- downsampling_results$model_frequencies %>%
        arrange(desc(frequency_in_sub2)) %>%
        head(15) %>%
        mutate(
          prop_sub2 = frequency_in_sub2 / length(unique(detailed_models_df$seed))
        )
      
      p_model_freq <- ggplot(model_freq, 
                            aes(x = reorder(model_name, frequency_in_sub2),
                                y = frequency_in_sub2)) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                    round(prop_sub2 * 100, 1), "%)")),
                 hjust = -0.1, size = 3) +
        coord_flip() +
        labs(
          title = "Top Models Across Downsampling Iterations",
          subtitle = paste("Models appearing in Δ CICc < 2 set across",
                          length(unique(detailed_models_df$seed)), 
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
    cat("No detailed models data available for plotting\n")
    p_violin <- NULL
    p_heatmap <- NULL
    p_model_freq <- NULL
  }
  
  # 4. Save all plots with proper file names
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
  
  if (save_pdf) {
    # Save individual PDFs for each plot
    if (!is.null(p_downsample_info)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " downsample_info.pdf")),
             p_downsample_info, width = 8, height = 6, device = "pdf")
    }
    if (!is.null(p_violin)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " violin.pdf")),
             p_violin, width = 10, height = 8, device = "pdf")
    }
    if (!is.null(p_heatmap)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " heatmap.pdf")),
             p_heatmap, width = 8, height = 8, device = "pdf")
    }
    if (!is.null(p_model_freq)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.pdf")),
             p_model_freq, width = 10, height = 8, device = "pdf")
    }
    
    cat("PDF files saved to:", output_dir, "\n")
  }
  
  return(list(
    downsample_info = p_downsample_info,
    violin = p_violin,
    heatmap = p_heatmap,
    model_freq = p_model_freq
  ))
}

# Test function to verify it works with different input types
test_create_downsampled_plots_complete <- function() {
  # Create dummy data for testing
  dummy_detailed_models <- data.frame(
    seed = rep(1:3, each = 2),
    model = rep(c("Model1", "Model2"), 3),
    CICc = rnorm(6, 100, 5),
    delta_CICc = c(0, 0.5, 0, 1.5, 0, 1),
    nSpecies = 100,
    FS_to_CB_est = rnorm(6, 0.2, 0.1),
    FS_to_CB_se = abs(rnorm(6, 0.05, 0.01)),
    FS_to_CB_p = runif(6, 0, 0.1)
  )
  
  dummy_results <- list(
    model_frequencies = data.frame(
      model_name = c("Model1", "Model2"),
      frequency_in_sub2 = c(3, 2)
    ),
    detailed_models = dummy_detailed_models
  )
  
  # Test Case A: Dataframe input
  cat("Testing with dataframe input...\n")
  test1 <- create_downsampled_plots_complete(
    downsampling_results = dummy_results,
    detailed_models_input = dummy_detailed_models,
    output_prefix = "test_dataframe",
    save_png = FALSE
  )
  
  # Test Case B: File path input
  cat("\nTesting with file path input...\n")
  temp_file <- tempfile(fileext = ".csv")
  write.csv(dummy_detailed_models, temp_file, row.names = FALSE)
  test2 <- create_downsampled_plots_complete(
    downsampling_results = dummy_results,
    detailed_models_input = temp_file,
    output_prefix = "test_filepath",
    save_png = FALSE
  )
  unlink(temp_file)
  
  # Test Case C: NULL input (will search for file)
  cat("\nTesting with NULL input (file search)...\n")
  test3 <- create_downsampled_plots_complete(
    downsampling_results = dummy_results,
    detailed_models_input = NULL,
    output_prefix = "test_null",
    save_png = FALSE
  )
  
  cat("\nAll tests completed\n")
  return(list(test1 = test1, test2 = test2, test3 = test3))
}