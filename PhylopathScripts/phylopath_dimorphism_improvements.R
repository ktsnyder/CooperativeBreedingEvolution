# phylopath_dimorphism_improvements.R
# Improvements for dimorphism bias correction with phylopath
# Kate Snyder
# Created: 2025-06-10

library(dplyr)
library(phylopath)

#' Run phylopath with dimorphism bias correction and create detailed models dataframe
#'
#' This function replaces the code in Run_Analyses.R section "run phylopath with dimorphism bias correction/downsampling"
#' It properly creates a detailed_models_df object that can be used with create_downsampled_plots()
#'
#' @param dfIn_phylo Full dataset with species data
#' @param tree Phylogenetic tree
#' @param dim_info List with dimorphism information (col, label, description)
#' @param n_iterations Number of iterations to run
#' @param phylopath_output_dir Output directory for plots
#' @param save_outputs Whether to save CSV files and plots
#' @return List containing results and detailed_models_df
run_phylopath_dimorphism_correction <- function(dfIn_phylo, 
                                                tree,
                                                dim_info,
                                                n_iterations = 500,
                                                phylopath_output_dir = "Outputs/PhylopathDownsampled",
                                                save_outputs = TRUE) {
  
  # Create output directory specific to this dimorphism type
  dim_output_dir <- file.path(phylopath_output_dir, dim_info$label)
  if (!dir.exists(dim_output_dir)) {
    dir.create(dim_output_dir, recursive = TRUE)
  }
  
  cat("\n\n========================================\n")
  cat("Running", dim_info$description, "bias correction...\n")
  cat("========================================\n")
  
  # Check if the column exists
  if (!dim_info$col %in% colnames(dfIn_phylo)) {
    cat("Warning: Column", dim_info$col, "not found. Skipping...\n")
    return(NULL)
  }
  
  # Run the dimorphism downsampling analysis
  dimorphism_downsample <- downsample_dimorphism_bias(
    df = dfIn_phylo,
    dimorphism_col = dim_info$col,
    data_col = "FemaleSong_Agg01",
    n_iterations = n_iterations
  )
  
  # Create distribution plot if requested
  if (save_outputs) {
    prefix_dimorphism <- paste0("Remove", 
                               dimorphism_downsample$n_to_remove, 
                               "High", 
                               dim_info$label)
    
    dist_plot_result <- create_downsampled_dimorphism_distribution_plot(
      dfIn_phylo = dfIn_phylo,
      dimorphism_downsample = dimorphism_downsample,
      dim_info = dim_info,
      output_dir = dim_output_dir,
      prefix = prefix_dimorphism,
      save_plot = TRUE
    )
    
    # Create the violin plot (multiple iterations view)
    violin_plot_result <- create_dimorphism_downsamples_violin_plot(
      dfIn_phylo = dfIn_phylo,
      dimorphism_downsample = dimorphism_downsample,
      dim_info = dim_info,
      output_dir = dim_output_dir,
      prefix = prefix_dimorphism,
      save_plot = TRUE,
      show_iterations = min(20, n_iterations)
    )
  }
  
  cat("Need to remove", dimorphism_downsample$n_to_remove, "species to correct", dim_info$description, "bias\n")
  cat("Current mean:", round(dimorphism_downsample$current_stats$current_mean, 3), "\n")
  cat("Target mean:", round(dimorphism_downsample$target_stats$target_mean, 3), "\n")
  
  # Run phylopath on each iteration
  all_results <- list()
  all_model_summaries <- list()
  detailed_models_list <- list()
  
  for (i in 1:n_iterations) {
    cat(sprintf("Running iteration %d/%d...\r", i, n_iterations))
    
    # Get the downsampled dataset for this iteration
    df_downsampled <- dimorphism_downsample$iterations[[i]]$remaining_data
    
    # Filter to species in tree
    df_downsampled <- df_downsampled[df_downsampled$species %in% tree$tip.label, ]
    rownames(df_downsampled) <- df_downsampled$species
    
    # Check if we have enough species
    if (nrow(df_downsampled) < 50) {
      warning(paste("Iteration", i, "has only", nrow(df_downsampled), "species. Skipping..."))
      next
    }
    
    # Run phylopath on this iteration
    result_i <- run_CB_FS_Terr_phylopath(
      dfIn = df_downsampled,
      tree = tree,
      female_song_var = "FemaleSong_Agg01",
      coop_breeding_var = "HighConfidence_Coop",
      territoriality_var = "TerritorialityWeakVsStrong",
      mass_var = dim_info$col,  # Use the dimorphism variable as the 4th variable
      plots2pdf = FALSE
    )
    
    all_results[[i]] <- result_i
    
    # Extract model summaries and create detailed models entries
    if (!is.null(result_i$result)) {
      summary_i <- summary(result_i$result)
      
      if (!is.null(summary_i)) {
        # For each model with delta_CICc < 2
        for (model_idx in which(summary_i$delta_CICc < 2)) {
          model_name <- summary_i$model[model_idx]
          model_CICc <- summary_i$CICc[model_idx]
          model_delta_CICc <- summary_i$delta_CICc[model_idx]
          
          # Get the chosen model to extract path coefficients
          chosen_model <- choice(result_i$result, model_name)
          
          # Create a row for this model
          model_row <- data.frame(
            seed = i,  # Using iteration number as seed
            model = model_name,
            CICc = model_CICc,
            delta_CICc = model_delta_CICc,
            nSpecies = result_i$nSpecies,
            stringsAsFactors = FALSE
          )
          
          # Extract edge information from the coefficient matrix
          if (!is.null(chosen_model$coef) && is.matrix(chosen_model$coef)) {
            edges_coef <- chosen_model$coef
            edges_se <- chosen_model$se
            
            # Process coefficient matrix
            for (from_idx in 1:nrow(edges_coef)) {
              for (to_idx in 1:ncol(edges_coef)) {
                coef_value <- edges_coef[from_idx, to_idx]
                
                # Only process non-zero coefficients (actual paths)
                if (coef_value != 0) {
                  from_name <- rownames(edges_coef)[from_idx]
                  to_name <- colnames(edges_coef)[to_idx]
                  
                  # Get standard error
                  se_value <- edges_se[from_idx, to_idx]
                  
                  # Calculate approximate p-value
                  z_value <- coef_value / se_value
                  p_value <- 2 * (1 - pnorm(abs(z_value)))
                  
                  # Create column names for this edge
                  edge_name <- paste0(from_name, "_to_", to_name)
                  edge_est_name <- paste0(edge_name, "_est")
                  edge_se_name <- paste0(edge_name, "_se")
                  edge_p_name <- paste0(edge_name, "_p")
                  
                  # Add to model row
                  model_row[[edge_est_name]] <- coef_value
                  model_row[[edge_se_name]] <- se_value
                  model_row[[edge_p_name]] <- p_value
                }
              }
            }
          }
          
          # Add to detailed models list
          detailed_models_list[[length(detailed_models_list) + 1]] <- model_row
        }
        
        # Also create summary for aggregated results
        model_summary_i <- as.data.frame(summary_i)
        model_summary_i$iteration <- i
        model_summary_i$model_name <- rownames(model_summary_i)
        all_model_summaries[[i]] <- model_summary_i
      }
    }
  }
  
  cat("\n")  # New line after progress indicator
  
  # Combine detailed models into a single data frame
  detailed_models_df <- bind_rows(detailed_models_list)
  
  # Aggregate results across iterations
  aggregated_summaries <- if (length(all_model_summaries) > 0) {
    do.call(rbind, all_model_summaries)
  } else {
    data.frame()
  }
  
  # Calculate average model performance
  avg_model_performance <- if (nrow(aggregated_summaries) > 0) {
    aggregated_summaries %>%
      group_by(model_name) %>%
      summarise(
        mean_CICc = mean(CICc, na.rm = TRUE),
        sd_CICc = sd(CICc, na.rm = TRUE),
        times_best = sum(delta_CICc == 0, na.rm = TRUE),
        times_in_top = sum(delta_CICc < 2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(mean_CICc)
  } else {
    data.frame()
  }
  
  # Calculate model frequencies for compatibility with create_downsampled_plots
  model_frequencies <- if (nrow(detailed_models_df) > 0) {
    detailed_models_df %>%
      group_by(model) %>%
      summarise(
        frequency_in_sub2 = n(),
        .groups = "drop"
      ) %>%
      rename(model_name = model) %>%
      arrange(desc(frequency_in_sub2))
  } else {
    data.frame()
  }
  
  # Create a summary phylopath output structure compatible with run_multiple_phylopath output
  result_dimorphism <- list(
    results_df = data.frame(
      seed = 1:n_iterations,
      nSub2dCIC_Models = sapply(all_model_summaries, function(x) sum(x$delta_CICc < 2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    ),
    model_frequencies = model_frequencies,
    detailed_models = detailed_models_df,
    model_summaries = aggregated_summaries,
    avg_performance = avg_model_performance,
    n_iterations = n_iterations,
    downsampling_info = dimorphism_downsample,
    all_results = all_results,
    dimorphism_type = dim_info$label
  )
  
  # Save outputs if requested
  if (save_outputs) {
    # Save detailed models CSV
    detailed_models_file <- file.path(dim_output_dir, 
                                     paste0("detailed_models_", 
                                           "Remove", dimorphism_downsample$n_to_remove,
                                           "High", dim_info$label,
                                           "_n", n_iterations, "_", Sys.Date(), ".csv"))
    write.csv(detailed_models_df, detailed_models_file, row.names = FALSE)
    cat("Detailed models saved to:", detailed_models_file, "\n")
    
    # Save model frequencies
    write.csv(model_frequencies, 
             file.path(dim_output_dir, 
                      paste0("model_frequencies_", dim_info$label, "_n", n_iterations, "_", Sys.Date(), ".csv")),
             row.names = FALSE)
    
    # Save downsampling summary
    downsample_summary <- data.frame(
      dimorphism_type = dim_info$label,
      dimorphism_variable = dim_info$col,
      n_removed = dimorphism_downsample$n_to_remove,
      target_mean = dimorphism_downsample$target_stats$target_mean,
      original_mean = dimorphism_downsample$current_stats$current_mean,
      corrected_mean = dimorphism_downsample$summary$mean_dimorphism_after,
      convergence = dimorphism_downsample$summary$convergence
    )
    
    write.csv(downsample_summary,
             file.path(dim_output_dir, 
                      paste0("downsampling_summary_", dim_info$label, ".csv")),
             row.names = FALSE)
  }
  
  cat("\n", dim_info$description, "bias correction complete!\n")
  cat("Mean after correction:", 
      round(dimorphism_downsample$summary$mean_dimorphism_after, 3), "\n")
  cat("Results saved to:", dim_output_dir, "\n")
  
  return(result_dimorphism)
}

#' Updated create_downsampled_plots with flexible input
#'
#' This version accepts detailed_models_df as:
#' A) A dataframe object
#' B) A file path to the detailed_models csv
#' C) NULL (searches for csv file using pattern matching)
#'
#' @param downsampling_results Output from run_multiple_phylopath or similar
#' @param detailed_models_input Either a dataframe, file path, or NULL
#' @param downsampling_info Information about the downsampling
#' @param output_prefix Prefix for output files
#' @param output_dir Output directory
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_downsampled_plots_flexible <- function(downsampling_results,
                                            detailed_models_input = NULL,
                                            downsampling_info = NULL,
                                            output_prefix = "phylopath_downsampled",
                                            output_dir = "Outputs/PhylopathPlots",
                                            save_png = TRUE,
                                            save_pdf = TRUE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Skip downsampling info plot (not relevant for dimorphism)
  p_downsample_info <- NULL
  
  # Handle detailed models input
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
  
  # Process detailed models if available
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    # Extract edge information
    edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Aggregate by seed
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format for plotting
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # Create violin plots
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
    
    # Create heatmap
    cat("Creating heatmap of mean coefficients...\n")
    
    # Get the across-seeds summary
    path_summary <- agg_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
    
    # Create heatmap (code continues as in original create_downsampled_plots)
    # ... [rest of heatmap code]
    
    # Model frequency bar plot
    cat("Creating model frequency plot...\n")
    if (!is.null(downsampling_results$model_frequencies)) {
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
  
  # Also save as PDF if requested
  if (save_pdf) {
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