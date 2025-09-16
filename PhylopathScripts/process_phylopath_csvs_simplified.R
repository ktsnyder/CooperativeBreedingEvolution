# process_phylopath_csvs_simplified.R
# Simplified runner script to process CSV files and generate PDF plots
# Kate Snyder
# Created: 2025-06-17

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)

# Source the aggregation functions
source("claude_code_sessions/PhylopathScripts/run_phylopath_fxns.R")

# Set the base directory
base_dir <- "/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathDownsampled"

# Find all detailed_models CSV files
detailed_models_files <- list.files(base_dir, 
                                   pattern = "detailed_models_.*\\.csv$", 
                                   recursive = TRUE, 
                                   full.names = TRUE)

cat("Found", length(detailed_models_files), "detailed models files to process\n\n")

# Process each file
for (file_path in detailed_models_files) {
  cat("\n=====================================\n")
  cat("Processing:", basename(file_path), "\n")
  
  # Extract info from filename
  filename <- basename(file_path)
  dir_name <- basename(dirname(file_path))
  
  # Extract the Remove pattern
  remove_pattern <- str_extract(filename, "Remove[0-9]+[A-Za-z0-9]+")
  
  # Extract number of iterations
  n_iterations <- str_extract(filename, "_([0-9]+)_[0-9]{4}-[0-9]{2}-[0-9]{2}", group = 1)
  if (is.na(n_iterations)) {
    n_iterations <- str_extract(filename, "_n([0-9]+)_", group = 1)
  }
  
  # Set output directory
  if (dir_name %in% c("PlumageDimorphism", "WingDimorphism")) {
    output_dir <- file.path(dirname(file_path), "Plots")
  } else {
    output_dir <- file.path(base_dir, "Plots", remove_pattern)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create output prefix
  output_prefix <- paste0("phylopath_", remove_pattern)
  if (!is.na(n_iterations)) {
    output_prefix <- paste0(output_prefix, "_", n_iterations, "iterations")
  }
  
  cat("Output directory:", output_dir, "\n")
  cat("Output prefix:", output_prefix, "\n")
  
  # Load the data
  detailed_models_df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check for model frequencies file
  freq_filename <- gsub("detailed_models_", "model_frequencies_", filename)
  freq_file <- file.path(dirname(file_path), freq_filename)
  
  model_frequencies <- NULL
  if (file.exists(freq_file)) {
    cat("Found model frequencies file\n")
    model_frequencies <- read.csv(freq_file, stringsAsFactors = FALSE)
  }
  
  # Process the data
  tryCatch({
    # Aggregate by seed
    cat("Aggregating results by seed...\n")
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # Create violin plot
    if (!is.null(seed_level_conditional) && nrow(seed_level_conditional) > 0) {
      cat("Creating violin plot...\n")
      
      # Clean path names
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
        
        # Save violin plot
        ggsave(file.path(output_dir, paste0(output_prefix, "_violin.pdf")),
               p_violin, width = 10, height = 8, device = "pdf")
        ggsave(file.path(output_dir, paste0(output_prefix, "_violin.png")),
               p_violin, width = 10, height = 8, dpi = 300)
        cat("  Saved violin plot\n")
      }
    }
    
    # Create model frequency plot
    if (!is.null(model_frequencies) && nrow(model_frequencies) > 0) {
      cat("Creating model frequency plot...\n")
      
      model_freq <- model_frequencies %>%
        arrange(desc(frequency_in_sub2)) %>%
        head(15)
      
      # Calculate proportion if not already present
      if (!"proportion_in_sub2" %in% names(model_freq)) {
        model_freq$proportion_in_sub2 <- model_freq$frequency_in_sub2 / 
                                         length(unique(detailed_models_df$seed))
      }
      
      p_model_freq <- ggplot(model_freq, 
                            aes(x = reorder(model_name, frequency_in_sub2),
                                y = frequency_in_sub2)) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                    round(proportion_in_sub2 * 100, 1), "%)")),
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
      
      # Save model frequency plot
      ggsave(file.path(output_dir, paste0(output_prefix, "_model_freq.pdf")),
             p_model_freq, width = 10, height = 8, device = "pdf")
      ggsave(file.path(output_dir, paste0(output_prefix, "_model_freq.png")),
             p_model_freq, width = 10, height = 8, dpi = 300)
      cat("  Saved model frequency plot\n")
    }
    
    # Create heatmap
    path_summary <- agg_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
    
    if (!is.null(path_summary) && nrow(path_summary) > 0) {
      cat("Creating heatmap...\n")
      
      # Create from-to matrix
      unique_from <- unique(path_summary$from)
      unique_to <- unique(path_summary$to)
      
      coef_matrix <- matrix(NA, 
                           nrow = length(unique_from), 
                           ncol = length(unique_to),
                           dimnames = list(unique_from, unique_to))
      
      for (i in 1:nrow(path_summary)) {
        coef_matrix[path_summary$from[i], path_summary$to[i]] <- path_summary$mean_coef_conditional[i]
      }
      
      # Convert to long format
      coef_long <- as.data.frame(as.table(coef_matrix))
      names(coef_long) <- c("From", "To", "Coefficient")
      coef_long <- coef_long[!is.na(coef_long$Coefficient), ]
      
      # Clean variable names
      clean_names <- function(x) {
        x <- gsub("FemaleSong_Agg01", "Female\\nSong", x)
        x <- gsub("HighConfidence_Coop", "Cooperation", x)
        x <- gsub("TerritorialityWeakVsStrong", "Territoriality", x)
        x <- gsub("Territory_12vs3", "Territory\\nType", x)
        x <- gsub("logMass_AVONET", "Body Mass", x)
        x <- gsub("logMaleFemalePlumageDiffAbs", "Plumage\\nDimorphism", x)
        x <- gsub("PercentAbsLogWingDimorphism", "Wing\\nDimorphism", x)
        return(x)
      }
      
      coef_long$From <- clean_names(coef_long$From)
      coef_long$To <- clean_names(coef_long$To)
      
      # Add significance info
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
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = paste0(round(Coefficient, 3), 
                                    ifelse(is.na(sig_label), "", sig_label))),
                 color = ifelse(abs(coef_long$Coefficient) > 0.3, "white", "black"),
                 size = 4) +
        scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2E86AB",
                           midpoint = 0, name = "Mean\\nCoefficient",
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
      
      # Save heatmap
      ggsave(file.path(output_dir, paste0(output_prefix, "_heatmap.pdf")),
             p_heatmap, width = 8, height = 8, device = "pdf")
      ggsave(file.path(output_dir, paste0(output_prefix, "_heatmap.png")),
             p_heatmap, width = 8, height = 8, dpi = 300)
      cat("  Saved heatmap\n")
    }
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  })
}

cat("\n\nProcessing complete!\n")
cat("PDFs have been saved to the respective output directories.\n")