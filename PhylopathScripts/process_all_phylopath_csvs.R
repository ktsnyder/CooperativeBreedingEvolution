# process_all_phylopath_csvs.R
# Runner script to process all CSV files in PhylopathDownsampled and generate PDF plots
# Kate Snyder
# Created: 2025-06-17

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Source the plotting functions
source("claude_code_sessions/PhylopathScripts/run_phylopath_fxns.R")
source("claude_code_sessions/PhylopathScripts/phylopath_plotting_comprehensive_fixed.R")
source("claude_code_sessions/PhylopathScripts/phylopath_dimorphism_improvements.R")

# Set the base directory
base_dir <- "/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathDownsampled"

# Find all CSV files
csv_files <- list.files(base_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

cat("Found", length(csv_files), "CSV files to process\n\n")

# Group files by analysis type based on file patterns
detailed_models_files <- grep("detailed_models_", csv_files, value = TRUE)
model_frequencies_files <- grep("model_frequencies_", csv_files, value = TRUE)
edge_summary_files <- grep("edge_summary_", csv_files, value = TRUE)
iteration_means_files <- grep("iteration_means\\.csv$", csv_files, value = TRUE)
summary_statistics_files <- grep("summary_statistics\\.csv$", csv_files, value = TRUE)

cat("File types found:\n")
cat("- Detailed models:", length(detailed_models_files), "\n")
cat("- Model frequencies:", length(model_frequencies_files), "\n")
cat("- Edge summaries:", length(edge_summary_files), "\n")
cat("- Iteration means:", length(iteration_means_files), "\n")
cat("- Summary statistics:", length(summary_statistics_files), "\n\n")

# Function to extract analysis info from filename
extract_analysis_info <- function(filepath) {
  filename <- basename(filepath)
  dir_name <- basename(dirname(filepath))
  
  # Extract the Remove pattern (e.g., "Remove155Terr3", "Remove167HighPlumageDimorphism")
  remove_pattern <- str_extract(filename, "Remove[0-9]+[A-Za-z0-9]+")
  
  # Extract number of iterations
  n_iterations <- str_extract(filename, "_([0-9]+)_[0-9]{4}-[0-9]{2}-[0-9]{2}", group = 1)
  if (is.na(n_iterations)) {
    n_iterations <- str_extract(filename, "_n([0-9]+)_", group = 1)
  }
  
  # Determine analysis type
  if (grepl("PlumageDimorphism|WingDimorphism", filepath)) {
    analysis_type <- dir_name
  } else if (grepl("Terr", remove_pattern)) {
    analysis_type <- "Territoriality"
  } else if (grepl("Coop", remove_pattern)) {
    analysis_type <- "Cooperation"
  } else {
    analysis_type <- "Other"
  }
  
  return(list(
    remove_pattern = remove_pattern,
    n_iterations = n_iterations,
    analysis_type = analysis_type,
    filename = filename,
    directory = dirname(filepath)
  ))
}

# Process detailed models files
cat("Processing detailed models files...\n")

for (i in seq_along(detailed_models_files)) {
  file_path <- detailed_models_files[i]
  info <- extract_analysis_info(file_path)
  
  cat("\nProcessing:", info$filename, "\n")
  cat("  Analysis type:", info$analysis_type, "\n")
  cat("  Remove pattern:", info$remove_pattern, "\n")
  cat("  Iterations:", info$n_iterations, "\n")
  
  # Create output prefix based on the file information
  output_prefix <- paste0("phylopath ", info$remove_pattern)
  if (!is.na(info$n_iterations)) {
    output_prefix <- paste0(output_prefix, " ", info$n_iterations, "iterations")
  }
  
  # Set output directory
  if (info$analysis_type %in% c("PlumageDimorphism", "WingDimorphism")) {
    output_dir <- file.path(info$directory, "Plots")
  } else {
    output_dir <- file.path(base_dir, "Plots", info$analysis_type)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  Created output directory:", output_dir, "\n")
  }
  
  # Load the detailed models data
  detailed_models_df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check if we have corresponding model frequencies file
  freq_pattern <- gsub("detailed_models_", "model_frequencies_", info$filename)
  freq_file <- file.path(info$directory, freq_pattern)
  
  model_frequencies <- NULL
  if (file.exists(freq_file)) {
    cat("  Found corresponding model frequencies file\n")
    model_frequencies <- read.csv(freq_file, stringsAsFactors = FALSE)
  }
  
  # Create a minimal downsampling results structure
  downsampling_results <- list(
    results_df = data.frame(
      seed = unique(detailed_models_df$seed),
      stringsAsFactors = FALSE
    ),
    model_frequencies = model_frequencies,
    detailed_models = detailed_models_df
  )
  
  # Try to create plots - use the comprehensive fixed function which is more complete
  tryCatch({
    cat("  Creating plots...\n")
    
    plots <- create_downsampled_plots(
      downsampling_results = downsampling_results,
      downsampling_info = NULL,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = TRUE,
      save_pdf = TRUE
    )
    
    cat("  Successfully created plots!\n")
    
  }, error = function(e) {
    cat("  ERROR creating plots:", e$message, "\n")
    
    # Try creating individual plots if the full function fails
    cat("  Trying to create individual plots...\n")
    
    tryCatch({
      # At minimum, try to create the model frequency plot
      if (!is.null(model_frequencies)) {
        model_freq <- model_frequencies %>%
          arrange(desc(frequency_in_sub2)) %>%
          head(15) %>%
          mutate(
            prop_sub2 = frequency_in_sub2 / max(frequency_in_sub2)
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
            subtitle = paste("Based on", info$n_iterations, "iterations"),
            x = NULL,
            y = "Frequency in Top Model Set"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 9)
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
        
        # Save the plot
        ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.png")),
               p_model_freq, width = 10, height = 8, dpi = 300)
        ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.pdf")),
               p_model_freq, width = 10, height = 8, device = "pdf")
        
        cat("    Created model frequency plot\n")
      }
      
    }, error = function(e2) {
      cat("    ERROR creating individual plots:", e2$message, "\n")
    })
  })
}

# Process dimorphism-specific files that might not have detailed_models prefix
cat("\n\nProcessing dimorphism-specific analyses...\n")

dimorphism_dirs <- c("PlumageDimorphism", "WingDimorphism")

for (dim_type in dimorphism_dirs) {
  dim_dir <- file.path(base_dir, dim_type)
  
  if (dir.exists(dim_dir)) {
    cat("\nProcessing", dim_type, "files...\n")
    
    # Look for iteration means files
    iter_means_files <- list.files(dim_dir, pattern = "iteration_means\\.csv$", full.names = TRUE)
    
    for (file_path in iter_means_files) {
      info <- extract_analysis_info(file_path)
      cat("  Found iteration means:", basename(file_path), "\n")
      
      # Create distribution plot
      tryCatch({
        iter_means_df <- read.csv(file_path, stringsAsFactors = FALSE)
        
        # Extract the remove pattern for output naming
        output_prefix <- info$remove_pattern
        output_dir <- file.path(dim_dir, "Plots")
        
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        # Create a simple distribution plot
        p <- ggplot(iter_means_df, aes(x = mean_dimorphism)) +
          geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
          geom_vline(aes(xintercept = mean(mean_dimorphism)), 
                     color = "red", linetype = "dashed", size = 1) +
          labs(
            title = paste(dim_type, "Distribution Across Iterations"),
            subtitle = paste(output_prefix, "-", length(iter_means_df$iteration), "iterations"),
            x = paste(dim_type, "Mean"),
            y = "Count"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12)
          )
        
        # Save as PNG and PDF
        ggsave(
          filename = file.path(output_dir, paste0(output_prefix, "_distribution.png")),
          plot = p,
          width = 8,
          height = 6,
          dpi = 300
        )
        
        ggsave(
          filename = file.path(output_dir, paste0(output_prefix, "_distribution.pdf")),
          plot = p,
          width = 8,
          height = 6,
          device = "pdf"
        )
        
        cat("    Created distribution plot\n")
        
      }, error = function(e) {
        cat("    ERROR creating distribution plot:", e$message, "\n")
      })
    }
  }
}

cat("\n\nProcessing complete!\n")
cat("Check the following directories for output PDFs:\n")
cat("- ", file.path(base_dir, "Plots"), "\n")
cat("- ", file.path(base_dir, "PlumageDimorphism", "Plots"), "\n")
cat("- ", file.path(base_dir, "WingDimorphism", "Plots"), "\n")