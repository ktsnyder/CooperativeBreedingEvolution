# process_phylopath_complete.R
# Complete runner script to process CSV files and generate all phylopath PDF plots
# Includes violin plots with full dataset comparison (red lines)
# Kate Snyder
# Created: 2025-06-17

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
library(phylopath)
library(ape)

# Source the aggregation functions
source("claude_code_sessions/PhylopathScripts/run_phylopath_fxns.R")

# Set the base directory
base_dir <- "/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathDownsampled"

# Load the full dataset and tree for comparison
cat("Loading full dataset and tree...\n")
dfIn_phylo <- read.csv("Data_R_2025-06-09.csv")

# Find tree file - prefer OscineSubset
tree_files <- list.files(".", pattern = "ConsensusPasserineTree.*OscineSubset\\.nex", recursive = TRUE, full.names = TRUE)
if (length(tree_files) == 0) {
  tree_files <- list.files(".", pattern = "ConsensusPasserineTree.*\\.nex", recursive = TRUE, full.names = TRUE)
}
if (length(tree_files) == 0) {
  stop("No consensus tree file found!")
}
tree_file <- tree_files[1]
cat("Using tree file:", tree_file, "\n")
tree_list <- read.nexus(tree_file)
if (class(tree_list) == "multiPhylo") {
  tree_phylo <- tree_list[[1]]
} else {
  tree_phylo <- tree_list
}
cat("Tree has", length(tree_phylo$tip.label), "tips\n")

# Find all detailed_models CSV files
detailed_models_files <- list.files(base_dir, 
                                   pattern = "detailed_models_.*\\.csv$", 
                                   recursive = TRUE, 
                                   full.names = TRUE)

cat("Found", length(detailed_models_files), "detailed models files to process\n\n")

# Function to run phylopath on full dataset
run_full_dataset_phylopath <- function(dfIn, tree, territoriality_var = "TerritorialityWeakVsStrong", 
                                      mass_var = "logMass_AVONET") {
  df_filtered <- dfIn[dfIn$species %in% tree$tip.label, ]
  rownames(df_filtered) <- df_filtered$species
  
  result <- run_CB_FS_Terr_phylopath(
    dfIn = df_filtered,
    tree = tree,
    female_song_var = "FemaleSong_Agg01",
    coop_breeding_var = "HighConfidence_Coop",
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    plots2pdf = FALSE
  )
  
  return(result)
}

# Cache for full dataset results
full_dataset_cache <- list()

# Process each file
for (file_path in detailed_models_files) {
  cat("\n=====================================\n")
  cat("Processing:", basename(file_path), "\n")
  
  # Extract info
  filename <- basename(file_path)
  dir_name <- basename(dirname(file_path))
  remove_pattern <- str_extract(filename, "Remove[0-9]+[A-Za-z0-9]+")
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
  
  output_prefix <- paste0("phylopath_", remove_pattern)
  if (!is.na(n_iterations)) {
    output_prefix <- paste0(output_prefix, "_", n_iterations, "iterations")
  }
  
  # Determine variables
  if (grepl("PlumageDimorphism", file_path)) {
    mass_var <- "logMaleFemalePlumageDiffAbs"
  } else if (grepl("WingDimorphism", file_path)) {
    mass_var <- "PercentAbsLogWingDimorphism"
  } else {
    mass_var <- "logMass_AVONET"
  }
  
  if (grepl("Territory_12vs3|Terr3", filename)) {
    territoriality_var <- "Territory_12vs3"
  } else {
    territoriality_var <- "TerritorialityWeakVsStrong"
  }
  
  # Get full dataset results
  cache_key <- paste(territoriality_var, mass_var, sep = "_")
  if (cache_key %in% names(full_dataset_cache)) {
    cat("Using cached full dataset results\n")
    full_result <- full_dataset_cache[[cache_key]]
  } else {
    cat("Running phylopath on full dataset\n")
    full_result <- run_full_dataset_phylopath(dfIn_phylo, tree_phylo, territoriality_var, mass_var)
    full_dataset_cache[[cache_key]] <- full_result
  }
  
  # Extract full dataset coefficients
  full_data_coefficients <- NULL
  if (!is.null(full_result) && !is.null(full_result$result)) {
    full_summary <- summary(full_result$result)
    best_models <- full_summary[full_summary$delta_CICc < 2, ]
    
    if (nrow(best_models) > 0) {
      all_paths <- list()
      
      for (i in 1:nrow(best_models)) {
        model_name <- best_models$model[i]
        chosen_model <- choice(full_result$result, model_name)
        
        if (!is.null(chosen_model$coef) && is.matrix(chosen_model$coef)) {
          coef_matrix <- chosen_model$coef
          
          for (from_idx in 1:nrow(coef_matrix)) {
            for (to_idx in 1:ncol(coef_matrix)) {
              coef_value <- coef_matrix[from_idx, to_idx]
              if (coef_value != 0) {
                from_name <- rownames(coef_matrix)[from_idx]
                to_name <- colnames(coef_matrix)[to_idx]
                path_name <- paste0(from_name, "_to_", to_name)
                all_paths[[path_name]] <- c(all_paths[[path_name]], coef_value)
              }
            }
          }
        }
      }
      
      if (length(all_paths) > 0) {
        full_data_coefficients <- data.frame(
          path = names(all_paths),
          full_data_mean = sapply(all_paths, mean),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Load data
  detailed_models_df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check for model frequencies
  freq_filename <- gsub("detailed_models_", "model_frequencies_", filename)
  freq_file <- file.path(dirname(file_path), freq_filename)
  model_frequencies <- NULL
  if (file.exists(freq_file)) {
    model_frequencies <- read.csv(freq_file, stringsAsFactors = FALSE)
  }
  
  # Create plots using the complete function from run_phylopath_fxns.R
  tryCatch({
    cat("Creating plots...\n")
    
    # Create downsampling results structure
    downsampling_results <- list(
      results_df = data.frame(seed = unique(detailed_models_df$seed)),
      model_frequencies = model_frequencies
    )
    
    # Call the complete plotting function
    plots <- create_downsampled_plots(
      downsampling_results = downsampling_results,
      detailed_models_input = detailed_models_df,
      model_frequencies_input = model_frequencies,
      full_data_phylopath_input = full_result,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = TRUE,
      save_pdf = TRUE
    )
    
    cat("Successfully created all plots!\n")
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  })
}

cat("\n\nProcessing complete!\n")
cat("All PDFs (violin plots with red lines, heatmaps, and model frequencies) have been created.\n")
cat("\nOutput locations:\n")
cat("- Main analyses:", file.path(base_dir, "Plots"), "\n")
cat("- Dimorphism analyses:", file.path(base_dir, "[PlumageDimorphism|WingDimorphism]", "Plots"), "\n")