# Extract summary statistics from all Simmap Overlap analyses
# Kate Snyder / Claude
# 2025-07-01

library(tidyverse)

# Function to extract info from folder name
parse_folder_name <- function(folder_name) {
  # Pattern: trait1_vs_trait2_nsims_real_nsims_dummy_label_timestamp
  parts <- strsplit(folder_name, "_vs_")[[1]]
  trait1 <- parts[1]
  
  # Split the rest
  remaining <- parts[2]
  parts2 <- strsplit(remaining, "_")[[1]]
  
  # Find where the numbers start (nsims_real)
  # Look for the first numeric value
  numeric_start <- which(grepl("^[0-9]+$", parts2))[1]
  
  if (is.na(numeric_start)) {
    return(list(trait1 = trait1, trait2 = NA, nsims_real = NA, nsims_dummy = NA))
  }
  
  # trait2 is everything before the first number
  trait2 <- paste(parts2[1:(numeric_start-1)], collapse = "_")
  
  # Extract nsims
  nsims_real <- as.numeric(parts2[numeric_start])
  nsims_dummy <- as.numeric(parts2[numeric_start + 1])
  
  return(list(
    trait1 = trait1,
    trait2 = trait2,
    nsims_real = nsims_real,
    nsims_dummy = nsims_dummy
  ))
}

# Set directory
base_dir <- "Simmap_Overlap_Outputs"

# Get all subdirectories
all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Initialize results dataframe
results_list <- list()

# Process each directory
for (i in seq_along(all_dirs)) {
  dir_path <- all_dirs[i]
  folder_name <- basename(dir_path)
  
  cat("Processing:", folder_name, "\n")
  
  # Look for summary CSV file
  summary_files <- list.files(file.path(dir_path, "summaries"), 
                             pattern = "analysis_summary.*\\.csv$", 
                             full.names = TRUE)
  
  if (length(summary_files) == 0) {
    cat("  No summary file found, skipping\n")
    next
  }
  
  # Read the summary file
  summary_file <- summary_files[1]  # Take first if multiple
  summary_data <- read.csv(summary_file)
  
  # Parse folder name
  folder_info <- parse_folder_name(folder_name)
  
  # Extract required information
  result_row <- data.frame(
    FolderName = folder_name,
    trait1 = folder_info$trait1,
    trait2 = folder_info$trait2,
    nsims_real = folder_info$nsims_real,
    nsims_dummy = folder_info$nsims_dummy,
    stringsAsFactors = FALSE
  )
  
  # Add number of species (if available)
  if ("Nspecies" %in% names(summary_data)) {
    result_row$NumberOfSpecies <- summary_data$Nspecies[1]
  } else {
    result_row$NumberOfSpecies <- NA
  }
  
  # Add empirical p-value
  if ("p_value" %in% names(summary_data)) {
    result_row$Empirical_pval <- summary_data$p_value[1]
  } else {
    result_row$Empirical_pval <- NA
  }
  
  # Add fractions dummy <= median real
  fraction_cols <- grep("FractionDummyLessThanMedianReal", names(summary_data), value = TRUE)
  
  if (length(fraction_cols) > 0) {
    # Add each fraction column
    for (col in fraction_cols) {
      result_row[[col]] <- summary_data[[col]][1]
    }
  } else {
    # If not found, try alternative naming
    # Check for individual columns
    for (state in c("0_0", "0_1", "1_0", "1_1")) {
      col_name <- paste0("FractionDummyLessThanMedianReal_", state)
      alt_col_name <- paste0("Overlap_", state, "_FractionDummyLTEMedianReal")
      
      if (alt_col_name %in% names(summary_data)) {
        result_row[[col_name]] <- summary_data[[alt_col_name]][1]
      } else {
        result_row[[col_name]] <- NA
      }
    }
  }
  
  # Add to results list
  results_list[[i]] <- result_row
}

# Combine all results
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Ensure all fraction columns exist
  expected_cols <- c("FractionDummyLessThanMedianReal_0_0",
                     "FractionDummyLessThanMedianReal_0_1", 
                     "FractionDummyLessThanMedianReal_1_0",
                     "FractionDummyLessThanMedianReal_1_1")
  
  for (col in expected_cols) {
    if (!(col %in% names(results_df))) {
      results_df[[col]] <- NA
    }
  }
  
  # Reorder columns
  col_order <- c("FolderName", "trait1", "trait2", "nsims_real", "nsims_dummy", 
                 "NumberOfSpecies", "Empirical_pval", expected_cols)
  
  results_df <- results_df[, col_order]
  
  # Sort by traits for easier reading
  results_df <- results_df %>%
    arrange(trait1, trait2, desc(nsims_real))
  
  # Save to CSV
  output_file <- paste0("simmap_overlap_results_summary_", format(Sys.Date(), "%Y%m%d"), ".csv")
  write.csv(results_df, output_file, row.names = FALSE)
  
  cat("\nResults saved to:", output_file, "\n")
  cat("Total analyses found:", nrow(results_df), "\n")
  
  # Display the table
  print(results_df)
  
} else {
  cat("No results found to compile.\n")
}

# Also create a simplified view with key statistics
if (exists("results_df") && nrow(results_df) > 0) {
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("Total analyses:", nrow(results_df), "\n")
  cat("Significant results (p < 0.05):", sum(results_df$Empirical_pval < 0.05, na.rm = TRUE), "\n")
  cat("Non-significant results:", sum(results_df$Empirical_pval >= 0.05, na.rm = TRUE), "\n")
  
  # Show fraction patterns
  cat("\n=== FRACTION PATTERNS ===\n")
  cat("(Values near 0 or 1 indicate strong deviation from independence)\n\n")
  
  # Calculate mean fractions for each state
  fraction_means <- results_df %>%
    summarise(
      Mean_0_0 = mean(FractionDummyLessThanMedianReal_0_0, na.rm = TRUE),
      Mean_0_1 = mean(FractionDummyLessThanMedianReal_0_1, na.rm = TRUE),
      Mean_1_0 = mean(FractionDummyLessThanMedianReal_1_0, na.rm = TRUE),
      Mean_1_1 = mean(FractionDummyLessThanMedianReal_1_1, na.rm = TRUE)
    )
  
  print(round(fraction_means, 3))
}