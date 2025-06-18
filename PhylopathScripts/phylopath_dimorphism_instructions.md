# Instructions for Implementing Phylopath Dimorphism Improvements

## Overview

This document provides instructions for:
1. Replacing the dimorphism bias correction code in Run_Analyses.R with a version that creates `detailed_models_df`
2. Updating `create_downsampled_plots()` and `create_all_phylopath_plots()` to accept alternative input formats

## 1. Replace Dimorphism Code in Run_Analyses.R

Replace the entire loop in the dimorphism bias correction section (lines 1099-1288) with:

```r
# Source the improvements file
source("claude_code_sessions/phylopath_dimorphism_improvements.R")

# Loop through each dimorphism variable
for (dim_type in names(dimorphism_vars)) {
  dim_info <- dimorphism_vars[[dim_type]]
  
  # Run the improved phylopath dimorphism correction
  result_dimorphism <- run_phylopath_dimorphism_correction(
    dfIn_phylo = dfIn_phylo,
    tree = tree_phylo,
    dim_info = dim_info,
    n_iterations = n_iterations,
    phylopath_output_dir = phylopath_output_dir,
    save_outputs = TRUE
  )
  
  # If results were obtained, create plots
  if (!is.null(result_dimorphism)) {
    # Create output prefix with dimorphism type
    prefix_dimorphism <- paste0("Remove", 
                               result_dimorphism$downsampling_info$n_to_remove, 
                               "High", 
                               dim_info$label)
    
    # Use the flexible plotting function
    plots_dimorphism <- create_downsampled_plots_flexible(
      downsampling_results = result_dimorphism,
      detailed_models_input = result_dimorphism$detailed_models,  # Pass dataframe directly
      downsampling_info = NULL,
      output_prefix = prefix_dimorphism,
      output_dir = file.path(phylopath_output_dir, dim_info$label),
      save_png = TRUE
    )
  }
}

cat("\n\nAll dimorphism bias corrections complete!\n")
```

## 2. Update create_downsampled_plots() in run_phylopath_fxns.R

Add this parameter to the function signature (around line 1940):

```r
create_downsampled_plots <- function(downsampling_results,
                                   detailed_models_input = NULL,  # NEW PARAMETER
                                   downsampling_info = NULL,
                                   output_prefix = "phylopath_downsampled",
                                   output_dir = "Outputs/PhylopathPlots",
                                   save_png = TRUE,
                                   save_pdf = TRUE) {
```

Then replace the section that loads detailed models (lines 1994-2018) with:

```r
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
    # Case C: Search for file using pattern (existing code)
    # Fix: Look for files with pattern that matches the output_prefix
    detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", 
                                                        gsub(" ", "_", output_prefix), 
                                                        ".*\\.csv$"))
    
    # If no files found, try a more general pattern
    if (length(detailed_models_files) == 0) {
      # Extract the key part of the prefix (e.g., "Remove66HighTerr")
      prefix_parts <- strsplit(output_prefix, " ")[[1]]
      key_pattern <- grep("Remove", prefix_parts, value = TRUE)
      if (length(key_pattern) > 0) {
        detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", key_pattern, ".*\\.csv$"))
      }
    }
    
    detailed_models_file <- if (length(detailed_models_files) > 0) detailed_models_files[1] else NA
    
    if (!is.na(detailed_models_file)) {
      cat("Loading detailed models from:", detailed_models_file, "\n")
      detailed_models_df <- read.csv(detailed_models_file)
    }
    
  } else if (!is.null(downsampling_results$detailed_models)) {
    # Check if detailed_models is included in the results object
    cat("Using detailed models from results object\n")
    detailed_models_df <- downsampling_results$detailed_models
  }
  
  # Continue with rest of function using detailed_models_df...
```

## 3. Update create_all_phylopath_plots() wrapper

In `create_all_phylopath_plots()` (around line 2205), add the new parameter and pass it through:

```r
create_all_phylopath_plots <- function(analysis_type = c("nondownsampled", "downsampled"),
                                      phylopath_output,
                                      detailed_models_input = NULL,  # NEW PARAMETER
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
      detailed_models_input = detailed_models_input,  # PASS THROUGH
      downsampling_info = downsampling_info,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = save_png,
      save_pdf = save_pdf
    )
  }
  
  return(plots)
}
```

## 4. Alternative: Use the Flexible Function Directly

Instead of modifying the existing functions, you can use `create_downsampled_plots_flexible()` from the improvements file directly:

```r
# After running dimorphism correction
plots <- create_downsampled_plots_flexible(
  downsampling_results = result_dimorphism,
  detailed_models_input = result_dimorphism$detailed_models,  # or a file path, or NULL
  output_prefix = "Remove50HighPlumageDimorphism",
  output_dir = "Outputs/PhylopathDownsampled/PlumageDimorphism",
  save_png = TRUE
)
```

## Key Improvements

1. **`run_phylopath_dimorphism_correction()`** properly creates a `detailed_models_df` in the same format as `run_multiple_phylopath()`
2. **Flexible input handling** allows passing detailed models as:
   - A dataframe object (fastest, no file I/O)
   - A file path (explicit control)
   - NULL (automatic search with improved pattern matching)
3. **Consistent output structure** ensures compatibility with existing plotting functions
4. **Better error handling** for cases where no models meet the criteria

## Testing

To test the changes:

```r
# Test with dimorphism data
dim_info <- list(
  col = "logMaleFemalePlumageDiffAbs",
  label = "PlumageDimorphism",
  description = "log Plumage Dimorphism (Absolute Value)"
)

result <- run_phylopath_dimorphism_correction(
  dfIn_phylo = dfIn_phylo,
  tree = tree_phylo,
  dim_info = dim_info,
  n_iterations = 5,  # Small number for testing
  save_outputs = TRUE
)

# Check that detailed_models_df was created
str(result$detailed_models)
```