# Implementation Guide for create_downsampled_plots_complete

## Overview

`create_downsampled_plots_complete()` is a fully-featured replacement for `create_downsampled_plots()` that works with all types of downsampling outputs:
- Geographic bias correction (Holarctic, Tropical, Global cooperative)
- Territoriality bias correction (TerritorialityWeakVsStrong, Territory_12vs3)
- Dimorphism bias correction (Plumage, Wing)

## Key Features

1. **Flexible detailed_models input**:
   - Accepts a dataframe directly (no file I/O)
   - Accepts a file path
   - Automatically searches for files when given NULL

2. **Complete plotting functionality**:
   - Downsampling information plot (for geographic/territoriality)
   - Violin plots of path coefficients
   - Heatmap of mean coefficients with significance indicators
   - Model frequency bar plot

3. **Robust file finding**:
   - Searches in current directory
   - Falls back to searching in output directory
   - Uses flexible pattern matching

## Usage Examples

### 1. Geographic Bias Correction (e.g., Holarctic Non-cooperative)

```r
# After running run_multiple_phylopath()
source("claude_code_sessions/create_downsampled_plots_complete.R")

# Create plots with downsampling info
plots_geo <- create_downsampled_plots_complete(
  downsampling_results = result_geo_holarctic,
  detailed_models_input = NULL,  # Will auto-search for file
  downsampling_info = geo_holarctic_info,  # Include for info plot
  output_prefix = "Remove83HolarcticNoncoop",
  output_dir = "Outputs/PhylopathDownsampled",
  save_png = TRUE
)
```

### 2. Territoriality Bias Correction

```r
# After calculating territoriality downsampling
plots_terr <- create_downsampled_plots_complete(
  downsampling_results = result_terr,
  detailed_models_input = NULL,
  downsampling_info = terr_downsample_info$report,  # From calculate_downsampling()
  output_prefix = paste0("Remove", terr_downsample_info$n_to_remove, "StrongTerr"),
  output_dir = "Outputs/PhylopathDownsampled",
  save_png = TRUE
)
```

### 3. Dimorphism Bias Correction

```r
# Using output from run_phylopath_dimorphism_correction()
plots_dimorphism <- create_downsampled_plots_complete(
  downsampling_results = result_dimorphism,
  detailed_models_input = result_dimorphism$detailed_models,  # Pass dataframe directly
  downsampling_info = NULL,  # No info plot needed for dimorphism
  output_prefix = paste0("Remove", n_to_remove, "HighPlumageDimorphism"),
  output_dir = file.path("Outputs/PhylopathDownsampled", "PlumageDimorphism"),
  save_png = TRUE
)
```

## Replacing create_downsampled_plots in run_phylopath_fxns.R

To use this as a complete replacement, you have two options:

### Option 1: Replace the function entirely

In `run_phylopath_fxns.R`, replace the entire `create_downsampled_plots()` function with the contents of `create_downsampled_plots_complete()`, renaming it to `create_downsampled_plots`.

### Option 2: Add as an alias

Add this at the end of `run_phylopath_fxns.R`:

```r
# Source the complete version
source("claude_code_sessions/create_downsampled_plots_complete.R")

# Create an alias
create_downsampled_plots <- create_downsampled_plots_complete
```

## Updating create_all_phylopath_plots()

Update the wrapper function to pass through the new parameter:

```r
create_all_phylopath_plots <- function(analysis_type = c("nondownsampled", "downsampled"),
                                      phylopath_output,
                                      detailed_models_input = NULL,  # NEW
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
    plots <- create_downsampled_plots(  # Now uses the complete version
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

## Differences from Original

The complete version includes:

1. **Full heatmap implementation** with proper significance labeling
2. **Better handling of missing data** - won't crash if no paths meet criteria
3. **More variable name cleaning** - includes dimorphism variables
4. **Improved file searching** - checks output directory if not found in current directory
5. **Consistent theming** - uses theme_cowplot throughout

## Testing

Run the test function to verify everything works:

```r
source("claude_code_sessions/create_downsampled_plots_complete.R")
test_results <- test_create_downsampled_plots_complete()
```

This will test all three input methods and confirm the function handles them correctly.