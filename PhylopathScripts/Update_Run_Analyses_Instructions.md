# Instructions for Updating Run_Analyses.R

## Changes to Make in Run_Analyses.R

### 1. Update the source statements (around lines 649-651)

Replace:
```r
source("run_phylopath_fxns.R")
source("phylopath_plotting_comprehensive_fixed.R")
source("create_phylopath_plots_enhanced_fixed.R")
```

With:
```r
source("claude_code_sessions/run_phylopath_fxns_consolidated.R")
source("create_phylopath_plots_enhanced_fixed.R")  # Keep this if still needed
```

### 2. Remove the redundant source statement (around line 801)

Remove:
```r
source("calculate_downsampling_function.R")  # probably obsolete because of calculate_stratified_downsampling()
```

This is no longer needed as the territoriality downsampling function is now included in the consolidated file.

### 3. Ensure calculate_stratified_downsampling_with_territoriality.R is still sourced (around line 735)

Keep:
```r
source("calculate_stratified_downsampling_with_territoriality.R")
```

This is still needed for the stratified downsampling calculations.

## What's Fixed

1. **Violin plots will now work** - The file pattern matching issue has been fixed, so the detailed_models CSV files will be found correctly.

2. **All functions are in one place** - No more scattered functions between multiple files.

3. **Territoriality downsampling is properly integrated** - The calculate_territoriality_downsampling function is included in the consolidated file.

## Testing

After making these changes, test by running a phylopath analysis with downsampling (e.g., the Holarctic non-cooperative downsampling). The violin plots should now show the actual data instead of being empty.

## Files That Can Be Removed (Optional)

Once you confirm everything is working, you may consider removing:
- `calculate_downsampling_function.R` (redundant with the territoriality downsampling in the consolidated file)

But keep:
- `calculate_stratified_downsampling_with_territoriality.R` (still needed for the main downsampling calculations)
- `create_phylopath_plots_enhanced_fixed.R` (if you're still using the enhanced plots functionality)