# PhyloGLM Batch Analysis Error Summary

## Issue Description
When running `run_phyloglm_batch()` with 20 analyses and 1000 bootstrap iterations, all analyses fail with the error: "cannot open the connection"

## Error Details
- Error occurs during individual analysis execution, not during result integration
- The error happens after data preparation but before model fitting completes
- Reduced bootstrap iterations (10) work fine, suggesting a resource issue
- The error is not related to the new analysis names (tested and confirmed)

## Relevant Files

### Main Scripts
- `phyloglm_framework/batch_runner.R` - Main batch runner function
- `phyloglm_framework/batch_runner_helpers.R` - Helper functions (already fixed integration bugs)
- `phyloglm_framework/config_builder.R` - Creates 20 analysis configurations with new naming
- `phyloglm_framework/formula_builder.R` - Builds model formulas (updated to support fourway interactions)

### Data Files
- `Data_R_2025-06-09.csv` - Main data file
- `2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex` - Phylogenetic tree

## What Has Been Tried

### 1. Tested Individual Analysis
```r
# This works fine with bootstrap_n = 0
result <- run_single_analysis(config, data, tree, output_dir, bootstrap_n = 0, save_outputs = TRUE)
```

### 2. Tested Small Batch
```r
# This works with bootstrap_n = 10
batch_results <- run_phyloglm_batch(configs[1], data, tree, output_dir, bootstrap_n = 10)
```

### 3. Fixed Integration Bugs
- Added checks for empty data frames in `integrate_batch_results()`
- Added row count validation before adding Analysis column
- These fixes are in place but error occurs before integration

### 4. Attempted Solutions That Failed
- Setting `save_intermediate = FALSE` - error still occurs
- Running with single config - works with low bootstrap, fails with high bootstrap

## Current Hypothesis
The "cannot open the connection" error appears when:
- Running multiple analyses (20) with high bootstrap iterations (1000)
- Likely caused by resource exhaustion (memory, file handles, or temp files)
- The phylolm package with bootstrap may be creating temporary files or connections

## Recommended Next Steps

### Option 1: Run in Smaller Batches
```r
# Run analyses 1-5 with full bootstrap
batch1 <- run_phyloglm_batch(
  analysis_configs = configs[1:5],
  data = data,
  tree = tree,
  output_dir = "Outputs/PhyloglmResults/batch1",
  bootstrap_n = 1000,
  save_intermediate = FALSE
)

# Continue with remaining batches...
```

### Option 2: Reduce Bootstrap Iterations
```r
# Start with 100-200 iterations
batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  output_dir = "Outputs/PhyloglmResults",
  bootstrap_n = 200,
  save_intermediate = FALSE
)
```

### Option 3: Debug phylolm Bootstrap
- Check if phylolm creates temp files during bootstrap
- Monitor system resources during execution
- Add garbage collection between analyses

## Analysis Configuration Summary
The batch includes 20 analyses:

**Female Song as Response (10 analyses):**
- FS_vs_CB_TerrWS_Mass
- FS_vs_CB_TerrWS_absLat
- FS_vs_CB_TerrWS_PlumDim
- FS_vs_CB_TerrWS_WingDim
- FS_vs_CB_Terr3_Mass
- FS_vs_CB_Migration_Mass
- FS_vs_CB_TerrWS_Region
- FS_vs_CB_TerrWS_Region_Mass (4 predictors)
- FS_vs_CB_Fam_TerrWS
- FS_vs_Fam_TerrWS_Mass

**Cooperative Breeding as Response (10 analyses):**
- CB_vs_FS_TerrWS_Mass
- CB_vs_FS_TerrWS_absLat
- CB_vs_FS_TerrWS_PlumDim
- CB_vs_FS_TerrWS_WingDim
- CB_vs_FS_Terr3_Mass
- CB_vs_FS_Migration_Mass
- CB_vs_FS_TerrWS_Region
- CB_vs_FS_TerrWS_Region_Mass (4 predictors)
- CB_vs_FS_Fam_TerrWS
- CB_vs_Fam_TerrWS_Mass

## Error Location
The error occurs in `run_single_analysis()` function, likely during:
1. Model fitting loop with phyloglm
2. Bootstrap refitting of top models
3. File I/O operations within phyloglm

The exact location needs further debugging to pinpoint.