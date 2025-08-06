# Stepwise Model Expansion for PhyloGLM

This framework allows you to take your best phylogenetic GLM models and iteratively test adding additional predictors to create optimal "megamodels".

## Quick Start

```r
source("phyloglm_framework/run_stepwise_simple.R")

# Expand best models for Female Song and Cooperative Breeding
results <- run_stepwise_simple()
```

## What It Does

1. **Starts with your best models** from the 20 predefined analyses
   - For Female Song → Cooperative Breeding 
   - For Cooperative Breeding → Female Song

2. **Tests additional predictors** one at a time:
   - Territory (1-3 continuous)
   - Migration (1-3 continuous)  
   - Geographic Region (binary)
   - Absolute Latitude
   - Plumage Dimorphism
   - Wing Dimorphism
   - Familial Living (optional - reduces sample size)

3. **Keeps predictors that improve model fit** (ΔAIC > 2)

4. **Creates visualizations**:
   - Forest plots showing all effects in expanded models
   - Step-by-step improvement plots
   - Predictor importance across analyses

## File Structure

```
phyloglm_framework/
├── stepwise_model_expansion.R    # Core expansion functions
├── stepwise_visualization.R      # Plotting functions
├── run_stepwise_simple.R         # Simple interface
└── example_stepwise_usage.R      # Usage examples
```

## Key Functions

### Simple Interface
- `run_stepwise_simple()` - One-line function to run everything

### Core Functions
- `stepwise_expand_model()` - Expand a single model
- `run_stepwise_batch()` - Expand multiple models
- `get_candidate_predictors()` - Define predictors to test

### Visualization
- `plot_expanded_model_forest()` - Forest plot of all effects
- `plot_stepwise_improvement()` - Show improvement process
- `create_stepwise_summary_figure()` - Combined visualization

## Options

```r
# Include Familial Living (fewer species)
results <- run_stepwise_simple(include_familial = TRUE)

# Skip bootstrap for speed
results <- run_stepwise_simple(bootstrap_n = 0)

# Custom AIC threshold (default = 2)
expanded <- stepwise_expand_model(
  base_result = result,
  aic_threshold = 3  # More stringent
)
```

## Output

Results are saved to: `Outputs/PhyloglmResults/stepwise_expansion_[date]/`

- `all_expansion_results.rds` - Complete R results
- `stepwise_summary.png/pdf` - Main visualization
- `model_comparison.csv` - Table comparing models
- `predictor_importance.png` - Importance across analyses
- Individual analysis files (`*_expanded.rds`, `*_summary.txt`)

## Understanding Results

The expansion process:
1. Starts with base model (e.g., FS ~ CB * Territoriality + Mass)
2. Tests each candidate predictor
3. Adds predictor if AIC improves by > 2
4. Continues with updated model
5. Bootstraps final model for confidence intervals

Example output:
```
FS_vs_CB_TerrWS_Mass
--------------------
Added 2 predictor(s):
  - Territory (AIC improved by 3.4)
  - GeographicRegion_Jetz (AIC improved by 2.8)
Total AIC improvement: 6.2
```

## Customization

See `example_stepwise_usage.R` for:
- Using custom predictor sets
- Expanding specific analyses
- Batch processing multiple models
- Accessing detailed results

## Notes

- Predictors already in the model (or redundant versions) are automatically skipped
- Sample size may decrease as predictors are added due to missing data
- Bootstrap provides confidence intervals for the final expanded model
- The framework checks for various naming patterns to avoid redundancy (e.g., "abs(Latitude)" won't be added if "absLat" is already in the model)