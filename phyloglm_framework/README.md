# PhyloGLM Analysis Framework

A flexible framework for conducting multiple phylogenetic generalized linear model (phyloGLM) analyses with automated formula generation, variable type handling, and comprehensive visualization.

## Overview

This framework provides a systematic approach to:
- Automatically detect and handle different variable types (binary, categorical, continuous)
- Generate appropriate model formulas based on variable types and complexity levels
- Run batch analyses with different configurations
- Create comparative visualizations across analyses
- Handle missing data and phylogenetic tree matching
- Provide robust error handling and progress tracking

## Components

### 1. Variable Classification (`variable_classification.R`)
- Detects variable types automatically
- Handles binary, categorical, and continuous variables
- Provides variable preparation and validation functions

### 2. Formula Builder (`formula_builder.R`)
- Generates model formulas based on variable types
- Supports different complexity levels (null, main, additive, interactions)
- Handles appropriate interaction terms for different variable combinations

### 3. Data Preparation (`data_preparation.R`)
- Enhanced data preparation with automatic variable type handling
- Tree-data matching with species name standardization
- Missing data handling and outlier detection
- Quality checks and data summaries

### 4. Batch Runner (`batch_runner.R` & `batch_runner_helpers.R`)
- Runs multiple analyses with error handling
- Bootstrap confidence intervals for top models
- Model comparison and selection
- Model averaging for parameter estimates
- Organized output directory structure

### 5. Configuration Builder (`config_builder.R`)
- Create analysis configurations programmatically
- Support for YAML configuration files
- Standard configuration templates
- Configuration validation

### 6. Visualization Framework (`visualization_framework.R`)
- Comparative forest plots across analyses
- Model selection heatmaps
- Effect size comparison matrices
- Sample size comparisons
- Interaction visualizations

## Quick Start

```r
# Load framework
source("phyloglm_framework/example_usage.R")

# Create a simple analysis configuration
config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  name = "FemaleSong_CoopBreeding_Analysis"
)

# Run single analysis
result <- run_single_analysis(
  config = config,
  data = data,
  tree = tree,
  output_dir = "outputs/my_analysis"
)

# Or run batch analysis with multiple configurations
configs <- create_standard_configs(include_sets = c("basic", "territorial"))
batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  output_dir = "outputs/batch_analysis"
)
```

## Variable Types

The framework automatically detects and handles:
- **Binary variables**: Presence/absence, yes/no outcomes
- **Categorical variables**: Multi-level factors (e.g., 3-state territory)
- **Continuous variables**: Numeric measurements
- **Special continuous**: Log-transformed variables with special handling

## Model Complexity Levels

- **null**: Intercept-only model
- **main**: Single predictor models
- **additive**: Multiple predictors without interactions
- **twoway**: Two-way interactions
- **threeway**: Three-way interactions

## Output Structure

Each analysis creates an organized output directory:
```
PhyloGLM_Batch_[timestamp]/
├── configs/                    # Saved configurations
├── individual_analyses/        # Results for each analysis
│   └── [analysis_name]/
│       ├── data_summary.csv
│       ├── model_comparison.csv
│       ├── coefficients.csv
│       ├── effect_sizes.csv
│       ├── model_averaged_effects.csv
│       └── plots/
├── integrated_results/         # Combined results across analyses
├── comprehensive_visualizations/
├── all_results.rds            # Complete R results object
└── batch_summary.txt          # Summary report
```

## Advanced Usage

### Using YAML Configuration

```yaml
analyses:
  - name: "FemaleSong_vs_CoopBreeding_Geographic"
    response: "FemaleSong_Agg01"
    predictors: 
      - "HighConfidence_Coop"
      - "GeographicRegion_Jetz"
    controls:
      - "logMass_AVONET"
    complexity_levels:
      - "main"
      - "additive"
      - "twoway"
```

### Custom Transformations

```r
config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop"),
  controls = c("Centroid.Latitude_AVONET"),
  transformations = list(
    Centroid.Latitude_AVONET = "abs"  # Use absolute latitude
  )
)
```

### Parallel Processing

```r
batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  output_dir = output_dir,
  parallel = TRUE,
  n_cores = 4  # Use 4 cores
)
```

## Key Functions

- `create_analysis_config()`: Create single analysis configuration
- `create_standard_configs()`: Generate standard analysis sets
- `run_single_analysis()`: Run one phyloGLM analysis
- `run_phyloglm_batch()`: Run multiple analyses
- `create_comparative_forest_plot()`: Compare effects across analyses
- `create_batch_summary_figure()`: Create comprehensive summary visualization

## Requirements

- R packages: `ape`, `phylolm`, `dplyr`, `ggplot2`, `tidyr`, `patchwork`, `viridis`, `yaml`
- Phylogenetic tree in Nexus format
- Data in CSV format with species names matching tree tips