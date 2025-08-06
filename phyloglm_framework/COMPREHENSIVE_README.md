# Phylogenetic Comparative Analysis of Cooperative Breeding and Female Song Evolution

## Table of Contents
1. [Overview](#overview)
2. [Conceptual Framework](#conceptual-framework)
3. [Methods](#methods)
4. [Implementation Guide](#implementation-guide)
5. [Key Scripts and Outputs](#key-scripts-and-outputs)
6. [Figure Descriptions](#figure-descriptions)
7. [Important Considerations](#important-considerations)
8. [Troubleshooting](#troubleshooting)

## Overview

This framework implements a comprehensive phylogenetic comparative analysis to investigate the evolutionary relationships between cooperative breeding and female song in passerine birds (order Passeriformes). The analysis addresses fundamental questions about social evolution and sexual selection by testing bidirectional evolutionary hypotheses while accounting for phylogenetic non-independence and territorial system moderation.

## Conceptual Framework

### Research Questions

1. **Primary Questions:**
   - Does cooperative breeding drive the evolution of female song in birds?
   - Does female song drive the evolution of cooperative breeding?
   
2. **Secondary Questions:**
   - How do territorial systems (weak/permissive vs. strong/exclusive) moderate these relationships?
   - What other ecological and life-history factors (migration, latitude, sexual dimorphism) influence these evolutionary patterns?
   - Does the directionality of causation differ based on territorial context?

### Theoretical Background

The framework tests competing hypotheses about the relationship between cooperative breeding and female song:

- **Cooperative Breeding → Female Song**: In cooperatively breeding species, females may use song for territory defense, coordination with helpers, or maintaining dominance hierarchies
- **Female Song → Cooperative Breeding**: Species with female song may be pre-adapted for the complex social coordination required in cooperative breeding systems
- **Territorial Moderation**: The strength of territorial boundaries may influence both traits and their interaction

### Analytical Philosophy

The framework employs:
- **Bidirectional testing** to avoid assumptions about evolutionary direction
- **Multi-model inference** to account for uncertainty in model specification
- **Phylogenetic correction** to account for evolutionary non-independence
- **Bootstrap resampling** to provide robust confidence intervals
- **Proper stepwise evaluation** accounting for varying sample sizes across predictors

## Methods

### Data Requirements

1. **Species Trait Database** (`Data_R_2025-06-09.csv`)
   - Binary coding for cooperative breeding presence (0/1)
   - Binary coding for female song presence (0/1)
   - Territorial system classification (weak/permissive vs. strong/exclusive)
   - Body mass (log-transformed)
   - Additional predictors: migration, latitude, dimorphism indices, familial living

2. **Phylogenetic Tree** (`2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex`)
   - Consensus tree for oscine passerines
   - Branch lengths representing evolutionary time
   - Tips matched to species names in trait database

3. **Configuration Files** (YAML format)
   - Define analysis specifications
   - Set model formulations
   - Specify predictor combinations

### Statistical Methods

#### Phylogenetic Generalized Linear Models (PhyloGLM)

We implement phylogenetic logistic regression following Ives & Garland (2010):

```
Y ~ Binomial(p)
logit(p) = Xβ + ε
ε ~ MVN(0, σ²V)
```

Where:
- Y = binary response (female song or cooperative breeding)
- X = design matrix of predictors
- β = regression coefficients
- V = phylogenetic variance-covariance matrix
- σ² = phylogenetic signal parameter

Implementation uses the `phylolm` package (Ho & Ané 2014) with:
- Method: "logistic_MPLE" (most penalized likelihood estimation)
- Bounds: btol = 50, log.alpha.bound = 4
- Optimization: bounded optimization to ensure convergence

#### Model Formulations

For each evolutionary direction, we test a hierarchy of models:

1. **Null Model**: `Response ~ 1`
2. **Main Effects**: `Response ~ Predictor`
3. **Controlled**: `Response ~ Predictor + log(Mass)`
4. **Territorial**: `Response ~ Predictor * Territory`
5. **Full**: `Response ~ Predictor * Territory + log(Mass)`

Where interaction terms test whether effects differ between territorial systems.

#### Model Selection and Averaging

- **AIC-based selection**: Models ranked by Akaike Information Criterion
- **Threshold**: ΔAIC > 2 indicates meaningful difference in support
- **Model averaging**: When multiple models have ΔAIC < 2, parameters averaged weighted by Akaike weights
- **Evidence ratios**: Calculated as exp(ΔAIC/2) to quantify relative support

#### Bootstrap Confidence Intervals

Phylogenetic bootstrap (default n = 1000 iterations):

1. **Resample**: Draw species with replacement maintaining sample size
2. **Prune tree**: Match phylogeny to resampled species  
3. **Refit models**: Estimate parameters on bootstrap sample
4. **Calculate CIs**: Use 2.5% and 97.5% percentiles for 95% CIs

Failed iterations (convergence issues) are excluded from CI calculations.

#### Stepwise Predictor Evaluation

To identify important additional predictors while properly handling varying sample sizes:

1. **Subset data** to species with predictor data available
2. **Refit base model** to the subset (critical for valid comparison)
3. **Add predictor** and fit expanded model on same subset
4. **Calculate AIC improvement**: ΔAICimprovement = AICbase_subset - AICexpanded
5. **Threshold**: Improvement > 2 AIC units considered meaningful

This approach ensures AIC comparisons are valid (same dataset) unlike naive comparison across different sample sizes.

### Predictors Tested

1. **Territory** (ordinal 1-3): Territorial system strength
2. **Migration** (ordinal 1-3): Migratory tendency  
3. **Geographic Region** (categorical): Biogeographic regions
4. **Latitude** (continuous): Absolute value of centroid latitude
5. **Plumage Dimorphism** (continuous): Log male-female plumage difference
6. **Wing Dimorphism** (continuous): Percent wing length dimorphism
7. **Familial Living** (binary): Extended family group living

## Implementation Guide

### Installation and Setup

```r
# Required packages
install.packages(c("phylolm", "ape", "dplyr", "ggplot2", "tidyr", 
                   "patchwork", "yaml", "parallel", "viridis"))

# Clone or download the framework
# Set working directory to CooperativeBreedingEvolution/
setwd("path/to/CooperativeBreedingEvolution")
```

### Running Basic Analysis

```r
# Source the framework
source("phyloglm_framework/batch_runner.R")

# Load data
data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Run analyses for cooperative breeding evolution
configs <- create_standard_configs("comprehensive")
results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = data,
  tree = tree,
  n_bootstrap = 1000,  # Default is 1000 iterations
  output_dir = "Outputs/PhyloglmResults/my_analysis"
)
```

### Running Stepwise Expansion (Corrected Method)

```r
# Source corrected stepwise expansion
source("phyloglm_framework/stepwise_expansion_corrected.R")

# Run with proper AIC comparison
stepwise_results <- run_stepwise_expansion_proper(
  output_dir = "Outputs/PhyloglmResults/stepwise_corrected"
)
```

### Combining Multiple Batch Results

```r
source("phyloglm_framework/combine_batch_results.R")

# Combine results from multiple runs
combine_all_batches(
  batch_dirs = c("batch1", "batch2", "batch3", "batch4"),
  results_path = "Outputs/PhyloglmResults",
  output_name = "combined_all_analyses_batch1-4_boot1000"
)
```

## Key Scripts and Outputs

### Core Framework Scripts

1. **`batch_runner.R`** - Main analysis execution
   - Runs multiple phyloGLM analyses
   - Handles bootstrap iterations (default n=1000)
   - Outputs: Complete results objects, model comparisons, coefficients

2. **`batch_runner_helpers.R`** - Helper functions
   - `prepare_analysis_data()`: Data-tree matching
   - `run_phyloglm_safely()`: Error-handled model fitting
   - `bootstrap_phyloglm()`: Bootstrap implementation
   - `calculate_effect_sizes()`: Odds ratio calculation

3. **`stepwise_expansion_corrected.R`** - Proper stepwise evaluation
   - Tests additional predictors with correct AIC comparison
   - Refits base model to each subset
   - Outputs: Improvement plots, forest plots, summary tables

4. **`create_phyloglm_plots.R`** - Visualization functions
   - Model comparison plots
   - Forest plots for effects
   - Outputs: Publication-quality figures

5. **`combine_batch_results.R`** - Results integration
   - Merges multiple batch runs
   - Creates summary tables
   - Outputs: Combined RDS files, integrated CSVs

### Configuration Files

Configuration files in YAML format define analyses:

```yaml
analyses:
  - name: "FS_vs_CB_TerrWS_Mass"
    response: "FemaleSong_Agg01"
    predictors: 
      - "HighConfidence_Coop"
      - "TerritorialityWeakVsStrong"
    controls:
      - "logMass_AVONET"
    model_complexity:
      - "main"
      - "additive"
      - "interaction"
```

### Output Directory Structure

```
Outputs/PhyloglmResults/
├── [analysis_name]_[timestamp]/
│   ├── all_results.rds                    # Complete R results object
│   ├── analysis_metadata.txt              # Analysis specifications
│   ├── batch_summary.txt                  # Human-readable summary
│   ├── configs/
│   │   └── analysis_configs.rds           # Saved configurations
│   ├── individual_analyses/
│   │   └── [analysis_name]/
│   │       ├── analysis_result.rds        # Individual analysis results
│   │       ├── bootstrap_results.rds      # Bootstrap iterations
│   │       ├── coefficients.csv           # Model coefficients
│   │       ├── convergence_info.csv       # Convergence diagnostics
│   │       ├── data_summary.csv           # Sample sizes and missingness
│   │       ├── effect_sizes.csv           # Odds ratios with CIs
│   │       ├── fitted_models.rds          # Fitted model objects
│   │       ├── model_averaged_effects.csv # Averaged parameters
│   │       ├── model_comparison.csv       # AIC comparison table
│   │       ├── model_formulas.csv         # Model specifications
│   │       ├── prepared_data.rds          # Processed data and tree
│   │       └── plots/
│   │           ├── best_model_coefficients.png
│   │           ├── effect_sizes.png
│   │           └── model_comparison.png
│   ├── integrated_results/
│   │   ├── all_coefficients.csv           # Coefficients across analyses
│   │   ├── all_model_comparisons.csv      # Combined model selection
│   │   ├── analysis_summary_stats.csv     # Overview statistics
│   │   └── effect_comparison.csv          # Comparative effect sizes
│   └── comprehensive_visualizations/
│       ├── all_analyses_summary.png       # Multi-panel summary
│       ├── cooperative_breeding_effects.png
│       └── model_selection_all_analyses.png
└── stepwise_corrected_[date]/
    ├── expansion_results_corrected.rds
    ├── expansion_summary_corrected.txt
    └── improvements_comparison_corrected.png
```

## Figure Descriptions

### Figure 1: Model Comparison Plot (`model_comparison.png`)
**Caption**: Model selection results for phylogenetic generalized linear models testing the relationship between [response variable] and [predictor variables]. Models are ranked by Akaike Information Criterion (AIC), with ΔAIC values calculated relative to the best-supported model (ΔAIC = 0). Points represent AIC values with error bars showing the range across 1000 bootstrap iterations. The shaded region (ΔAIC < 2) indicates models with substantial empirical support that should be considered competitive. Model abbreviations: Null = intercept only; Main = single predictor effect; Add = additive effects; Inter = including interaction term; Mass = controlling for log body mass; Terr = including territorial system; Full = complete model with all predictors and interactions. Sample size: n = [N] species with complete data for all variables.

### Figure 2: Forest Plot of Effect Sizes (`effect_sizes.png`, `best_model_coefficients.png`)
**Caption**: Odds ratios and 95% confidence intervals from phylogenetic generalized linear models of [response variable]. Effect sizes are shown on a logarithmic scale where OR = 1 (vertical dashed line) indicates no effect, OR > 1 indicates positive association, and OR < 1 indicates negative association. Points represent maximum likelihood estimates with horizontal bars showing 95% bootstrap confidence intervals based on 1000 phylogenetic bootstrap iterations. Red symbols indicate statistically significant effects (p < 0.05) while gray symbols indicate non-significant effects. For interaction terms, effects represent the change in the focal relationship between territorial system categories. Estimates are from the best-supported model or model-averaged across models with ΔAIC < 2 when multiple models showed equivalent support.

### Figure 3: Stepwise Predictor Improvement (`improvements_comparison_corrected.png`)
**Caption**: Comparison of model improvements when adding single predictors to base phylogenetic models testing evolutionary relationships between female song and cooperative breeding. Bars show AIC improvement (ΔAIC = AICbase - AICexpanded) where positive values indicate improved model fit. Crucially, base models were refitted to the subset of species with data for each predictor to ensure valid AIC comparison, as AIC values are not comparable across different datasets. Numbers adjacent to bars indicate the sample size (n species) for each comparison. The horizontal dashed line at ΔAIC = 2 represents the threshold for meaningful model improvement. Predictors tested include: Territory (ordinal 1-3 scale of territorial exclusivity), Migration (ordinal 1-3 scale), Wing Dimorphism (percent sexual dimorphism in wing length), Latitude (absolute value), and Familial Living (binary presence of extended family groups). Note the substantially reduced sample size for Familial Living (n ≈ 407) compared to other predictors (n ≈ 825-875).

### Figure 4: Comprehensive Analysis Summary (`all_analyses_summary.png`)
**Caption**: Integrated results from phylogenetic comparative analyses testing bidirectional evolutionary relationships between cooperative breeding and female song across multiple model specifications and datasets. (A) Sample sizes showing the number of species with complete data for each analysis configuration, highlighting variation due to data availability for different predictor combinations. (B) Model selection results showing relative support (ΔAIC) for increasingly complex model formulations, from null models to full models including interactions and controls. Darker colors indicate stronger model support. (C) Primary effect sizes (odds ratios) for the focal predictor across all analyses, with 95% bootstrap confidence intervals. Points are colored by statistical significance and sized by sample size. (D) Interaction effects showing how the relationship between cooperative breeding and female song differs between weak/permissive and strong/exclusive territorial systems, demonstrating context-dependency in evolutionary associations. All analyses account for phylogenetic non-independence using generalized linear models with 1000 bootstrap iterations for uncertainty estimation.

### Figure 5: Master Stepwise Improvements Comparison (`master_improvements_comparison.png`)
**Caption**: Comprehensive evaluation of single-predictor additions across four base model configurations in the phylogenetic analysis framework. Grouped bars compare AIC improvements when adding each predictor to different base models testing: (1) Female Song evolution as a function of Cooperative Breeding with weak vs. strong territory distinction (FS~CB*TerrWS+Mass), (2) Female Song evolution with 3-level territory (FS~CB*Terr3+Mass), (3) Cooperative Breeding evolution as a function of Female Song with territory interaction (CB~FS*TerrWS), and (4) Cooperative Breeding evolution with 3-level territory and mass control (CB~FS*Terr3+Mass). Bar heights represent AIC improvement with positive values indicating better model fit. Only predictors not already included in each base model were tested. Note that earlier versions of this analysis incorrectly compared models fitted to different subsets of species; this figure shows results from naive AIC comparison and should be interpreted with caution. See corrected stepwise analysis for valid comparisons.

### Figure 6: Corrected vs. Uncorrected Stepwise Comparison (if created)
**Caption**: Comparison of stepwise model selection results using incorrect (left) versus correct (right) AIC calculation methods, illustrating the critical importance of proper model comparison in phylogenetic analyses. The incorrect method compares AIC values between models fitted to different numbers of species, creating spurious improvements when predictors have limited data availability. For example, Familial Living appears to improve models by >400 AIC units in the incorrect analysis (left) because it reduces the dataset from 875 to 407 species. The correct method (right) refits the base model to the same subset of species before adding each predictor, revealing that Familial Living actually provides minimal improvement for Female Song models (ΔAIC ≈ -0.3) but substantial improvement for Cooperative Breeding models (ΔAIC ≈ 60.9). This methodological distinction is crucial for valid inference about which predictors truly improve model fit versus those that merely reduce sample size.

## Important Considerations

### Bootstrap Iterations
The framework uses **1000 bootstrap iterations** by default for calculating confidence intervals. This provides robust uncertainty estimates while maintaining computational feasibility. Bootstrap parameters can be adjusted:

```r
results <- run_phyloglm_batch(
  analysis_configs = configs,
  n_bootstrap = 1000,  # Default value
  parallel = TRUE,     # Use parallel processing
  n_cores = 4         # Number of cores for parallel execution
)
```

### Sample Size Effects on AIC
**Critical methodological issue**: AIC values are only comparable between models fitted to identical datasets. The framework includes two approaches:

1. **Incorrect approach** (comprehensive stepwise script): 
   - Compares models with different numbers of species
   - Produces spurious "improvements" when predictors have missing data
   - Example: Familial Living shows 447 AIC improvement (artifact of n=875→407)

2. **Correct approach** (corrected stepwise script):
   - Refits base model to subset of species with predictor data
   - Compares models on identical datasets
   - Example: Familial Living shows -0.3 AIC improvement (true effect)

Always use the corrected approach for valid statistical inference.

### Phylogenetic Signal
Models estimate phylogenetic signal (α parameter) during optimization. Higher α values indicate stronger phylogenetic conservatism in the trait. The bounds (log.alpha.bound = 4) prevent extreme values that can cause convergence issues.

### Convergence Warnings
Some bootstrap iterations may produce convergence warnings:
- "boundary of linear predictor reached" - usually not problematic if infrequent
- Monitor the proportion of successful bootstrap iterations
- Consider increasing btol if warnings are frequent

### Missing Data Handling
- Complete-case analysis: species with missing data for any model variable are excluded
- Sample sizes vary by analysis depending on data availability
- Always report sample sizes for transparency

### Interaction Interpretation
For models with territorial interactions:
- Main effect = effect in reference category (weak/permissive territory)
- Interaction term = difference in effect between territorial systems
- Use model predictions to visualize interaction patterns

## Troubleshooting

### Common Issues and Solutions

1. **"Object 'phy' is not of class 'phylo'"**
   - Check tree format and loading
   - Ensure species names match between tree and data
   - Verify tree is properly imported with `read.nexus()`

2. **"System is computationally singular"**
   - Indicates perfect collinearity
   - Check for redundant predictors
   - Consider simpler model formulations

3. **Bootstrap taking too long**
   - Reduce n_bootstrap for testing (e.g., 100)
   - Enable parallel processing
   - Check for model convergence issues

4. **Different results between runs**
   - Set random seed for reproducibility: `set.seed(12345)`
   - Save complete results objects (.rds files)
   - Document package versions

### Getting Help

1. Check error messages in individual analysis folders
2. Review convergence_info.csv for diagnostic information
3. Examine data_summary.csv for sample size issues
4. Create minimal reproducible examples for debugging

## References

- Ho, L. S. T., & Ané, C. (2014). A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. Systematic Biology, 63(3), 397-408.
- Ives, A. R., & Garland, T. (2010). Phylogenetic logistic regression for binary dependent variables. Systematic Biology, 59(1), 9-26.
- Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526-528.

## Version Information

Framework version: 2.0
Last updated: 2025-07-15
Tested with:
- R version: 4.3.1
- phylolm version: 2.6.2
- ape version: 5.7

## Contact

For questions, bug reports, or contributions, please contact: [contact information]

---

*Note: This framework represents a comprehensive approach to phylogenetic comparative analysis with emphasis on proper statistical methodology, particularly regarding model comparison across datasets of different sizes. Users should carefully consider the methodological implications of their analytical choices.*