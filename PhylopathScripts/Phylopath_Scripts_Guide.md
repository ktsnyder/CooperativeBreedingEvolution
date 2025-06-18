# Phylopath Scripts Guide

**Date:** June 8, 2025  
**Purpose:** Complete reference of all scripts used in phylopath analyses, including functions defined and used in each script.

---

## Core Phylopath Functions

### 1. `run_phylopath_fxns.R`
**Location:** Main directory  
**Purpose:** Core functions for running phylopath analyses and downsampling

#### Functions Defined:
- `run_CB_FS_Terr_phylopath()` - Run single phylopath analysis with CB, FS, Territoriality, and optional mass variable
- `run_multiple_phylopath()` - Run multiple phylopath analyses with downsampling
- `downsample_run_phylopath()` - Run single downsampled phylopath analysis
- `aggregate_by_seed()` - Aggregate detailed model results by seed
- `convert_all_seed_results_to_dataframes()` - Convert aggregated results to dataframes
- `test_sampleRegion_function()` - Test function for region sampling
- `summarize_models_by_seed()` - Summarize model selection by seed
- `average_models_conditional()` - Calculate conditional model averages
- `average_models_full()` - Calculate full model averages
- `extract_model_info()` - Extract model information from phylopath results

#### Functions Used:
- `phylopath::define_model_set()` - Define phylopath models
- `phylopath::phylo_path()` - Run phylopath analysis
- `phylopath::best()` - Get best model
- `phylopath::average()` - Average models
- `phylopath::plot()` - Plot DAGs
- `phytools::keep.tip()` - Subset phylogenetic tree
- `ggplot2` functions for plotting
- `dplyr` functions for data manipulation

---

### 2. `execute_run_phylopath_fxns.R`
**Location:** Main directory (identical copy in Box_Scripts/)  
**Purpose:** Example execution code for running phylopath analyses

#### Functions Defined:
- None (execution script)

#### Functions Used:
- `run_multiple_phylopath()` - From run_phylopath_fxns.R
- `aggregate_by_seed()` - From run_phylopath_fxns.R
- `convert_all_seed_results_to_dataframes()` - From run_phylopath_fxns.R
- Various plotting and data manipulation functions

---

## Enhanced/Fixed Functions

### 3. `phylopath_plotting_comprehensive_fixed.R`
**Location:** claude_code_sessions/  
**Purpose:** Fixed plotting functions with proper file naming, DAG layouts, and colors

#### Functions Defined:
- `create_nondownsampled_plots()` - Create all plots for non-downsampled analyses
- `create_downsampled_plots()` - Create all plots for downsampled analyses
- `create_all_phylopath_plots()` - Wrapper function for all plot types

#### Functions Used:
- `phylopath::best()` - Get best model
- `phylopath::average()` - Average models
- `phylopath::plot()` - Plot DAGs
- `phylopath::summary()` - Summarize phylopath results
- `ggplot2` functions for custom plotting
- `cowplot::plot_grid()` - Combine plots
- `aggregate_by_seed()` - From run_phylopath_fxns.R
- `convert_all_seed_results_to_dataframes()` - From run_phylopath_fxns.R

---

### 4. `downsample_run_phylopath_fixed.R`
**Location:** claude_code_sessions/  
**Purpose:** Fixed version of downsample_run_phylopath that handles HaveData column correctly

#### Functions Defined:
- `downsample_run_phylopath_fixed()` - Fixed downsampling function

#### Functions Used:
- `run_CB_FS_Terr_phylopath()` - From run_phylopath_fxns.R

---

### 5. `create_phylopath_plots_enhanced_fixed.R`
**Location:** claude_code_sessions/  
**Purpose:** Enhanced plotting for detailed model files without renaming

#### Functions Defined:
- `create_enhanced_phylopath_plots()` - Create enhanced plots from detailed_models CSV files
- `process_all_detailed_models()` - Process all detailed model files in directory

#### Functions Used:
- `ggplot2` functions for plotting
- `dplyr` functions for data manipulation
- `cowplot::plot_grid()` - Combine plots

---

## Downsampling Calculation Scripts

### 6. `calculate_downsampling_function.R`
**Location:** claude_code_sessions/  
**Purpose:** Calculate territoriality bias downsampling numbers

#### Functions Defined:
- `calculate_downsampling()` - Calculate number of species to remove for territoriality bias

#### Functions Used:
- Base R functions only

---

### 7. `calculate_stratified_downsampling.R`
**Location:** Main directory  
**Purpose:** Calculate all stratified downsampling scenarios (geographic, cooperative)

#### Functions Defined:
- `calculate_stratified_downsampling()` - Calculate downsampling for specific trait combinations
- `generate_downsampling_report()` - Generate markdown report

#### Functions Used:
- `dplyr` functions for data manipulation
- `phytools::read.nexus()` - Read phylogenetic tree

---

### 8. `calculate_stratified_downsampling_with_territoriality.R`
**Location:** claude_code_sessions/  
**Purpose:** Extended version including territoriality calculations

#### Functions Defined:
- `calculate_territoriality_downsampling()` - Calculate territoriality-specific downsampling
- `calculate_stratified_downsampling()` - Extended version with territoriality
- `generate_downsampling_report()` - Enhanced report generator

#### Functions Used:
- `dplyr` functions for data manipulation
- `phytools::read.nexus()` - Read phylogenetic tree

---

## Test and Integration Scripts

### 9. `test_phylopath_comprehensive_fixed.R`
**Location:** claude_code_sessions/  
**Purpose:** Comprehensive test of all phylopath analyses with n=50 iterations

#### Functions Defined:
- None (test script)

#### Functions Used:
- `run_CB_FS_Terr_phylopath()` - From run_phylopath_fxns.R
- `run_multiple_phylopath()` - From run_phylopath_fxns.R
- `create_all_phylopath_plots()` - From phylopath_plotting_comprehensive_fixed.R
- `calculate_downsampling()` - From calculate_downsampling_function.R
- `calculate_stratified_downsampling()` - From calculate_stratified_downsampling_with_territoriality.R
- `create_enhanced_phylopath_plots()` - From create_phylopath_plots_enhanced_fixed.R
- `downsample_run_phylopath_fixed()` - From downsample_run_phylopath_fixed.R

---

### 10. `phylopath_additions_for_Run_Analyses.R`
**Location:** claude_code_sessions/  
**Purpose:** Code snippets to be added to Run_Analyses.R

#### Functions Defined:
- None (code snippets)

#### Functions Used:
- All functions from above scripts (demonstrates usage)

---

## Function Call Hierarchy

### For Non-Downsampled Analysis:
```
Run_Analyses.R
├── run_CB_FS_Terr_phylopath() [from run_phylopath_fxns.R]
│   ├── phylopath::define_model_set()
│   ├── phylopath::phylo_path()
│   └── phylopath plotting functions
└── create_all_phylopath_plots() [from phylopath_plotting_comprehensive_fixed.R]
    └── create_nondownsampled_plots()
        ├── phylopath::best()
        ├── phylopath::average()
        └── phylopath::plot()
```

### For Downsampled Analysis:
```
Run_Analyses.R
├── calculate_downsampling() [from calculate_downsampling_function.R]
├── run_multiple_phylopath() [from run_phylopath_fxns.R]
│   └── downsample_run_phylopath_fixed() [overrides downsample_run_phylopath]
│       └── run_CB_FS_Terr_phylopath()
└── create_all_phylopath_plots() [from phylopath_plotting_comprehensive_fixed.R]
    └── create_downsampled_plots()
        ├── aggregate_by_seed()
        ├── convert_all_seed_results_to_dataframes()
        └── ggplot2 plotting functions
```

### For Report Generation:
```
Run_Analyses.R
└── calculate_stratified_downsampling() [from calculate_stratified_downsampling_with_territoriality.R]
    ├── calculate_territoriality_downsampling()
    └── generate_downsampling_report()
```

---

## Required Source Order in Run_Analyses.R

```r
# 1. Core functions
source("run_phylopath_fxns.R")

# 2. Fixed/enhanced functions
source("claude_code_sessions/phylopath_plotting_comprehensive_fixed.R")
source("claude_code_sessions/downsample_run_phylopath_fixed.R")
source("claude_code_sessions/create_phylopath_plots_enhanced_fixed.R")
source("claude_code_sessions/calculate_downsampling_function.R")

# 3. Override function
downsample_run_phylopath <- downsample_run_phylopath_fixed

# 4. For report generation (optional)
source("claude_code_sessions/calculate_stratified_downsampling_with_territoriality.R")
```

---

## Output Files Generated

### From Non-Downsampled Analyses:
- `phylopath [variables] [n]species best_model.png`
- `phylopath [variables] [n]species full_avg.png`
- `phylopath [variables] [n]species cond_avg.png`
- `phylopath [variables] [n]species cicbar.png`
- `phylopath [variables] [n]species summary.png`
- `phylopath [variables] [n]species dags_combined.png`

### From Downsampled Analyses:
- `[variables] Remove[n][group] downsample_info.png`
- `[variables] Remove[n][group] violin.png`
- `[variables] Remove[n][group] heatmap.png`
- `[variables] Remove[n][group] model_freq.png`
- `detailed_models_[scenario]_n[iterations]_[date].csv`
- `conditional_average_plots_[scenario]_[iterations]_[date].pdf`

### From Enhanced Plotting:
- `enhanced_[scenario]_models.png`
- `enhanced_[scenario]_edges.png`
- `enhanced_[scenario]_combined.png`

### From Report Generation:
- `Stratified_Downsampling_Calculations.md`

---

## OBSOLETE Scripts (DO NOT USE)

### Main Directory Obsolete Scripts:
1. **`scratch PhyloPath.R`** - Early testing/scratch work
2. **`run_phylopath_downsampled.R`** - Replaced by run_multiple_phylopath() in run_phylopath_fxns.R
3. **`test_phylopath_downsampled.R`** - Old test script
4. **`run_phylopath_additional.R`** - Replaced by comprehensive functions
5. **`test_phylopath_additional.R`** - Old test script
6. **`downsample_run_phylopath_fixed.R`** - USE claude_code_sessions/ version instead
7. **`create_phylopath_plots_enhanced.R`** - USE claude_code_sessions/create_phylopath_plots_enhanced_fixed.R instead

### Box_Scripts Obsolete Scripts:
1. **`Box_Scripts/phylopath_continuousTerritory.R`** - Old implementation
2. **`Box_Scripts/PhyloPath_tests.R`** - Old tests
3. **`Box_Scripts/phylopath_replaceBinaryResponseVariable.R`** - Not needed with current implementation
4. **`Box_Scripts/phylopath_multipleSets.R`** - Functionality integrated into main functions
5. **`Box_Scripts/run_phylopath_fxns.R`** - Identical to main directory version (use main)
6. **`Box_Scripts/execute_run_phylopath_fxns.R`** - Identical to main directory version (use main)

### claude_code_sessions Obsolete Scripts:
1. **`phylopath_integration_code.R`** - Early integration attempt
2. **`test_phylopath_integration.R`** - Replaced by test_phylopath_comprehensive_fixed.R
3. **`phylopath_plotting_comprehensive.R`** - USE phylopath_plotting_comprehensive_fixed.R instead
4. **`test_phylopath_comprehensive_n50.R`** - USE test_phylopath_comprehensive_fixed.R instead

---

## Detailed Source Explanations for Run_Analyses.R

### 1. `source("run_phylopath_fxns.R")`
**What it does:**
- Provides the core phylopath analysis functions
- `run_CB_FS_Terr_phylopath()`: Runs a single phylopath analysis testing relationships between cooperative breeding (CB), female song (FS), territoriality, and an optional continuous variable (mass/dichromatism/dimorphism)
- `run_multiple_phylopath()`: Manages multiple iterations of downsampled analyses, calling downsample_run_phylopath() n times with different random seeds
- Contains helper functions for model aggregation and summarization across multiple runs
- Handles all the phylopath model definitions (48 models testing different causal relationships)

### 2. `source("claude_code_sessions/phylopath_plotting_comprehensive_fixed.R")`
**What it does:**
- Provides enhanced plotting functions that fix label cutoff issues and implement proper file naming
- `create_nondownsampled_plots()`: Creates 6 plot types for single analyses including DAGs with proper margins and spacing
- `create_downsampled_plots()`: Creates 4 plot types for downsampled analyses including violin plots and heatmaps
- Implements the file naming convention: non-downsampled files include all variables and species count; downsampled files include "Remove[n][group]" format
- Sets proper plot dimensions (12x10 for individual DAGs, 24x10 for combined)

### 3. `source("claude_code_sessions/downsample_run_phylopath_fixed.R")`
**What it does:**
- Provides a fixed version of the downsampling function that correctly handles the HaveData column
- Ensures compatibility with different column naming conventions (HaveData, HaveFSData, HaveCBData)
- Adds geographic regions if not present in the data
- Performs the actual random sampling to remove specified numbers of species from target groups
- This MUST be followed by `downsample_run_phylopath <- downsample_run_phylopath_fixed` to override the original

### 4. `source("claude_code_sessions/create_phylopath_plots_enhanced_fixed.R")`
**What it does:**
- Processes the detailed_models CSV files generated by downsampled analyses
- Creates enhanced visualizations: model frequency plots, edge coefficient plots, and combined views
- Preserves original model names (no renaming) in frequency plots
- Summarizes which models appear most frequently across iterations
- Shows edge coefficients with significance and direction information

### 5. `source("claude_code_sessions/calculate_downsampling_function.R")`
**What it does:**
- Provides the `calculate_downsampling()` function specifically for territoriality bias correction
- Calculates how many species need to be removed from high vs low territoriality groups to balance data availability proportions
- Works with any binary territoriality variable (TerritorialityWeakVsStrong or Territory_12vs3)
- Returns the number to remove and which group to downsample from
- Used during phylopath runs to determine territoriality bias correction parameters

### 6. `source("claude_code_sessions/calculate_stratified_downsampling_with_territoriality.R")` (OPTIONAL)
**What it does:**
- Generates the comprehensive Stratified_Downsampling_Calculations.md report
- Calculates all downsampling scenarios: geographic (83 species), tropical cooperative (24 species), global cooperative (15 species), and territoriality biases
- Extended version that includes territoriality calculations for both TerritorialityWeakVsStrong and Territory_12vs3
- Only needed when you want to regenerate the downsampling calculations report
- Not required for running the actual phylopath analyses

---

## Output Files - Which Script Generates Each

### From Non-Downsampled Analyses:
**Generated by:** `create_nondownsampled_plots()` in `phylopath_plotting_comprehensive_fixed.R`
- `phylopath [variables] [n]species best_model.png`
- `phylopath [variables] [n]species full_avg.png`
- `phylopath [variables] [n]species cond_avg.png`
- `phylopath [variables] [n]species cicbar.png`
- `phylopath [variables] [n]species summary.png`
- `phylopath [variables] [n]species dags_combined.png`

### From Downsampled Analyses:
**Generated by:** `create_downsampled_plots()` in `phylopath_plotting_comprehensive_fixed.R`
- `[variables] Remove[n][group] downsample_info.png`
- `[variables] Remove[n][group] violin.png`
- `[variables] Remove[n][group] heatmap.png`
- `[variables] Remove[n][group] model_freq.png`

**Generated by:** `run_multiple_phylopath()` in `run_phylopath_fxns.R`
- `detailed_models_[scenario]_n[iterations]_[date].csv`
- `conditional_average_plots_[scenario]_[iterations]_[date].pdf`

### From Enhanced Plotting:
**Generated by:** `create_enhanced_phylopath_plots()` in `create_phylopath_plots_enhanced_fixed.R`
- `enhanced_[scenario]_models.png`
- `enhanced_[scenario]_edges.png`
- `enhanced_[scenario]_combined.png`

### From Report Generation:
**Generated by:** `calculate_stratified_downsampling()` in `calculate_stratified_downsampling_with_territoriality.R`
- `Stratified_Downsampling_Calculations.md`

---

*This guide covers all essential scripts for phylopath analyses as of June 8, 2025*