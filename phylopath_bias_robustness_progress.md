# Phylopath Bias Robustness Visualization Progress

## Overview
This document tracks the progress of creating publication-ready figures showing the robustness of the cooperative breeding-female song association to various bias corrections.

## Key Scripts and Functions

### Main Script
- `create_phylopath_bias_robustness_figure.R` - Contains all visualization functions

### Key Functions Created
1. `extract_phylopath_results()` - Extracts results from CSV files in PhylopathDownsampled directory
2. `create_bias_robustness_figure()` - Original 3-panel figure with forest plot and model frequencies
3. `create_coefficient_violin_plot()` - Violin plots showing coefficient distributions
4. `create_coefficient_heatmap()` - Heatmap of coefficients across models and bias corrections
5. `create_model_frequency_figure()` - Faceted bar plots of top models
6. `create_coefficient_ridge_plot()` - Ridge plots showing distributions by model structure
7. `create_comprehensive_robustness_figure()` - Combined multi-panel figure

## What Has Worked Well
1. **Data extraction** - Successfully extracts data from all bias corrections including dimorphism subdirectories
2. **Forest plot** - Shows coefficients with CIs clearly (though missing some error bars)
3. **Violin plots** - Good visualization of coefficient distributions
4. **Color scheme** - Consistent blue theme (#2E86AB) works well

## Current Issues and Areas for Improvement

### High Priority Issues
1. **Panel positioning** - Panel B frequently cut off on right side
2. **Missing DAG plot** - Need to show the actual path model structure
3. **Species counts missing** - No sample sizes shown for each bias correction
4. **Violin plot reference lines** - Single line spanning all violins doesn't make sense
5. **Missing error bars** - Some bias corrections show no error bars in forest plot

### Medium Priority Issues
1. **Heatmap problems** - Many NaN values, shows non-top models
2. **Panel redundancy** - Forest plot and violin plots show similar information
3. **Y-axis spacing** - Weird gaps between labels and plots
4. **Panel C not impactful** - Model frequency plot needs improvement

### Data Issues
- Wing and Plumage Dimorphism data exists but may have different coefficient columns
- Need to verify actual sample sizes after downsampling

## Existing Plots Analysis

### What's Already Available in PhylopathDownsampled
1. **violin-box_coefficients.png** - Shows ALL path coefficients with proper reference lines
   - Horizontal violin/box plots for each path
   - Red vertical lines show full dataset conditional averages
   - Much better than my implementation!

2. **model_freq.png** - Clear model frequency visualization  
   - Horizontal bar chart
   - Shows count and percentage
   - Only shows models with ΔCICc < 2

3. **downsample_info.png** - Data availability visualization
   - Stacked bar chart showing data vs missing
   - Shows sample sizes and percentages
   - Clear target proportion annotation

4. **conditional_average_plots PDFs** - Likely contains DAG visualizations

### Key Insights from Existing Plots
- They show ALL paths, not just CB→FS
- Reference lines are correctly placed as vertical lines on horizontal plots
- Clear labeling with human-readable path names
- Consistent blue color scheme
- Sample sizes included

## Ideas for Improvement

### Based on User's Figure Caption Example
1. **Panel A**: DAG with conditional-averaged path coefficients
   - Show all paths: CB→FS, T→FS, T→CB, Mass effects
   - Arrow thickness proportional to coefficient strength
   - Only include models with ΔCICc < 2

2. **Panel B**: Model selection results
   - Bar plot of CICc values for top 20 models
   - Color code: blue for ΔCICc < 2, gray for ≥ 2
   - Show ΔCICc values as labels

3. **Panels C-E**: Path coefficient distributions
   - Separate violin plots for each bias correction
   - Show ALL path coefficients, not just CB→FS
   - Add reference lines for full dataset values

### New Visualization Ideas
1. **Sankey diagram** - Show flow of model selection across bias corrections
2. **Paired plots** - Before/after downsampling compositions
3. **Network plot** - Show relationships between variables with edge weights
4. **Coefficient stability plot** - Line plot showing how coefficients change across iterations

## Recent Progress (Update 1)

### New Functions Added
1. **`create_phylopath_dag()`** - Creates DAG plot from phylopath results
   - Uses conditional averaging
   - Proper layout for 4 variables
   - Can load from RDS file or use existing result object

2. **`get_species_count()`** - Extracts species count from detailed models

3. **`create_all_paths_coefficient_plot()`** - Shows ALL path coefficients
   - Faceted by bias correction
   - Horizontal violin plots like existing examples
   - Shows distributions for all paths, not just CB→FS

### Key Discoveries
- Existing plots in PhylopathDownsampled are much better than my initial attempts
- They show ALL path coefficients, not just CB→FS
- Use horizontal orientation for better readability
- Include proper reference lines as vertical lines

### Issues to Address
1. Need to find actual phylopath result RDS files for full dataset
2. Species counts still missing from data
3. Panel sizing issues persist
4. Need to create proper model selection bar plot

## Recent Progress (Update 2)

### Successfully Created Functions
1. **`create_forest_plot_with_counts()`** - Forest plot WITH species counts!
   - Shows n = xxx for each bias correction
   - Clean layout with reference line at 0.56
   - All bias corrections now display properly including dimorphism

2. **`create_model_selection_plot()`** - Bar plot of top models
   - Shows CICc values and Δ CICc
   - Color codes models with Δ CICc < 2
   - Clean horizontal layout

3. **`create_all_paths_coefficient_plot()`** - Shows ALL paths (not tested yet)

### Key Improvements
- Species counts now extracted from nSpecies column in detailed_models
- Dimorphism data now displays correctly
- Better function modularity - can create individual plots
- Consistent styling across all plots

### Still To Do
1. Test DAG function with actual phylopath results
2. Create final multi-panel figure with best combination
3. Fix panel spacing issues
4. Test all paths coefficient plot

## Recent Progress (Update 3)

### Major Accomplishments
1. **Created publication-ready figure function** - `create_publication_robustness_figure()`
   - Matches user's figure caption structure
   - Panel A: DAG placeholder (needs full dataset result)
   - Panel B: Model selection bar plot  
   - Panels C-E: Violin plots for specific bias corrections
   - Red reference lines for full dataset values

2. **All individual plot functions now working**:
   - Forest plot with species counts ✓
   - Model selection bar plot ✓
   - Violin plots with proper reference lines ✓
   - DAG function created (needs testing) ✓

### Complete Function List

#### Main Functions
1. `extract_phylopath_results()` - Extracts results from CSV files, including dimorphism subdirectories
2. `create_bias_robustness_figure()` - Original 3-panel figure (needs updating)
3. `create_publication_robustness_figure()` - NEW! Matches user's caption structure

#### Individual Plot Functions  
1. `create_forest_plot_with_counts()` - Forest plot with species counts
2. `create_model_selection_plot()` - Bar plot of model CICc values
3. `create_coefficient_violin_plot()` - Violin plots for single coefficient
4. `create_all_paths_coefficient_plot()` - Shows ALL path coefficients
5. `create_phylopath_dag()` - Creates DAG from phylopath results
6. `create_coefficient_heatmap()` - Heatmap of coefficients (needs fixing)
7. `create_coefficient_ridge_plot()` - Ridge plots by model structure
8. `create_comprehensive_robustness_figure()` - Multi-panel with various elements

#### Helper Functions
1. `get_species_count()` - Extracts species count from data

### What Works Well
- Species counts now display correctly (n = xxx)
- All bias corrections including dimorphism show data
- Model selection plot matches existing style
- Individual functions allow flexible plotting
- Consistent color scheme (#2E86AB)

### Remaining Issues
1. Panel B still gets cut off in some multi-panel layouts
2. Need actual phylopath result for DAG
3. Some panel spacing issues persist
4. Heatmap shows too many NaN values

### Usage Examples
```r
# Load the script
source('create_phylopath_bias_robustness_figure.R')

# Create publication figure
pub_fig <- create_publication_robustness_figure()

# Create individual plots
forest <- create_forest_plot_with_counts()
models <- create_model_selection_plot(detailed_models_df)

# Extract results for custom analysis
results <- extract_phylopath_results()
```

## Recent Progress (Update 4 - Final)

### Latest Addition
- **`create_all_paths_horizontal_plot()`** - Replicates existing style perfectly!
  - Shows ALL path coefficients in horizontal violin/box format
  - Matches the style of plots in PhylopathDownsampled directory
  - Ordered by mean coefficient value
  - Clean, professional appearance

### Summary of Available Visualizations

1. **Forest Plots**
   - `create_forest_plot_with_counts()` - Shows CB→FS coefficients with species counts

2. **Model Selection** 
   - `create_model_selection_plot()` - Bar plot of CICc values for top models

3. **Coefficient Distributions**
   - `create_coefficient_violin_plot()` - Vertical violins for single coefficient
   - `create_all_paths_horizontal_plot()` - Horizontal violins for ALL paths (best!)
   - `create_all_paths_coefficient_plot()` - Faceted version

4. **Multi-panel Figures**
   - `create_publication_robustness_figure()` - Matches user's caption structure
   - `create_comprehensive_robustness_figure()` - Alternative layout
   - `create_bias_robustness_figure()` - Original 3-panel version

5. **Other Visualizations**
   - `create_phylopath_dag()` - DAG plot (needs phylopath result)
   - `create_coefficient_heatmap()` - Heatmap (needs improvement)
   - `create_model_frequency_figure()` - Faceted model frequencies

### Best Practices Discovered
1. Use horizontal orientation for path coefficient plots
2. Show ALL paths, not just CB→FS
3. Include species counts (n = xxx)
4. Use consistent color scheme (#2E86AB)
5. Add reference lines for full dataset values
6. Keep titles and labels clear and concise

### Final Recommendations
For publication-ready figures, use:
1. `create_all_paths_horizontal_plot()` for showing coefficient distributions
2. `create_model_selection_plot()` for model comparison
3. `create_forest_plot_with_counts()` for focused CB→FS comparison
4. `create_publication_robustness_figure()` for multi-panel figure

## Recent Progress (Update 5 - Latest Session)

### Major Accomplishments

1. **Created Model Consistency Heatmap Function** - `create_model_consistency_heatmap()`
   - Shows how often each model appears in top models across bias corrections
   - Colors indicate percentage of iterations where model had ΔCICc < 2
   - Numbers show count of iterations where model was best
   - Automatically sizes plot based on number of models
   - Clean, publication-ready output

2. **Successfully Ran Full Dataset Phylopath Analysis**
   - Created `run_full_dataset_phylopath.R` script
   - Generated phylopath results for 875 species
   - Extracted key coefficients:
     - **CB→FS: 0.556** (this is our reference value!)
     - FS→CB: 0.391
     - TERR→FS: 1.215
     - TERR→CB: 0.354
   - Saved results as RDS object and created DAG visualization

### Key Findings from Model Consistency Heatmap
- Model Y2_TERR→FS_FS→COOP_MASS→FS_MASS→COOP appears most consistently across bias corrections
- Models X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS and Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS also very consistent
- Shows clear pattern of which model structures are robust to different biases

## Recent Progress (Update 6 - Final Session Accomplishments)

### Improvements to Model Consistency Heatmap
1. **Added stars (★)** to mark models that were best in full dataset analysis
2. **Fixed color scale** - now uses 2-color gradient (light blue to dark blue) for better readability
3. **Scale now goes to 500** as requested
4. **Models ordered** with full dataset best models at the top
5. **Cleaner appearance** with black text on lighter backgrounds

### Created Enhanced DAG Function
- `create_enhanced_dag()` - Produces cleaner, more professional DAG visualizations
- Uses human-readable labels (e.g., "Cooperative Breeding" instead of "HighConfidence_Coop")
- Shows path coefficients with proper arrow directions
- Edge thickness represents coefficient magnitude
- Blue arrows for positive effects
- Includes sample size in subtitle

### Created Final Publication Figures
1. **Individual high-quality plots**:
   - `final_enhanced_dag.png` - Clean DAG from full dataset
   - `final_model_consistency_heatmap.png` - Improved heatmap with all requested features
   - `final_forest_plot.png` - Forest plot with species counts
   - `final_all_paths_plot.png` - Horizontal violin plots for all paths

2. **Comprehensive multi-panel figure**:
   - `final_comprehensive_robustness_figure.png/pdf` - 4-panel publication figure
   - Top row: DAG and model consistency heatmap
   - Bottom row: Forest plot and path coefficient distributions
   - Shows robustness of CB→FS association across all analyses

### Key Findings Summary
- CB→FS coefficient in full dataset: **0.556**
- Across bias corrections: ranges from 0.353 to 0.619
- Models X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS and Y2_TERR→FS_FS→COOP_MASS→FS_MASS→COOP consistently selected
- Association remains positive and substantial regardless of bias correction applied

## Recent Progress (Update 7 - DAG Improvements)

### Enhanced DAG Visualization Issues Fixed
1. **Label positioning for bidirectional edges** - Fixed sign cancellation issue
   - Problem: Both labels appeared at same position due to perpendicular vector calculations
   - Solution: Ensured perpendicular vectors point in same direction for bidirectional edge pairs
   - Now both CB→FS (0.56) and FS→CB (0.39) labels are visible and properly positioned

2. **Standalone DAG script created**
   - `create_enhanced_DAG.R` - Dedicated script for DAG visualization
   - Fixed curvature values for bidirectional arrows (both now use -0.2)
   - Improved arrow head visibility with larger size (0.4 cm)
   - Increased node radius to 0.7 to ensure arrows don't overlap nodes
   - Fixed Delta symbol encoding for PDF output using `bquote()` and plotmath

### Technical Details of Fixes
- **Label positioning**: Added code to check if perpendicular vectors point in opposite directions and flip one if needed
- **PDF encoding**: Uses `bquote()` with mathematical Delta symbol to avoid Unicode issues
- **Output formats**: Saves both PDF (600 dpi) and PNG (300 dpi) versions

### Current Limitations
- **Not fully flexible**: The `create_enhanced_DAG.R` script is functional for this specific analysis but has hardcoded assumptions:
  - Expects exactly these four variables: logMass_AVONET, TerritorialityWeakVsStrong, HighConfidence_Coop, FemaleSong_Agg01
  - Node positions are fixed for this specific 4-variable layout
  - Bidirectional edge handling is specifically coded for CB↔FS relationship
  - Would need modification to handle alternative variables or different bidirectional edge configurations
- For analyses with different variables or edge patterns, the node positions, variable names, and bidirectional edge logic would need to be updated

## File Locations
- Results: `/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathDownsampled/`
- Scripts: `/Users/kate/Desktop/CooperativeBreedingEvolution/`
  - Main visualization functions: `create_phylopath_bias_robustness_figure.R`
  - Standalone DAG function: `create_enhanced_DAG.R`
  - Full dataset analysis: `run_full_dataset_phylopath.R`
  - Final figure creation: `create_final_publication_figure.R`
- Test outputs: Various `test_*.png` files in main directory
- Full dataset results:
  - `phylopath_full_dataset_result.rds` - R object with full analysis
  - `phylopath_full_dataset_dag.png` - DAG visualization
  - `phylopath_full_dataset_coefficients.csv` - Path coefficients

## Recent Progress (Update 8 - Model Consistency Heatmap Refinement)

### Standalone Model Consistency Heatmap Script Created
1. **`create_bias_model_consistency_heatmap.R`** - Dedicated script for model consistency visualization
   - Moved from main visualization functions file to standalone script
   - Produces publication-ready heatmap showing model selection consistency across bias corrections
   - Automatically extracts full dataset model information and delta CICc values
   - Orders models by delta CICc from full dataset analysis

### Key Features of Final Implementation
1. **Three-panel layout**:
   - Top panel: Heatmap showing iterations where each model had ΔCICc < 2
   - Middle panel: ΔCICc values from full dataset analysis
   - Bottom panel: Rate presence/absence matrix showing model structure

2. **Visual improvements**:
   - White cells for models never selected (ΔCICc never < 2)
   - Dynamic text color (white on dark cells, black on light)
   - Stars (★) mark models from full dataset analysis
   - Comprehensive legend with all symbols explained
   - Clean margins and spacing

3. **Technical refinements**:
   - Changed from `expression()` to `bquote()` for cleaner Delta symbol rendering
   - Used dummy aesthetic method to add white square, dot, and star to legend
   - Optimized margins and spacing for compact layout
   - Subtitle and axis labels positioned for clarity

### Publication Readiness
- The function produces high-quality outputs that require only minor editing in Illustrator
- Both PDF (600 dpi) and PNG (300 dpi) versions are generated
- Figure dimensions automatically adjust based on number of models and bias corrections
- All text elements properly render in both formats using `bquote()`

### Figure Caption
"Model consistency across bias corrections. Heatmap showing the frequency with which each phylogenetic path model was selected (ΔCICc < 2) across 500 iterations under seven different bias correction schemes. Numbers in cells indicate iterations where that model had the lowest CICc. Models are ordered by their ΔCICc values from the full dataset analysis (shown between panels). Bottom panel shows the presence/absence of seven possible evolutionary transitions in each model. Stars indicate models that were within ΔCICc < 2 in the full dataset analysis (n = 875 species). White cells indicate models never selected under that bias correction."

## Recent Progress (Update 9 - Forest Plot Function Separated)

### New Standalone Script: `create_forest_plot_with_counts.R`

A dedicated forest plot function has been created and moved to its own script for better modularity. This function creates publication-ready forest plots showing path coefficients across bias corrections.

#### Key Features:
1. **Proper conditional averaging** - Uses phylopath's conditional averaging approach with CICc weights
2. **Automatic reference value loading** - Loads values from `phylopath_full_dataset_result.rds` if not provided
3. **Flexible rate selection** - Supports all 7 possible paths (COOP→FS, FS→COOP, TERR→FS, etc.)
4. **Smart file naming** - Automatically generates descriptive filenames based on selected rate
5. **Custom ordering** - Bias corrections ordered to match other figure panels
6. **Clear labeling** - Shows "n ~ XXX species" with resampling clarification
7. **PDF compatibility** - Uses tilde (~) instead of ≈ for proper PDF rendering

#### Usage Examples:
```r
# Load required functions
source('create_phylopath_bias_robustness_figure.R')  # For extract_phylopath_results()
source('create_forest_plot_with_counts.R')

# Simplest usage - all defaults
create_forest_plot_with_counts()
# Creates: Outputs/Figures/forest_plot_COOP_FS_with_counts.pdf

# Different rate
create_forest_plot_with_counts(rate_to_plot = "FS->COOP")
# Creates: Outputs/Figures/forest_plot_FS_COOP_with_counts.pdf

# Custom output path
create_forest_plot_with_counts(
  rate_to_plot = "TERR->FS",
  output_file = "my_custom_forest_plot.pdf"
)

# Override reference value if needed
create_forest_plot_with_counts(
  rate_to_plot = "COOP->FS",
  reference_value = 0.6,
  full_dataset_n = 900
)
```

#### Function Parameters:
- `bias_results_list`: List of bias correction results (default: NULL, auto-loads)
- `rate_to_plot`: Which path to plot (default: "COOP->FS")
- `reference_value`: Full dataset coefficient (default: NULL, auto-loads from RDS)
- `full_dataset_n`: Full dataset sample size (default: NULL, auto-loads from RDS)
- `output_file`: Output path (default: NULL, auto-generates based on rate)

#### Visual Design:
- Darker blue points (#1B4F72) for contrast with reference line (#2E86AB)
- Subtitle clarifies "Mean and 95% CI of conditional averages from 500 downsampling iterations"
- Bias corrections ordered: Plumage, Wing, Territoriality (Level 3), Territoriality (Strong), Geographic (Global), Geographic (Tropical), Geographic (Holarctic)
- Creates both PDF (600 dpi) and PNG (300 dpi) versions

#### Dependencies:
- Requires `extract_phylopath_results()` from `create_phylopath_bias_robustness_figure.R`
- Requires `phylopath_full_dataset_result.rds` for automatic reference values
- Uses libraries: ggplot2, dplyr, phylopath

## Recent Progress (Update 10 - Jackknife Species Analysis)

### New Jackknife Species Analysis Function: `run_jackknife_species_phylopath.R`

A new function has been created to perform jackknife analysis by iteratively removing each species one at a time. This analysis helps identify influential species that significantly affect the phylopath model results.

#### Key Features:
1. **Iterative species removal** - Removes each of 875 species one at a time
2. **Tracks removed species** - Adds "RemovedSpecies" column to all outputs
3. **Maintains seed column** - Uses iteration number as seed for reproducibility
4. **Comprehensive outputs** - Saves detailed models, summaries, and edge statistics
5. **Progress tracking** - Shows elapsed time and estimated completion
6. **Error handling** - Continues analysis even if individual iterations fail

#### Output Files:
All files include "JackknifeSpecies" keyword and are saved to `Outputs/PhylopathJackknife/`:
- `results_JackknifeSpecies_n875_[date].csv` - Main results with best models
- `detailed_models_JackknifeSpecies_n875_[date].csv` - All models with ΔCICc < 2
- `model_frequencies_JackknifeSpecies_n875_[date].csv` - Model selection frequency
- `edge_summary_JackknifeSpecies_n875_[date].csv` - Path coefficient statistics
- `conditional_average_plots_JackknifeSpecies_n875_[date].pdf` - Optional DAG plots

#### Usage Example:
```r
# Source required functions
source("run_phylopath_fxns.R")
source("run_jackknife_species_phylopath.R")

# Load data
df <- read.csv("Data_R_2025-06-09.csv", stringsAsFactors = FALSE)
tree <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Run jackknife analysis (will take several hours)
jackknife_results <- run_jackknife_species_phylopath(
  dfIn = df,
  tree = tree,
  female_song_var = "FemaleSong_Agg01",
  coop_breeding_var = "HighConfidence_Coop", 
  territoriality_var = "TerritorialityWeakVsStrong",
  mass_var = "logMass_AVONET",
  save_conditional_plots = FALSE,  # Set TRUE for DAG plots
  save_path_coefficients = TRUE
)

# Analyze results to find influential species
influential_species <- jackknife_results$results_df %>%
  filter(best_model != mode(best_model)) %>%  # Different from most common
  select(RemovedSpecies, best_model, nSub2dCIC_Models)
```

#### Analysis Applications:
- Identify species that strongly influence model selection
- Test robustness of CB→FS association to individual species
- Find potential outliers or data quality issues
- Assess phylogenetic influence of particular clades

#### Computational Considerations:
- Full analysis takes approximately 3-4 hours for 875 species
- Progress updates every 10 iterations
- Can be interrupted and restarted (results saved incrementally)
- Memory efficient - processes one iteration at a time

## Notes
- User wants plots to be INFORMATIVE, COMPLETE, NEAT, and CLEAR
- Generate individual plots first before combining
- Consider the biological interpretation in all visualizations
- Panel B in multi-panel figures should NOT compare CICc values across different datasets
- DAG should show bidirectional connections with curved arrows and both rate labels visible
- Model consistency heatmap function moved to standalone script for better organization
- Forest plot function also moved to standalone script for modularity
- Jackknife analysis provides species-level sensitivity testing