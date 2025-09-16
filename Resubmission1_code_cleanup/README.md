# Cooperative Breeding Evolution Analysis

This repository contains the complete computational pipeline used in the manuscript "Territoriality modulates the coevolution of cooperative breeding and female song in songbirds", including phylogenetic analyses, statistical modeling, and data visualization scripts.

## System Requirements

### Software Dependencies
- **R version 4.3.1 or higher** (tested on R 4.3.1 "Beagle Scouts")
- **Operating Systems**: macOS, Windows, Linux
- **RStudio** (recommended for interactive use)

### Required R Packages
The analysis requires 29 R packages from CRAN and Bioconductor:

**CRAN packages:**
- ape, cowplot, dplyr, emmeans, flextable, ggplot2, ggpubr, ggrepel
- grid, gridExtra, mnormt, patchwork, phylolm, phylopath, phytools
- R.utils, stringr, tidyr, tidyverse

**Bioconductor packages:**
- graph, RBGL

*All package dependencies are automatically checked and can be installed using the provided `check_required_packages.R` script.*

### Hardware Requirements
- **RAM**: Minimum 8GB, 16GB+ recommended for full analyses; for demo (using the pre-set values in Runner_Script.R)
- **Storage**: ~5GB free space for outputs and temporary files

## Installation Guide

### 1. Install R and RStudio
- Download R from [CRAN](https://cran.r-project.org/)
- Download RStudio from [Posit](https://posit.co/download/rstudio-desktop/)

### 2. Install Required Packages
Run the automated package checker:
```r
source("check_required_packages.R")
```
When prompted, type `y` and hit Enter/Return to install missing packages.

### 3. Verify Setup
Ensure your working directory contains:
- All R scripts from this repository
- `Source Data Process_CB/` subdirectory
- `Data_R.csv` (main dataset)

**Typical installation time:** 10-30 minutes on a standard desktop computer (depending on internet speed and package compilation requirements).

## Quick Demo

To run a abbreviated version of the analyses for demonstration:

### Instructions
1. Ensure all packages are installed via `check_required_packages.R`
2. Run the demo script:
```r
source("Runner_Script.R")
```

### Expected Output
The demo will generate:
- Console output showing progress and timing for each analysis step
- Multiple output subdirectories (`Outputs/`)
- CSV files with statistical results
- PDF plots and figures  

### Expected Runtime
- **2019 MacBook Pro**: ~60 minutes
- **High-performance Windows desktop**: ~28 minutes
- **Standard desktop computer**: 30-90 minutes (depending on specifications)

*Note: Demo uses reduced simulation parameters. Full manuscript analyses may take several hours to days.*

## Instructions for Use

### Running Full Analyses
1. **Set working directory** to the main project folder containing this README
2. **Modify parameters** in `Runner_Script.R` as needed:
   - `nsim`: Number of simulations (increase for publication-quality results)
   - `nBoot`: Bootstrap iterations (set to 500 for manuscript replication)
   - `n_iterations`: Phylopath iterations (set to 500 for full analyses)
   - `nTreesToSample`: Number of trees for multitree analyses (set to 200 for full)

3. **Execute the pipeline**:
```r
source("Runner_Script.R")
```

### Key Analysis Components

#### Core Analyses (always run):
- **PhylANOVA** (`run_01_phylANOVA.R`): Tests for differences in male song features between non-cooperative and cooperative species
- **Brownie** (`run_02_brownie.R`): Song feature evolution rate analyses  
- **Binary transition rates** (`run_03_binary_transition_rates.R`): ARD vs ER model comparison for binary traits
- **Simmap overlap** (`run_04_simmap_overlap.R`): Testing whether two binary or categorical traits are evolutionarily associated
- **Bias analyses** (`run_05_bias_analyses.R`): Data sampling bias assessment
- **Phylopath** (`run_06_phylopath_main.R`, `run_07_phylopath_downsampling.R`, `run_08_phylopath_altCoop.R`): Causal modeling
- **PhyloGLM** (`run_09_phyloglm.R`): Phylogenetic regression models
- **Multitree tests** (`run_10_multitree_tests.R`): Robustness of analyses across multiple trees
- **Species counts** (`run_11_generate_species_counts_tables.R`): Descriptive statistics

#### Optional Analyses:
- **Consensus trees** (`run_00_make various consensus trees.R`): Generate new consensus phylogenies

### Output Structure
Results are organized in `Outputs/` with subdirectories for each analysis type:
- `Brownie_outputs/`
- `PhyloGLM_outputs/`
- `Phylopath_outputs/`
- `Simmap_outputs/`

## Reproduction Instructions

### Manuscript Results Replication
To reproduce all quantitative results from the manuscript:

1. **Use full parameter settings** in `Runner_Script.R`:
```r
nsim = 1000              # Full simulations
nBoot = 500              # Full bootstrap replicates  
n_iterations = 500       # Full phylopath iterations
nTreesToSample = 200     # Full multitree sampling
jackknife_families_above = 3  # Complete jackknife analysis
```

2. **Enable all analyses**:
```r
run_all_sociality_traits = TRUE    # All trait combinations
run_alt_coop = TRUE                # Alternative cooperative breeding classifications  
count_transitions = TRUE           # Transition count analyses
```

3. **Expected runtime**: 2-7 days depending on computing resources

## File Organization

### Main Scripts
- `Runner_Script.R`: Complete analysis pipeline
- `Runner_Script_timed.R`: Pipeline with timing measurements
- `check_required_packages.R`: Automated dependency management

### Analysis Functions
- `Brownie_functions/`: Character evolution rate analysis functions
- `Phylopath_functions/`: Causal modeling and visualization functions  
- `PhyloGLM_functions/`: Phylogenetic regression utilities
- `Simmap_Overlap_functions/`: Correlated evolution analysis functions
- `Bias_Test_functions/`: Data bias assessment functions

### Data Processing
- `Source Data Process_CB/`: Data compilation and cleaning scripts

## Citation

If you use this code, please cite: Snyder, Loughran-Pierce, and Creanza (2025) bioRxiv.

## Contact

For questions about the code or analyses, please visit the repository page on GitHub.

## AI Acknowledgement

We used Claude Sonnet 4 and Claude Opus 4.1 for occasional coding assistance. All code was reviewed, validated, and implemented by author KTS. 

---

*This analysis pipeline was developed and tested on R 4.3.1. While it may work on earlier R versions, we recommend using R 4.3.1 or later for full compatibility.*