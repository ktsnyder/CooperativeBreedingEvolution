# Box_Scripts vs Main Directory Comparison Report - Updated
**Date:** 2025-06-01
**Purpose:** Document all Box scripts and their dependencies

## Core Script Comparisons

### 1. subsettreedata.R
- **Main directory version**: 9,016 bytes (Last modified: 9-14-2023)
- **Box_Scripts version**: 6,469 bytes (Last modified: 6-18-2024)
- **Recommendation**: ✅ Use Box_Scripts version (newer, cleaner, updated dependencies)

### 2. findQrates.R
- **Main directory version**: 6,313 bytes (Last modified: 5-4-2024)
- **Box_Scripts version**: 6,313 bytes (Last modified: 5-4-2024)
- **Recommendation**: ✅ Either version (identical files)

### 3. test_trait_overlap_simmaps.R
- **Main directory version**: 32,288 bytes (Last modified: 6-2-2025)
- **Box_Scripts version**: 33,069 bytes (Last modified: 2-3-2025)
- **Recommendation**: ✅ Use Box_Scripts version (recent functional improvements)

### 4. find transition counts by state for 2 Discrete traits.R
- **Main directory version**: Present
- **Box_Scripts version**: Present
- **Status**: Identical files (already in both locations)

### 5. transition_plot.R
- **Main directory version**: Present
- **Box_Scripts version**: Present
- **Status**: Identical files (already in both locations)

## New Box Scripts Added

### Priority 1: Core Infrastructure
1. **bias tests.R** ✅
   - Geographic/territorial bias testing framework
   - Tests for Holarctic vs Tropical biases
   - Integrates with downsampling strategies

2. **run_phylopath_fxns.R** ✅
   - Contains `run_multiple_phylopath()` - Runs 500+ iterations with aggregation
   - Contains `run_CB_FS_Terr_phylopath()` - Core phylopath analysis
   - Contains `downsample_run_phylopath()` - Bias correction through downsampling
   - Helper function `build_formula()` for model construction

3. **execute_run_phylopath_fxns.R** ✅
   - Script to execute the phylopath analyses
   - Sets up parameters and runs the iterative analyses

### Priority 2: 3-Trait Analysis Scripts
4. **TransitionCounts_3trait_flexTerr_fxns.R** ✅
   - Main 3-trait analysis functions
   - Contains `plot_transition_counts_3trait()` function
   - Includes `getLabels()` function (resolves dependency)
   - Flexible territoriality handling

5. **TransitionCounts_3trait.R** ✅
   - Original 3-trait transition counting
   - Contains another version of `getLabels()`

6. **TransitionCounts_3trait_flexTerr.R** ✅
   - Flexible territoriality version
   - Contains `getLabels()` function

7. **TransitionCounts_3trait_FlexClaude.R** ✅
   - Claude-modified version for flexibility
   - Contains `getLabels()` function

8. **MapOverlapThree.R** ✅
   - 3-trait overlap mapping functionality
   - Works with transition count scripts

9. **run simmap overlap 3traits.R** ✅
   - Executes 3-trait simmap overlap analyses
   - Coordinates multiple trait analyses

### Priority 3: Supporting Scripts
10. **scratch testing data biases.R** ✅
    - Testing and exploration of bias patterns
    - Development/testing script

11. **summarize data columns.R** ✅
    - Data exploration and summary functions
    - Helps understand data structure

### Data File
12. **Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv** ✅
    - Updated dataset with Tobias territoriality classifications
    - Includes weak/strong territorial distinctions
    - Essential for bias-corrected analyses

### Master Analysis Script
13. **Run_Analyses.R** ✅
    - Master script for running all manuscript analyses
    - Box version (2/3/2025) vs Sandbox version (6/19/2024)
    - Box version includes:
      - New Cockburn-inferred cooperative breeding classifications
      - Tobias territory integration (added 12/9/2024)
      - Territory-specific simmap overlap tests (added 2/3/2025)
      - More comprehensive analysis pipeline
    - **Recommendation**: Use Box version (more recent with additional analyses)

## Dependency Resolution

### Previously Missing Dependencies - NOW RESOLVED:
1. ✅ `getLabels()` function - Found in multiple TransitionCounts scripts
2. ✅ `downsample_run_phylopath()` function - Found in run_phylopath_fxns.R
3. ✅ `build_formula()` helper function - Found in run_phylopath_fxns.R
4. ✅ Tobias territoriality data - Included in the CSV file

### Key Dependencies Between Scripts:
```
bias tests.R 
    ↓
run_phylopath_fxns.R (uses bias testing results)
    ↓
execute_run_phylopath_fxns.R (runs the analyses)

TransitionCounts_3trait_flexTerr_fxns.R
    ├─ uses findQrates.R ✅
    ├─ uses subsettreedata.R ✅
    └─ works with MapOverlapThree.R ✅

test_trait_overlap_simmaps.R
    ├─ uses findQrates.R ✅
    ├─ uses subsettreedata.R ✅
    └─ uses find transition counts by state for 2 Discrete traits.R ✅
```

## Integration Status

### Ready for Integration:
1. **Bias correction pipeline**: All components present
   - bias tests.R → run_phylopath_fxns.R → execute_run_phylopath_fxns.R

2. **3-trait analysis pipeline**: All components present
   - TransitionCounts_3trait variants + MapOverlapThree.R + supporting functions

3. **Core utilities**: Updated versions ready
   - subsettreedata.R (Box version)
   - test_trait_overlap_simmaps.R (Box version)

### Path Updates Needed:
- `bias tests.R` contains hardcoded Box paths that need updating:
  - Line 16: Box cloud storage path
  - Line 22: Another Box cloud storage path

## Recommendations

1. **Immediate Actions**:
   - Update hardcoded paths in `bias tests.R`
   - Use Box versions of `subsettreedata.R` and `test_trait_overlap_simmaps.R`
   - Test integration of bias correction pipeline

2. **Version Management**:
   - Archive current main directory versions before replacement
   - Document which TransitionCounts_3trait variant to use as primary

3. **Next Steps**:
   - Run test analyses with small datasets to verify integration
   - Check for any additional Box scripts that might enhance browniefunction.R
   - Update Run_Analyses.R to incorporate new bias correction pipeline