# Box_Scripts vs Main Directory Comparison Report

## Overview
This report compares three R scripts found in both the `Box_Scripts/` folder and the main directory to identify version differences and determine which should be used.

## 1. subsettreedata.R

### File Information
- **Main directory version**: 9,016 bytes (Last modified: 9-14-2023)
- **Box_Scripts version**: 6,469 bytes (Last modified: 6-18-2024)

### Key Differences
1. **Package versions**:
   - Main: `ape_5.3  phytools_0.5-38   maps_3.1.0  btw_V1.0`
   - Box_Scripts: `ape_5.7-1  phytools_1.9-16`

2. **Default data and tree files**:
   - Main: Uses `SnyderCreanza_NatComms2019_SupplementalData_R.csv` and `birdzillatreeMaybeConsensus.nex`
   - Box_Scripts: Uses `Data_R.csv` and `ConsensusPasserineTreeHackett4_1000_OscineSubset.nex`

3. **Code cleanup**:
   - Box_Scripts version has removed extensive comments about update history
   - Box_Scripts version has removed commented-out code for directory creation
   - Box_Scripts version is cleaner and more streamlined

### Recommendation
**Use the Box_Scripts version** - It's more recent (2024 vs 2023), has updated package dependencies, uses newer data files, and has cleaner code without sacrificing functionality.

## 2. findQrates.R

### File Information
- **Main directory version**: 6,313 bytes (Last modified: 5-4-2024)
- **Box_Scripts version**: 6,313 bytes (Last modified: 5-4-2024)

### Key Differences
**None** - The files are identical (verified with diff command).

### Recommendation
**Either version can be used** - They are the same file.

## 3. test_trait_overlap_simmaps.R

### File Information
- **Main directory version**: 32,288 bytes (Last modified: 6-2-2025)
- **Box_Scripts version**: 33,069 bytes (Last modified: 2-3-2025)

### Key Differences
1. **Recent updates in Box_Scripts version**:
   - Added `setQratesData` parameter to `CharacterSimmaps` function (2/3/2025)
   - Modified `calcHuelflex` function to dynamically find numeric columns (12/11/2024)
   - More robust handling of column detection using `numericColumnStart`

2. **Functional improvements**:
   - Box_Scripts version allows setting both tree and data for Q rate calculations
   - Better handling of different data formats in `calcHuelflex`

### Recommendation
**Use the Box_Scripts version** - It has more recent updates (February 2025) with important functional improvements, including better flexibility for Q rate calculations and more robust column handling.

## Summary

1. **subsettreedata.R**: Use Box_Scripts version (newer, cleaner)
2. **findQrates.R**: Use either version (identical)
3. **test_trait_overlap_simmaps.R**: Use Box_Scripts version (more recent with functional improvements)

## Additional Files in Box_Scripts
The Box_Scripts folder also contains:
- `Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv` (data file)
- `TransitionCounts_3trait_flexTerr_fxns.R`
- `bias tests.R`
- `run_phylopath_fxns.R`

These appear to be additional analysis scripts not present in the main directory.