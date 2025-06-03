# Script Comparison Report: Main Directory vs Sandbox
**Date:** January 6, 2025  
**Project:** CooperativeBreedingEvolution

## Executive Summary

This report compares R scripts that exist in both the main directory (`/Users/kate/Desktop/CooperativeBreedingEvolution/`) and the Sandbox directory (`/Users/kate/Desktop/CooperativeBreedingEvolution/Sandbox/`). The analysis identifies which versions should be used going forward based on file modifications, code differences, and functionality.

## Overview of Findings

### Total Scripts Analyzed
- **Main directory:** 54 R scripts
- **Sandbox directory:** 11 R scripts
- **Scripts in both locations:** 9 scripts
- **Unique to Sandbox:** 2 scripts (`BrownieMultistate.R`, `Run_Analyses.R`)

## Detailed Comparison of Duplicate Scripts

### 1. **plotbrownie.R**
- **File size:** 6,192 bytes (both versions)
- **Last modified:** March 7, 2022 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 2. **jackknifingbrownie.R**
- **File size:** 7,516 bytes (both versions)
- **Last modified:** November 15, 2023 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 3. **find transition counts by state for 2 Discrete traits.R**
- **File size:** 7,506 bytes (both versions)
- **Last modified:** November 27, 2023 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 4. **findQrates.R**
- **File size:** 6,313 bytes (both versions)
- **Last modified:** May 4, 2024 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 5. **transition_plot.R**
- **File size:** 14,338 bytes (both versions)
- **Last modified:** May 23, 2024 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 6. **browniefunction.R**
- **File size:** 5,928 bytes (both versions)
- **Last modified:** May 24, 2024 (both versions)
- **Status:** IDENTICAL
- **Recommendation:** Use either version (they are the same)

### 7. **subsettreedata.R**
- **File size:** Main: 9,016 bytes | Sandbox: 6,469 bytes
- **Last modified:** Main: September 14, 2023 | Sandbox: June 18, 2024
- **Status:** SIGNIFICANTLY DIFFERENT
- **Key differences:**
  - Sandbox version is more recent (June 2024 vs September 2023)
  - Sandbox version is significantly shorter (removed extensive comments and version history)
  - Sandbox version has updated package versions (ape_5.7-1, phytools_1.9-16)
  - Sandbox version uses different default files:
    - Data: `Data_R.csv` instead of `SnyderCreanza_NatComms2019_SupplementalData_R.csv`
    - Tree: `ConsensusPasserineTreeHackett4_1000_OscineSubset.nex` instead of `birdzillatreeMaybeConsensus.nex`
  - Core functionality appears the same but Sandbox version is cleaner
- **Recommendation:** ⚠️ **Use MAIN version** - Although Sandbox is newer, it appears to be a simplified version for specific analyses. The main version has more complete documentation and flexibility.

### 8. **brownie relative rates.R**
- **File size:** Main: 13,211 bytes | Sandbox: 3,289 bytes
- **Last modified:** June 18, 2024 (both versions, but different times)
- **Status:** SIGNIFICANTLY DIFFERENT
- **Key differences:**
  - Main version is a complete script with full implementation
  - Sandbox version is converted into a function `BrownieRelativeRates()`
  - Sandbox version removes hardcoded file paths and trait lists
  - Sandbox version is more modular and reusable
- **Recommendation:** ✅ **Use SANDBOX version** - It's a cleaner, more modular implementation as a function that can be called with parameters.

### 9. **test_trait_overlap_simmaps.R**
- **File size:** Main: 32,149 bytes | Sandbox: 32,288 bytes
- **Last modified:** Main: June 7, 2024 | Sandbox: June 19, 2024
- **Status:** SLIGHTLY DIFFERENT
- **Key differences:**
  - Sandbox version is 12 days newer
  - Sandbox version has additional functionality (75 more lines)
  - Main version has extensive comments about version history
  - Sandbox version includes additional function `calcHuelflex` for multistate categorical traits
  - Sandbox version has cleaner header (removed detailed change log)
- **Recommendation:** ✅ **Use SANDBOX version** - It's more recent and includes additional functionality for multistate analysis.

## Scripts Unique to Sandbox

### 1. **BrownieMultistate.R**
- **Purpose:** Handles Brownie analysis for multistate discrete traits
- **Status:** Unique to Sandbox, appears to be a specialized function
- **Recommendation:** Keep in Sandbox for multistate analyses

### 2. **Run_Analyses.R**
- **Purpose:** Master script to run all analyses for manuscript generation
- **Last modified:** June 19, 2024
- **Status:** Unique to Sandbox, coordinates analysis workflow
- **Recommendation:** Keep in Sandbox as it appears to be the main analysis runner

## Summary Recommendations

### Scripts to Use from Main Directory:
1. `subsettreedata.R` - More complete with documentation

### Scripts to Use from Sandbox Directory:
1. `brownie relative rates.R` - Better modularized as a function
2. `test_trait_overlap_simmaps.R` - More recent with additional features

### Identical Scripts (use either):
1. `plotbrownie.R`
2. `jackknifingbrownie.R`
3. `find transition counts by state for 2 Discrete traits.R`
4. `findQrates.R`
5. `transition_plot.R`
6. `browniefunction.R`

### Action Items:
1. Consider consolidating the identical scripts to avoid confusion
2. Review whether the simplified `subsettreedata.R` in Sandbox should replace the main version
3. Consider moving the unique Sandbox scripts to main if they're part of the standard workflow
4. Update file paths in scripts as needed when consolidating

## File Organization Recommendation

The Sandbox appears to contain a curated subset of scripts for running specific analyses (particularly for manuscript generation), while the main directory contains the full suite of scripts including various test and scratch files. This separation is logical, but identical files should be consolidated to avoid version control issues.