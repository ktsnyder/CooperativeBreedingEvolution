# Folder Organization Plan - Pre-Commit
**Date:** 2025-06-02

## Current State

### 1. Modified Files (5)
- `Generate_Pub_Figs.R` - Added PNG functionality
- `brownie relative rates.R` - Refactored into function
- `multiTip tree plots.R` - Added PNG output
- `plotACEtree.R` - Added PNG output
- `test_trait_overlap_simmaps.R` - Fixed calcHuel PNG generation

### 2. Untracked Files from Box_Scripts Integration (11)
**Already copied to main directory:**
- `bias tests.R`
- `execute_run_phylopath_fxns.R`
- `run_phylopath_fxns.R`
- `MapOverlapThree.R`
- `Run_Analyses.R`
- `TransitionCounts_3trait_flexTerr_fxns.R`

**Documentation files:**
- `Box_Scripts_Comparison_Report.md`
- `Box_Scripts_Comparison_Report_Updated_2025-06-01.md`
- `Box_Scripts_Integration_Status_2025-06-02.md`
- `Script_Comparison_MainDir_vs_Sandbox_2025-06-01.md`

**Data file:**
- `justlongevitydata.csv`

### 3. Box_Scripts Folder (still contains originals)

### 4. claude_code_sessions Folder
- Session documentation files
- Updated Generate_Pub_Figs.R
- Multiple test scripts (7 files)
- Session 003 subfolder with PNG-related fixes

## Recommended Organization Before Commit

### Option 1: Clean Integration (Recommended)
1. **Move test files** from claude_code_sessions to a test subfolder
2. **Archive Box_Scripts** folder after integration is complete
3. **Keep important documentation** in main directory
4. **Create CHANGELOG.md** to document the integration

### Option 2: Preserve Current Structure
1. Keep everything as-is for full traceability
2. Add README in claude_code_sessions explaining contents
3. Keep Box_Scripts as reference

## Proposed Actions

### 1. Create Test Directory
```bash
mkdir -p tests/integration_tests
mv claude_code_sessions/test_*.R tests/integration_tests/
```

### 2. Move Key Files
```bash
# Move the updated Generate_Pub_Figs to replace the original
cp claude_code_sessions/Generate_Pub_Figs_Updated.R Generate_Pub_Figs_BoxIntegrated.R

# Move the helper function to a utilities folder
mkdir -p R/utilities
mv claude_code_sessions/fix_tree_node_labels.R R/utilities/
```

### 3. Archive Box_Scripts
```bash
# Create archive with date stamp
mv Box_Scripts Archive_Box_Scripts_2025-06-02
```

### 4. Consolidate Documentation
```bash
mkdir -p documentation/integration
mv Box_Scripts_*.md documentation/integration/
mv Script_Comparison_*.md documentation/integration/
mv claude_code_sessions/session_*.md documentation/integration/
```

### 5. Create CHANGELOG
Document all changes in a CHANGELOG.md file

## Benefits of Reorganization

1. **Cleaner structure** - Easier to navigate
2. **Clear separation** - Tests, documentation, and code are organized
3. **Preserved history** - Everything is archived, not deleted
4. **Ready for commit** - Professional repository structure

## Alternative: Minimal Changes

If you prefer minimal changes before commit:
1. Just add the untracked files
2. Commit the modified files
3. Leave reorganization for a future commit

This preserves the exact state of the integration work.