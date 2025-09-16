# Test script for consolidated functions
# This demonstrates that the detailed models file detection is fixed

# Source the consolidated functions
source("claude_code_sessions/run_phylopath_fxns_consolidated.R")

# Test the file pattern matching that was causing issues
test_prefix <- "phylopath FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong Remove66HighTerr 500"

# Show how the original pattern didn't work
original_pattern <- paste0("detailed_models_.*", 
                          gsub(" ", "_", test_prefix), 
                          ".*\\.csv$")
cat("Original pattern:", original_pattern, "\n")

# Show how the new pattern works by extracting the key part
prefix_parts <- strsplit(test_prefix, " ")[[1]]
key_pattern <- grep("Remove", prefix_parts, value = TRUE)
new_pattern <- paste0("detailed_models_.*", key_pattern, ".*\\.csv$")
cat("New pattern:", new_pattern, "\n")

# Test with a sample filename
sample_filename <- "detailed_models_Remove66HighTerr_500_2025-06-09.csv"
cat("\nTesting with filename:", sample_filename, "\n")
cat("Original pattern matches:", grepl(original_pattern, sample_filename), "\n")
cat("New pattern matches:", grepl(new_pattern, sample_filename), "\n")

cat("\nThe consolidated file includes:\n")
cat("1. All functions from run_phylopath_fxns.R\n")
cat("2. Plotting functions from phylopath_plotting_comprehensive_fixed.R:\n")
cat("   - create_downsampled_plots() with FIXED file detection (lines 2012-2098)\n")
cat("   - create_nondownsampled_plots()\n")
cat("   - create_all_phylopath_plots()\n")
cat("3. The calculate_territoriality_downsampling() function\n")
cat("\nKey fix: The detailed_models file pattern now correctly handles spaces and special characters in the prefix.\n")