# Quick test to verify everything is working

source("create_phylopath_bias_robustness_figure_updated.R")

# Test 1: Extract all results for TerritorialityWeakVsStrong
cat("=== Test 1: All results for TerritorialityWeakVsStrong ===\n")
results1 <- extract_phylopath_results(
  trait_set = "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET",
  verbose = TRUE
)

cat("\n\nSummary:\n")
cat("Total bias corrections found:", length(results1), "\n")
for (bias in names(results1)) {
  cat("  -", bias, ":", results1[[bias]]$filename, "\n")
}

# Test 2: Extract only 500-iteration results for Territory_12vs3
cat("\n\n=== Test 2: 500-iteration results for Territory_12vs3 ===\n")
results2 <- extract_phylopath_results(
  trait_set = "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET",
  n_iterations = 500,
  verbose = FALSE  # Less verbose this time
)

cat("\nSummary:\n")
cat("Total bias corrections found:", length(results2), "\n")
for (bias in names(results2)) {
  n_iter <- length(unique(results2[[bias]]$detailed_models$seed))
  cat("  -", bias, ":", n_iter, "iterations\n")
}

# Test 3: Check if full dataset results are found
cat("\n\n=== Test 3: Full dataset results ===\n")
source("phylopath_helper_functions.R")

full1 <- find_or_create_full_dataset_result(
  "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET"
)
if (!is.null(full1)) {
  cat("✓ Found full dataset result for TerritorialityWeakVsStrong\n")
}

full2 <- find_or_create_full_dataset_result(
  "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET"
)
if (!is.null(full2)) {
  cat("✓ Found full dataset result for Territory_12vs3\n")
}

cat("\n=== All tests complete ===\n")