# Test comprehensive pipeline with just one analysis

source("phyloglm_framework/stepwise_expansion_comprehensive.R")

# Run with just one analysis to test
results <- run_comprehensive_stepwise_expansion()

# Check what we got
if (exists("results")) {
  cat("\nAnalyses completed:\n")
  print(names(results))
}