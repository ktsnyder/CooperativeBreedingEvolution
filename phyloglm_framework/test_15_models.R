# Test script to verify that all 15 models are generated correctly

# Source framework
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/config_builder.R")
source("phyloglm_framework/variable_classification.R")

# Create a test configuration with 2 predictors + 1 control
test_config <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("logMass_AVONET"),
  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
  max_interactions = 3,
  name = "Test_15_Models"
)

# Load dummy data for testing
data <- data.frame(
  species = paste0("Species", 1:100),
  FemaleSong_Agg01 = rbinom(100, 1, 0.5),
  HighConfidence_Coop = rbinom(100, 1, 0.3),
  TerritorialityWeakVsStrong = rbinom(100, 1, 0.6),
  logMass_AVONET = rnorm(100, 3, 0.5)
)

# Build model set
formulas <- build_model_set_from_config(test_config, data)

# Expected 15 model names based on reference
expected_names <- c(
  "Null",
  "MainPred",
  "Terr", 
  "Mass",
  "MainPred_Terr",
  "MainPred_Mass",
  "Terr_Mass",
  "MainPredxTerr",
  "MainPredxMass",
  "TerrxMass",
  "MainPred_Terr_Mass",
  "MainPredxTerr_Mass",
  "MainPred_TerrxMass",
  "Terr_MainPredxMass",
  "FullInteractions"
)

# Check results
cat("=== Formula Generation Test ===\n")
cat("Number of formulas generated:", length(formulas), "\n")
cat("Expected:", length(expected_names), "\n\n")

# Print all formulas
cat("Generated formulas:\n")
for (i in seq_along(formulas)) {
  cat(sprintf("%2d. %-20s: %s\n", 
              i, 
              names(formulas)[i], 
              as.character(formulas[[i]])[3]))
}

# Check if all expected names are present
cat("\n=== Name Matching ===\n")
missing_names <- setdiff(expected_names, names(formulas))
if (length(missing_names) == 0) {
  cat("✓ All expected model names are present!\n")
} else {
  cat("✗ Missing model names:", paste(missing_names, collapse = ", "), "\n")
}

extra_names <- setdiff(names(formulas), expected_names)
if (length(extra_names) > 0) {
  cat("✗ Extra model names:", paste(extra_names, collapse = ", "), "\n")
}

# Verify exact formula structures for key models
cat("\n=== Formula Structure Verification ===\n")

check_formula <- function(name, expected_rhs) {
  if (name %in% names(formulas)) {
    actual_rhs <- as.character(formulas[[name]])[3]
    if (actual_rhs == expected_rhs) {
      cat(sprintf("✓ %-20s: %s\n", name, actual_rhs))
    } else {
      cat(sprintf("✗ %-20s: Expected '%s', got '%s'\n", name, expected_rhs, actual_rhs))
    }
  } else {
    cat(sprintf("✗ %-20s: Not found\n", name))
  }
}

# Check key formulas
check_formula("Null", "1")
check_formula("MainPred", "HighConfidence_Coop")
check_formula("MainPred_Terr", "HighConfidence_Coop + TerritorialityWeakVsStrong")
check_formula("MainPredxTerr", "HighConfidence_Coop * TerritorialityWeakVsStrong")
check_formula("FullInteractions", "HighConfidence_Coop * TerritorialityWeakVsStrong * logMass_AVONET")

# Test with different control variables
cat("\n\n=== Testing with Different Control Variables ===\n")

# Test with wing dimorphism
test_config2 <- create_analysis_config(
  response = "FemaleSong_Agg01",
  predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
  controls = c("PercentAbsLogWingDimorphism"),
  name = "Test_WingDim"
)

# Add wing dimorphism to data
data$PercentAbsLogWingDimorphism <- rnorm(100, 0, 1)

formulas2 <- build_model_set_from_config(test_config2, data)
cat("\nWing dimorphism config generated", length(formulas2), "formulas\n")

if (length(formulas2) == 15) {
  cat("✓ Correct number of formulas for wing dimorphism analysis\n")
} else {
  cat("✗ Incorrect number of formulas for wing dimorphism analysis\n")
}