# Verify that our conditional averaging matches phylopath::average()
library(phylopath)
library(dplyr)

# First, let's check if we have saved phylopath results to compare with
source('create_phylopath_bias_robustness_figure.R')

# Get one example dataset
results <- extract_phylopath_results()
result <- results[['Geographic (Holarctic)']]

# Look for a specific seed to test
test_seed <- 1001
seed_data <- result$detailed_models %>%
  filter(seed == test_seed)

cat("Checking if we have phylopath result objects saved...\n")

# Check if there are RDS files with phylopath results
rds_files <- list.files("Outputs/PhylopathDownsampled", 
                       pattern = ".*HolarcticNoncoop.*\\.rds$", 
                       full.names = TRUE, 
                       recursive = TRUE)

if (length(rds_files) > 0) {
  cat(paste("\nFound", length(rds_files), "RDS files\n"))
  cat("First few files:\n")
  print(head(rds_files))
  
  # Try to load one and check its structure
  test_result <- readRDS(rds_files[1])
  cat("\nStructure of saved result:\n")
  print(names(test_result))
  
  if ("result" %in% names(test_result)) {
    cat("\nChecking if this is a phylopath result object...\n")
    print(class(test_result$result))
    
    # If it's a phylopath result, we can use average() on it
    if (inherits(test_result$result, "phylopath")) {
      cat("\nThis is a phylopath object! Testing average() function...\n")
      
      # Get conditional average using phylopath
      avg_result <- average(test_result$result, avg_method = "conditional")
      
      # Extract the CB->FS coefficient
      if (!is.null(avg_result$coef)) {
        cb_fs_coef <- avg_result$coef["HighConfidence_Coop", "FemaleSong_Agg01"]
        cat(paste("\nPhylopath conditional average for CB->FS:", round(cb_fs_coef, 4), "\n"))
      }
    }
  }
} else {
  cat("\nNo RDS files found. Let me check the CSV structure more carefully...\n")
  
  # Check if we have all the necessary information to replicate phylopath averaging
  cat("\nColumns available in detailed_models:\n")
  print(names(seed_data))
  
  # For conditional averaging, we need:
  # 1. Models with delta_CICc < 2
  # 2. CICc weights (or we calculate them)
  # 3. Path coefficients
  
  models_under_2 <- seed_data %>%
    filter(delta_CICc < 2)
  
  cat(paste("\nModels with delta_CICc < 2 for seed", test_seed, ":", nrow(models_under_2), "\n"))
  
  # Show the calculation step by step
  coef_col <- 'HighConfidence_Coop_to_FemaleSong_Agg01_est'
  
  cat("\nStep-by-step conditional averaging:\n")
  cat("1. Models and their delta_CICc values:\n")
  print(models_under_2[, c("model", "delta_CICc", coef_col)])
  
  cat("\n2. Calculate weights = exp(-0.5 * delta_CICc):\n")
  models_under_2$weight <- exp(-0.5 * models_under_2$delta_CICc)
  print(models_under_2[, c("model", "delta_CICc", "weight")])
  
  cat("\n3. Normalize weights to sum to 1:\n")
  models_under_2$norm_weight <- models_under_2$weight / sum(models_under_2$weight)
  print(models_under_2[, c("model", "weight", "norm_weight")])
  
  # Only include models that have the coefficient
  valid_models <- models_under_2[!is.na(models_under_2[[coef_col]]), ]
  
  if (nrow(valid_models) > 0) {
    cat("\n4. Models with CB->FS coefficient:\n")
    print(valid_models[, c("model", coef_col, "norm_weight")])
    
    # Need to renormalize weights for just the models with coefficients
    valid_models$final_weight <- valid_models$norm_weight / sum(valid_models$norm_weight)
    
    cat("\n5. Renormalized weights for models with coefficient:\n")
    print(valid_models[, c("model", "norm_weight", "final_weight")])
    
    # Calculate conditional average
    cond_avg <- sum(valid_models[[coef_col]] * valid_models$final_weight)
    cat(paste("\n6. Conditional average =", round(cond_avg, 4), "\n"))
  }
}