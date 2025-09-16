# Compare the two averaging methods to see why values didn't change
source('create_phylopath_bias_robustness_figure.R')
library(dplyr)

results <- extract_phylopath_results()
result <- results[['Geographic (Holarctic)']]

# Let's compare both methods for multiple seeds
cat("Comparing simple mean vs conditional average for multiple seeds:\n\n")

coef_col <- 'HighConfidence_Coop_to_FemaleSong_Agg01_est'

comparison_results <- data.frame()

for (seed_val in unique(result$detailed_models$seed)[1:20]) {
  seed_data <- result$detailed_models %>%
    filter(seed == seed_val & delta_CICc < 2)
  
  if (nrow(seed_data) > 0) {
    # Method 1: Simple mean (incorrect)
    simple_mean <- mean(seed_data[[coef_col]], na.rm = TRUE)
    
    # Method 2: Conditional average (correct)
    all_weights <- exp(-0.5 * seed_data$delta_CICc)
    all_weights <- all_weights / sum(all_weights)
    
    valid_rows <- !is.na(seed_data[[coef_col]])
    
    if (sum(valid_rows) > 0) {
      valid_weights <- all_weights[valid_rows]
      valid_weights <- valid_weights / sum(valid_weights)
      cond_avg <- sum(seed_data[[coef_col]][valid_rows] * valid_weights)
      
      # Check if all non-NA coefficients have the same value
      unique_coefs <- unique(seed_data[[coef_col]][valid_rows])
      
      comparison_results <- rbind(comparison_results, data.frame(
        seed = seed_val,
        simple_mean = simple_mean,
        cond_avg = cond_avg,
        difference = abs(cond_avg - simple_mean),
        n_models = sum(valid_rows),
        n_unique_coefs = length(unique_coefs),
        all_same = length(unique_coefs) == 1
      ))
    }
  }
}

# Show results
print(comparison_results)

cat("\n\nSummary:\n")
cat(paste("Seeds where all coefficients are the same:", sum(comparison_results$all_same), "out of", nrow(comparison_results), "\n"))
cat(paste("Average difference between methods:", round(mean(comparison_results$difference), 6), "\n"))
cat(paste("Max difference between methods:", round(max(comparison_results$difference), 6), "\n"))

# Let's look at one seed where coefficients might vary
varied_seed <- comparison_results$seed[which.max(comparison_results$n_unique_coefs)]
if (length(varied_seed) > 0) {
  cat(paste("\n\nDetailed look at seed", varied_seed, "with", 
            comparison_results$n_unique_coefs[comparison_results$seed == varied_seed], 
            "unique coefficient values:\n"))
  
  seed_data <- result$detailed_models %>%
    filter(seed == varied_seed & delta_CICc < 2)
  
  valid_rows <- !is.na(seed_data[[coef_col]])
  print(seed_data[valid_rows, c("model", "delta_CICc", coef_col)])
}