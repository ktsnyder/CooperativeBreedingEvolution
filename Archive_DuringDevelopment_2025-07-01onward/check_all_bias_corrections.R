# Check coefficient variation across all bias corrections
source('create_phylopath_bias_robustness_figure.R')
library(dplyr)

results <- extract_phylopath_results()

cat("Checking coefficient variation across all bias corrections:\n\n")

for (bias_name in names(results)) {
  result <- results[[bias_name]]
  
  if (!is.null(result$detailed_models)) {
    # Find the COOP->FS coefficient column
    coef_col <- grep("Coop.*to.*FemaleSong.*est", names(result$detailed_models), 
                     value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.na(coef_col)) {
      # Count seeds with varying coefficients
      seed_summary <- result$detailed_models %>%
        filter(delta_CICc < 2 & !is.na(.data[[coef_col]])) %>%
        group_by(seed) %>%
        summarise(
          n_unique_coefs = n_distinct(.data[[coef_col]]),
          n_models = n(),
          .groups = "drop"
        )
      
      seeds_with_variation <- sum(seed_summary$n_unique_coefs > 1)
      total_seeds <- nrow(seed_summary)
      
      cat(paste(bias_name, ":\n"))
      cat(paste("  Seeds with coefficient variation:", seeds_with_variation, "out of", total_seeds, 
                "(", round(100 * seeds_with_variation / total_seeds, 1), "%)\n"))
      
      # Calculate the average difference between methods
      if (seeds_with_variation > 0) {
        # Just for seeds with variation, calculate the difference
        varied_seeds <- seed_summary$seed[seed_summary$n_unique_coefs > 1]
        
        differences <- numeric()
        for (s in varied_seeds[1:min(10, length(varied_seeds))]) {
          seed_data <- result$detailed_models %>%
            filter(seed == s & delta_CICc < 2 & !is.na(.data[[coef_col]]))
          
          # Simple mean
          simple_mean <- mean(seed_data[[coef_col]])
          
          # Conditional average
          all_weights <- exp(-0.5 * seed_data$delta_CICc)
          all_weights <- all_weights / sum(all_weights)
          cond_avg <- sum(seed_data[[coef_col]] * all_weights)
          
          differences <- c(differences, abs(cond_avg - simple_mean))
        }
        
        cat(paste("  Average difference when coefficients vary:", 
                  round(mean(differences), 6), "\n"))
      }
      cat("\n")
    }
  }
}