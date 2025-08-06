# Script to view stepwise expansion results

# Load results
results <- readRDS("Outputs/PhyloglmResults/stepwise_20250715/stepwise_results.rds")

cat("STEPWISE MODEL EXPANSION RESULTS\n")
cat("================================\n\n")

# For each analysis
for (name in names(results)) {
  cat(name, "\n")
  cat(paste(rep("-", nchar(name)), collapse=""), "\n")
  
  result <- results[[name]]
  
  # Show improvements
  cat("\nPredictors tested:\n")
  print(result$improvements)
  
  # If we have a final model, show the coefficients
  if ("final_effects" %in% names(result)) {
    cat("\nFinal model effects:\n")
    effects <- result$final_effects
    effects$OR_text <- paste0(round(effects$OddsRatio, 2), 
                             " (", round(effects$OR_CI_lower, 2), 
                             "-", round(effects$OR_CI_upper, 2), ")")
    effects$Sig <- ifelse(effects$P.Value < 0.001, "***",
                         ifelse(effects$P.Value < 0.01, "**",
                               ifelse(effects$P.Value < 0.05, "*", "")))
    
    print(effects[, c("Parameter", "OR_text", "P.Value", "Sig")])
  }
  
  cat("\n\n")
}

# Show the images
cat("Visualizations saved:\n")
cat("- Outputs/PhyloglmResults/stepwise_20250715/stepwise_improvements.png\n")
cat("- Outputs/PhyloglmResults/stepwise_20250715/expanded_model_effects.png\n")