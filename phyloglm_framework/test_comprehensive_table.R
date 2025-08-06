# Test the comprehensive table creation with just one analysis

source("phyloglm_framework/stepwise_expansion_comprehensive_table.R")

# Run with limited scope to test
results <- run_stepwise_expansion_comprehensive_table()

# Check if output was created
output_dir <- paste0("Outputs/PhyloglmResults/stepwise_comprehensive_table_", format(Sys.Date(), "%Y%m%d"))
if (file.exists(file.path(output_dir, "comprehensive_improvement_comparison_table.csv"))) {
  cat("\nTable created successfully!\n")
  
  # Read and show first few rows
  table_data <- read.csv(file.path(output_dir, "comprehensive_improvement_comparison_table.csv"))
  cat("\nFirst few rows of comprehensive table:\n")
  print(head(table_data[, c("Analysis", "Predictor", "N_Species", "AIC_Improvement", "Predictor_p_value", "Note")]))
  
  cat("\nColumn names in table:\n")
  print(names(table_data))
} else {
  cat("\nTable was not created. Check for errors above.\n")
}