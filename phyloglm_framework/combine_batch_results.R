# Script to combine batch results for plotting
# Combines batches 1+2 (Female Song as response) and batches 3+4 (Cooperative Breeding as response)

library(dplyr)

# Function to combine two batch results
combine_batch_results <- function(batch1_path, batch2_path, output_name) {
  
  cat("Loading batch results...\n")
  # Load the all_results.rds files
  batch1 <- readRDS(batch1_path)
  batch2 <- readRDS(batch2_path)
  
  # Combine the results lists
  combined_results <- c(batch1, batch2)
  
  cat("Combined", length(batch1), "and", length(batch2), "analyses into", 
      length(combined_results), "total analyses\n")
  
  # Create output directory
  output_dir <- file.path("Outputs/PhyloglmResults", output_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save combined results
  saveRDS(combined_results, file.path(output_dir, "all_results.rds"))
  
  # Also combine the CSV files for convenience
  combine_csv_files <- function(file1_dir, file2_dir, output_dir, filename) {
    file1 <- file.path(file1_dir, "integrated_results", filename)
    file2 <- file.path(file2_dir, "integrated_results", filename)
    
    if (file.exists(file1) && file.exists(file2)) {
      df1 <- read.csv(file1)
      df2 <- read.csv(file2)
      combined <- rbind(df1, df2)
      
      # Create integrated_results directory
      int_dir <- file.path(output_dir, "integrated_results")
      dir.create(int_dir, showWarnings = FALSE)
      
      write.csv(combined, file.path(int_dir, filename), row.names = FALSE)
      cat("  Combined", filename, "-", nrow(combined), "rows\n")
    }
  }
  
  # Get batch directories
  batch1_dir <- dirname(batch1_path)
  batch2_dir <- dirname(batch2_path)
  
  # Combine each type of CSV file
  cat("\nCombining CSV files...\n")
  combine_csv_files(batch1_dir, batch2_dir, output_dir, "all_coefficients.csv")
  combine_csv_files(batch1_dir, batch2_dir, output_dir, "all_model_comparisons.csv")
  combine_csv_files(batch1_dir, batch2_dir, output_dir, "effect_comparison.csv")
  combine_csv_files(batch1_dir, batch2_dir, output_dir, "analysis_summary_stats.csv")
  
  # Create a summary of what's in the combined results
  cat("\nCreating analysis summary...\n")
  summary_info <- data.frame(
    Analysis = names(combined_results),
    Response = sapply(combined_results, function(x) x$config$response),
    N_Models = sapply(combined_results, function(x) 
      if (!is.null(x$comparison)) nrow(x$comparison$comparison) else 0),
    Success = sapply(combined_results, function(x) x$success),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_info, file.path(output_dir, "combined_analysis_summary.csv"), 
            row.names = FALSE)
  
  cat("\nSummary of combined results:\n")
  print(table(summary_info$Response))
  cat("\nSuccessful analyses:", sum(summary_info$Success), "out of", nrow(summary_info), "\n")
  
  return(combined_results)
}

# Combine batches 1 & 2 (Female Song as response)
cat("=== Combining Batches 1 & 2 (Female Song as Response) ===\n")
fs_results <- combine_batch_results(
  batch1_path = "Outputs/PhyloglmResults/batch1/PhyloGLM_Batch_20250711_034017_boot1000/all_results.rds",
  batch2_path = "Outputs/PhyloglmResults/batch2/PhyloGLM_Batch_20250711_034718_boot1000/all_results.rds",
  output_name = "combined_FS_response"
)

cat("\n")

# Combine batches 3 & 4 (Cooperative Breeding as response)
cat("=== Combining Batches 3 & 4 (Cooperative Breeding as Response) ===\n")
cb_results <- combine_batch_results(
  batch1_path = "Outputs/PhyloglmResults/batch3/PhyloGLM_Batch_20250711_035451_boot1000/all_results.rds",
  batch2_path = "Outputs/PhyloglmResults/batch4/PhyloGLM_Batch_20250711_040021_boot1000/all_results.rds",
  output_name = "combined_CB_response"
)

cat("\n")

# Combine all batches for complete analysis
cat("=== Combining All Batches (Complete Dataset) ===\n")
all_results <- combine_batch_results(
  batch1_path = "Outputs/PhyloglmResults/combined_FS_response/all_results.rds",
  batch2_path = "Outputs/PhyloglmResults/combined_CB_response/all_results.rds",
  output_name = "combined_all_analyses"
)

cat("\n=== Complete! ===\n")
cat("\nCombined results saved to:\n")
cat("- Female Song analyses: Outputs/PhyloglmResults/combined_FS_response/\n")
cat("- Cooperative Breeding analyses: Outputs/PhyloglmResults/combined_CB_response/\n")
cat("- All analyses: Outputs/PhyloglmResults/combined_all_analyses/\n")

cat("\nYou can now create plots using:\n")
cat("source('phyloglm_framework/phyloglm_summary_plots_simple.R')\n")
cat("# For all analyses:\n")
cat("all_results <- readRDS('Outputs/PhyloglmResults/combined_all_analyses/all_results.rds')\n")
cat("create_simple_summary_figure(all_results, output_file = 'phyloglm_results_all_1000boot.png')\n")
cat("\n# For Female Song analyses only:\n")
cat("fs_results <- readRDS('Outputs/PhyloglmResults/combined_FS_response/all_results.rds')\n")
cat("# Create custom plots for FS results...\n")
cat("\n# For Cooperative Breeding analyses only:\n")
cat("cb_results <- readRDS('Outputs/PhyloglmResults/combined_CB_response/all_results.rds')\n")
cat("# Create custom plots for CB results...\n")