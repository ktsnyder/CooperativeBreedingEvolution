# Run stepwise expansion analysis with 500 bootstrap iterations
# Using the updated bootstrap script

source("phyloglm_framework/stepwise_expansion_bootstrap_updated.R")
#source("/Users/kate/Downloads/stepwise_expansion_bootstrap_updated.r") # moved working version to phyloglm_framework

cat("================================================================================\n")
cat("Starting stepwise expansion analysis with 500 bootstrap iterations\n")
cat("This analysis will take several hours to complete.\n")
cat("Bootstrap matrices will be saved separately to manage memory.\n")
cat("Interim results will be saved to prevent data loss.\n")
cat("================================================================================\n\n")


start100time <- Sys.time()
# Run analysis with 500 bootstraps
results <- run_stepwise_expansion_bootstrap(
  n_bootstrap = 100,          # 100 bootstrap iterations
  test_mode = TRUE,          # test mode to start
  save_interim = TRUE,        # Save interim results
  save_boot_matrices = TRUE   # Save bootstrap matrices to disk
)
end100time <- Sys.time()
Total100time = end100time - start100time
print(paste("100boot took", Total100time))

# Run analysis with 500 bootstraps
results <- run_stepwise_expansion_bootstrap(
  n_bootstrap = 500,          # 500 bootstrap iterations
  test_mode = FALSE,          # Not test mode - full analysis
  save_interim = TRUE,        # Save interim results
  save_boot_matrices = TRUE   # Save bootstrap matrices to disk
)
end500time <- Sys.time()
Total500time = end500time - end100time
TotalTime = end500time - start100time
print(paste("500boot took", Total500time))
print(paste("Total time:", TotalTime))

cat("\n================================================================================\n")
cat("Analysis complete!\n")
cat("Results saved to: Outputs/PhyloglmResults/stepwise_boot500_", format(Sys.Date(), "%Y%m%d"), "\n")
cat("================================================================================\n")


source("/Users/kate/Downloads/create_bootstrap_forest_plot.r")
plot_data <- create_bootstrap_forest_plot(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "main_coefficients_comparison1.png"
)

detailed_data <- create_detailed_forest_plot(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "all_coefficients_forest_plot1.png"
)
