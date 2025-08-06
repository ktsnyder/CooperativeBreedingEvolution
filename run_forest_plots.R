# Run the fixed forest plot functions

source("create_bootstrap_forest_plot_fixed.R")

# Create main coefficients comparison plot
plot_data <- create_bootstrap_forest_plot(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "main_coefficients_comparison.png"
)

# Create detailed forest plot with all coefficients
detailed_data <- create_detailed_forest_plot(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "all_coefficients_forest_plot.png"
)

# You can also create filtered plots for specific predictors
# For example, just Geographic Region:
geo_data <- create_detailed_forest_plot(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "geographic_region_coefficients.png",
  predictor_filter = "Geographic Region"
)