library(phylopath)
source('create_phylopath_bias_robustness_figure.R')
source('create_forest_plot_with_counts.R')
source('create_bias_model_consistency_heatmap.R')

# Create forest plot with new labels
cat('Creating forest plot with updated labels...\n')
forest_plot <- create_forest_plot_with_counts(
  rate_to_plot = 'COOP->FS',
  output_file = 'forest_plot_new_labels.png'
)

# Create heatmap with new labels
cat('\nCreating heatmap with updated labels...\n')
heatmap <- create_bias_model_consistency_heatmap(
  output_file = 'heatmap_new_labels.png'
)

cat('\nBoth plots created with updated labels!\n')