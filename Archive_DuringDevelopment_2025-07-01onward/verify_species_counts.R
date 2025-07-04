source('create_phylopath_bias_robustness_figure.R')
library(dplyr)

# Check that species counts are consistent across iterations
results <- extract_phylopath_results()

cat('Verifying species counts are constant across iterations:\n\n')

for (bias_name in names(results)) {
  result <- results[[bias_name]]
  if (!is.null(result$detailed_models)) {
    species_counts <- result$detailed_models %>%
      group_by(seed) %>%
      summarise(n_species = first(nSpecies), .groups = 'drop')
    
    unique_counts <- unique(species_counts$n_species)
    
    cat(paste(bias_name, ':', unique_counts, 'species in every iteration\n'))
  }
}