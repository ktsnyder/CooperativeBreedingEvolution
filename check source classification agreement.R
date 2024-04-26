# Load necessary libraries
library(dplyr)

data = CoopSubset

# List of binary classification columns
classification_cols <- c(
  'DunnCoop', 'BiagoliniCoop', 'DowningCoop', 'JetzCoop', 'RubensteinCoop',
  'CockburnCoop', 'ReihlCoop', 'Griesser2017Coop', 'BOWCoop', 'DaleCoop', 'CornwallisCoop'
)

# Initialize vectors to store the results
pairs <- c()
agreement_rates <- c()
disagreement_rates <- c()
num_species_with_classification <- c()  # New vector to store the count of species with classifications


# Function to calculate agreement, disagreement, and count of species with classification
calculate_agreement <- function(data, col1, col2) {
  valid_data <- na.omit(data[c(col1, col2)])
  agreement <- sum(valid_data[[col1]] == valid_data[[col2]])
  disagreement <- sum(valid_data[[col1]] != valid_data[[col2]])
  total <- nrow(valid_data)
  
  if (total > 0) {
    agreement_rate <- agreement / total
    disagreement_rate <- disagreement / total
  } else {
    agreement_rate <- NA
    disagreement_rate <- NA
  }
  
  list(agreement_rate = agreement_rate, disagreement_rate = disagreement_rate, num_species = total)
}

# Calculate pairwise agreement, disagreement, and count for all combinations
combn(classification_cols, 2, function(cols) {
  result <- calculate_agreement(data, cols[1], cols[2])
  pairs <- c(pairs, paste(cols[1], cols[2], sep = "-"))
  agreement_rates <- c(agreement_rates, result$agreement_rate)
  disagreement_rates <- c(disagreement_rates, result$disagreement_rate)
  num_species_with_classification <- c(num_species_with_classification, result$num_species)
  
  # Explicitly assign the updated vectors back to the global environment
  assign("pairs", pairs, envir = .GlobalEnv)
  assign("agreement_rates", agreement_rates, envir = .GlobalEnv)
  assign("disagreement_rates", disagreement_rates, envir = .GlobalEnv)
  assign("num_species_with_classification", num_species_with_classification, envir = .GlobalEnv)
}, simplify = FALSE)

# Create a data frame to hold the results
results_df <- data.frame(
  Pair = pairs,
  AgreementRate = agreement_rates,
  DisagreementRate = disagreement_rates,
  NumSpeciesWithClassification = num_species_with_classification  # Add the new column to the data frame
)

# Display all pairwise comparisons
print(results_df)
