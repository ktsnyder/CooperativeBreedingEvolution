# evaluate overlap in species data - song and life history

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)

boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')

# Read your data
bird_data <- read.csv(newdata)

# 1. Check data coverage for key variables
coverage_analysis <- function(data, main_vars, additional_vars) {
  # Function to analyze data coverage and overlap
  
  # Count non-NA values for each variable
  var_counts <- sapply(c(main_vars, additional_vars), function(var) {
    sum(!is.na(data[[var]]))
  })
  
  # Create a matrix showing overlap between main variables and additional variables
  overlap_matrix <- matrix(0, nrow = length(main_vars), ncol = length(additional_vars))
  rownames(overlap_matrix) <- main_vars
  colnames(overlap_matrix) <- additional_vars
  
  for (i in 1:length(main_vars)) {
    for (j in 1:length(additional_vars)) {
      overlap_matrix[i, j] <- sum(!is.na(data[[main_vars[i]]]) & !is.na(data[[additional_vars[j]]]))
    }
  }
  
  # Create overlap between key variables
  key_overlap <- matrix(0, nrow = length(main_vars), ncol = length(main_vars))
  rownames(key_overlap) <- main_vars
  colnames(key_overlap) <- main_vars
  
  for (i in 1:length(main_vars)) {
    for (j in 1:length(main_vars)) {
      key_overlap[i, j] <- sum(!is.na(data[[main_vars[i]]]) & !is.na(data[[main_vars[j]]]))
    }
  }
  
  return(list(var_counts = var_counts, 
              overlap_matrix = overlap_matrix,
              key_overlap = key_overlap))
}

# 2. Analyze correlations between variables
correlation_analysis <- function(data, vars_to_check) {
  # Subset data to complete cases for the specified variables
  complete_data <- data %>% 
    select(all_of(vars_to_check)) %>%
    filter(complete.cases(.))
  
  # Calculate correlations
  # Convert any factor variables to numeric first
  numeric_data <- complete_data %>% 
    mutate_if(is.factor, as.numeric) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric)
  
  # Calculate correlation matrix
  cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
  
  # Return results
  return(list(
    n_complete = nrow(complete_data),
    correlation_matrix = cor_matrix
  ))
}

# Define variable sets
main_vars <- c("FemaleSong_Agg01", "HighConfidence_Coop", "Territory", 
               "Social.bond", "Duet", "Chorus")

ecological_vars <- c("sedentariness_Griesser2023", "food_energy_Griesser2023", 
                     "food_h_level_Griesser2023", "grouping_Griesser2023", "Griesser2023.Colonial01", "egg_mass_Griesser2023", "weight_Griesser2023", 
                     "brain_Griesser2023", "time_fed_Griesser2023", 
                     "clutch_size_Griesser2023", "caretakers_Griesser2023", 
                     "longevity_Griesser2023", "social_bonds_Griesser2023", "social_system_incl_nk_coop_Griesser2017")

life_history_vars <- c("Final.polygyny", "Final.EPP")

song_vars <- c("Song.rep.final", "Syllable.rep.final", 
               "Duration.final")

# Run analyses
coverage_results <- coverage_analysis(bird_data, main_vars, 
                                      c(ecological_vars, life_history_vars, song_vars))

# Create subsets for correlation analyses
# Subset for ecological variables
eco_vars_to_check <- c(main_vars, ecological_vars)
eco_corr <- correlation_analysis(bird_data, eco_vars_to_check)

# Subset for life history variables
life_vars_to_check <- c(main_vars, life_history_vars)
life_corr <- correlation_analysis(bird_data, life_vars_to_check)

# Subset for song variables
song_vars_to_check <- c(main_vars, song_vars)
song_corr <- correlation_analysis(bird_data, song_vars_to_check)

# 3. Create territory binary variables as suggested
process_territory <- function(data) {
  data$Territory <- as.character(data$Territory)
  data$Territory <- factor(data$Territory, levels = c("1", "2", "3"))
  
  # Create binary variables as suggested
  if (!"Territory_12vs3" %in% colnames(data) || !"Territory_1vs23" %in% colnames(data)) {
    data$Territory_12vs3 <- ifelse(data$Territory == 3, 1, 0)
    data$Territory_1vs23 <- ifelse(data$Territory == 1, 0, 1)
  }
  
  # Create individual binary variables for each territory level
  if (!"Territory1" %in% colnames(data)) {
    territory_data <- data[!is.na(data$Territory), ]
    unique(territory_data$Territory)
    Territory.matrix <- model.matrix(~ Territory - 1, data = territory_data)
    
    # Add these back to the original dataframe
    data$Territory1 <- NA
    data$Territory2 <- NA
    data$Territory3 <- NA
    
    data$Territory1[!is.na(data$Territory)] <- Territory.matrix[, 1]
    data$Territory2[!is.na(data$Territory)] <- Territory.matrix[, 2]
    data$Territory3[!is.na(data$Territory)] <- Territory.matrix[, 3]
  }
  
  return(data)
}

bird_data <- process_territory(bird_data)

# 4. Explore the relationship between HighConfidence_Coop and Kin_NK
kin_coop_table <- bird_data %>%
  group_by(HighConfidence_Coop, Kin_NK) %>%
  summarize(count = n(), .groups = "drop")

# Print results
print("Variable coverage counts:")
print(coverage_results$var_counts)

print("Overlap between main variables and ecological variables:")
print(coverage_results$overlap_matrix[, 1:length(ecological_vars)])

print("Overlap between main variables and life history variables:")
print(coverage_results$overlap_matrix[, (length(ecological_vars)+1):(length(ecological_vars)+length(life_history_vars))])

print("Overlap between main variables and song variables:")
print(coverage_results$overlap_matrix[, (length(ecological_vars)+length(life_history_vars)+1):ncol(coverage_results$overlap_matrix)])

print("Overlap among key variables:")
print(coverage_results$key_overlap)

print("Number of complete cases for ecological variables:")
print(eco_corr$n_complete)

print("Correlations for ecological variables:")
print(round(eco_corr$correlation_matrix, 2))

print("Number of complete cases for life history variables:")
print(life_corr$n_complete)

print("Correlations for life history variables:")
print(round(life_corr$correlation_matrix, 2))

print("Number of complete cases for song variables:")
print(song_corr$n_complete)

print("Correlations for song variables:")
print(round(song_corr$correlation_matrix, 2))

print("Relationship between HighConfidence_Coop and Kin_NK:")
print(kin_coop_table)

# 5. Plot key correlations (optional)
# This creates correlation plots for visual inspection
plot_correlations <- function(corr_matrix, title) {
  corrplot(corr_matrix, method = "circle", type = "upper", 
           tl.col = "black", tl.srt = 45, 
           title = title, mar = c(0, 0, 2, 0), tl.cex = 0.3)
}

par(mfrow = c(2, 2))
plot_correlations(eco_corr$correlation_matrix, "Ecological Variables")
plot_correlations(life_corr$correlation_matrix, "Life History Variables")
plot_correlations(song_corr$correlation_matrix, "Song Variables")


#### column-wise data summary ----

bird_data = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET.csv")

# Load libraries
library(dplyr)

# Create a function to summarize each column
summarize_columns <- function(data) {
  summary_data <- data.frame(
    Column = character(),
    Class = character(),
    NumUnique = numeric(),
    NumNotNA = numeric(),
    UniqueValues = character(),
    stringsAsFactors = FALSE
  )
  
  for (col_name in names(data)) {
    col_data <- data[[col_name]]
    col_class <- class(col_data)[1]
    num_unique <- length(unique(col_data))
    num_not_na <- sum(!is.na(col_data))
    
    # Get unique values for display (if not too many)
    if (num_unique <= 10) {
      unique_vals <- paste(sort(unique(na.omit(col_data))), collapse = ", ")
    } else if (col_class %in% c("numeric", "integer", "double")) {
      # For numeric with many values, show range
      unique_vals <- paste("Range:", min(col_data, na.rm = TRUE), "to", max(col_data, na.rm = TRUE))
    } else {
      # For categorical with many values, show truncated list
      unique_vals <- paste(head(sort(unique(na.omit(col_data))), 5), collapse = ", ") 
      unique_vals <- paste0(unique_vals, ", ... (", num_unique, " total unique values)")
    }
    
    # Add to summary data frame
    summary_data <- rbind(summary_data, data.frame(
      Column = col_name,
      Class = col_class,
      NumUnique = num_unique,
      NumNotNA = num_not_na,
      UniqueValues = unique_vals,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_data)
}

# Apply function to your dataset
column_summary <- summarize_columns(bird_data)

# Print result in a format easy to share with LLM
print_summary <- function(summary_df) {
  cat("DATASET SUMMARY:\n\n")
  
  for (i in 1:nrow(summary_df)) {
    cat("Column:", summary_df$Column[i], "\n")
    cat("- Class:", summary_df$Class[i], "\n")
    cat("- Number of unique values:", summary_df$NumUnique[i], "\n")
    cat("- Number of non-NA values:", summary_df$NumNotNA[i], "\n")
    cat("- Unique values:", summary_df$UniqueValues[i], "\n\n")
  }
}

# Print the summary
print_summary(column_summary)

# Write summary to a file for easier sharing
writeLines(capture.output(print_summary(column_summary)), "Data_with_AVONET_column_summary.txt")
