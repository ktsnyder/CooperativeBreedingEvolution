# calculate_downsampling_function.R
# Function extracted from bias tests.R for territoriality downsampling
# This should be sourced when running territoriality bias correction

# Function to calculate number of species to downsample
# Based on a specified territoriality column and data presence columns
calculate_downsampling <- function(df, 
                                   territoriality_col,
                                   territory_value_high = 1,  # Value for high territoriality (e.g., year-round)
                                   territory_value_low = 0,   # Value for low territoriality (e.g., seasonal) 
                                   data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop"),
                                   target_proportion = NULL) {
  
  # Validate inputs
  if (!territoriality_col %in% colnames(df)) {
    stop(paste("Territoriality column", territoriality_col, "not found in dataframe"))
  }
  
  for (col in data_cols) {
    if (!col %in% colnames(df)) {
      stop(paste("Data column", col, "not found in dataframe"))
    }
  }
  
  print(territoriality_col)
  
  df[,territoriality_col] <- as.character(df[,territoriality_col])
  
  # Create a variable for having complete data across specified columns
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  # Filter to rows with non-NA territory values
  df_filtered <- df[!is.na(df[[territoriality_col]]), ]
  
  # Calculate data availability for high and low territoriality
  terr_high <- df_filtered[df_filtered[[territoriality_col]] == territory_value_high, ]
  terr_low <- df_filtered[df_filtered[[territoriality_col]] == territory_value_low, ]
  
  # Count species with complete data for each territoriality type
  n_complete_high <- sum(terr_high$has_complete_data)
  n_complete_low <- sum(terr_low$has_complete_data)
  
  # Total species counts for each territoriality type
  n_total_high <- nrow(terr_high)
  n_total_low <- nrow(terr_low)
  
  # Calculate current proportions
  prop_high <- n_complete_high / n_total_high
  prop_low <- n_complete_low / n_total_low
  
  # If target proportion is not specified, use the lower of the two
  if (is.null(target_proportion)) {
    target_proportion <- min(prop_high, prop_low)
  }
  
  # Calculate how many species need to be removed from the group with higher proportion
  if (prop_high > prop_low) {
    # Need to remove species from the high territoriality group
    target_n_complete_high <- target_proportion * n_total_high
    n_to_remove <- n_complete_high - target_n_complete_high
  } else if (prop_low > prop_high) {
    # Need to remove species from the low territoriality group
    target_n_complete_low <- target_proportion * n_total_low
    n_to_remove <- n_complete_low - target_n_complete_low
  } else {
    # Proportions are already equal
    n_to_remove <- 0
  }
  
  # samplingMatrixOut <- create_downsampling_matrix(df = df, territoriality_col = territoriality_col, territory_value_high = territory_value_high, territory_value_low = territory_value_low, data_cols = data_cols)
  
  # Create a detailed report
  report <- list(
    proportions = data.frame(
      territoriality = c("High", "Low"),
      species_with_data = c(n_complete_high, n_complete_low),
      total_species = c(n_total_high, n_total_low),
      proportion = c(prop_high, prop_low)
    ),
    target_proportion = target_proportion,
    downsample_info = if (prop_high > prop_low) {
      list(
        group_to_downsample = "High territoriality",
        territory_value = territory_value_high,
        n_to_remove = round(n_to_remove),
        original_proportion = prop_high,
        new_proportion = target_proportion
      )
    } else if (prop_low > prop_high) {
      list(
        group_to_downsample = "Low territoriality",
        territory_value = territory_value_low, 
        n_to_remove = round(n_to_remove),
        original_proportion = prop_low,
        new_proportion = target_proportion
      )
    } else {
      list(
        group_to_downsample = "None - already balanced",
        n_to_remove = 0
      )
    }
  )
  
  # Return the number to remove and the detailed report
  return(list(
    n_to_remove = ceiling(n_to_remove),
    report = report,
    filter_statement = if (n_to_remove > 0) {
      if (prop_high > prop_low) {
        paste0("To balance, remove ", ceiling(n_to_remove), " species with ", 
               territoriality_col, " = ", territory_value_high, 
               " and complete data for ", paste(data_cols, collapse = " and "))
      } else {
        paste0("To balance, remove ", ceiling(n_to_remove), " species with ", 
               territoriality_col, " = ", territory_value_low, 
               " and complete data for ", paste(data_cols, collapse = " and "))
      }
    } else {
      "No downsampling needed - proportions are already balanced"
    }
  ))
}