# run_bias_tests.R
# Functions for testing and correcting geographic and territorial biases in comparative analyses
# Kate Snyder
# Created: 2025-06-08
# 
# This script provides functions to:
# 1. Calculate data availability biases across geographic regions and territoriality levels
# 2. Determine appropriate downsampling to correct for biases
# 3. Execute bias-corrected analyses
#
# Usage from Run_Analyses.R:
#   source("claude_code_sessions/run_bias_tests.R")
#   bias_results <- run_bias_tests(df = dfOs, tree = tree, output_dir = "Outputs/BiasTests")

library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)

#' Calculate downsampling requirements for bias correction
#'
#' @param df Data frame with species data
#' @param bias_type Either "geographic" or "territoriality" 
#' @param grouping_col Column name for grouping (e.g., "GeographicRegion_Jetz", "TerritorialityWeakVsStrong")
#' @param group_value_high Value for the group expected to have more data
#' @param group_value_low Value for the group expected to have less data
#' @param data_cols Vector of column names to check for data completeness
#' @param target_proportion Optional target proportion; if NULL, uses the lower of the two groups
#' @return List with downsampling report and parameters
calculate_downsampling <- function(df, 
                                 bias_type = c("geographic", "territoriality"),
                                 grouping_col,
                                 group_value_high,
                                 group_value_low,
                                 data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop"),
                                 target_proportion = NULL) {
  
  bias_type <- match.arg(bias_type)
  
  # Validate inputs
  if (!grouping_col %in% colnames(df)) {
    stop(paste("Grouping column", grouping_col, "not found in dataframe"))
  }
  
  for (col in data_cols) {
    if (!col %in% colnames(df)) {
      stop(paste("Data column", col, "not found in dataframe"))
    }
  }
  
  # Create a variable for having complete data across specified columns
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  # Filter to rows with non-NA group values
  df_filtered <- df[!is.na(df[[grouping_col]]), ]
  
  # Calculate data availability for high and low groups
  group_high <- df_filtered[df_filtered[[grouping_col]] == group_value_high, ]
  group_low <- df_filtered[df_filtered[[grouping_col]] == group_value_low, ]
  
  # Count species with complete data for each group
  n_complete_high <- sum(group_high$has_complete_data)
  n_complete_low <- sum(group_low$has_complete_data)
  
  # Total species counts for each group
  n_total_high <- nrow(group_high)
  n_total_low <- nrow(group_low)
  
  # Calculate current proportions
  prop_high <- n_complete_high / n_total_high
  prop_low <- n_complete_low / n_total_low
  
  # If target proportion is not specified, use the lower of the two
  if (is.null(target_proportion)) {
    target_proportion <- min(prop_high, prop_low)
  }
  
  # Calculate how many species need to be removed from the group with higher proportion
  if (prop_high > prop_low) {
    # Need to remove species from the high group
    target_n_complete_high <- target_proportion * n_total_high
    n_to_remove <- round(n_complete_high - target_n_complete_high)
    group_to_downsample <- group_value_high
  } else if (prop_low > prop_high) {
    # Need to remove species from the low group
    target_n_complete_low <- target_proportion * n_total_low
    n_to_remove <- round(n_complete_low - target_n_complete_low)
    group_to_downsample <- group_value_low
  } else {
    # Proportions are already equal
    n_to_remove <- 0
    group_to_downsample <- NA
  }
  
  # Create detailed report
  report <- list(
    bias_type = bias_type,
    grouping_col = grouping_col,
    proportions = data.frame(
      group = c(group_value_high, group_value_low),
      label = c("High", "Low"),
      species_with_data = c(n_complete_high, n_complete_low),
      total_species = c(n_total_high, n_total_low),
      proportion = c(prop_high, prop_low),
      stringsAsFactors = FALSE
    ),
    target_proportion = target_proportion,
    downsample_info = list(
      group_to_downsample = group_to_downsample,
      n_to_remove = n_to_remove,
      original_proportion = ifelse(is.na(group_to_downsample), NA, 
                                 ifelse(group_to_downsample == group_value_high, prop_high, prop_low)),
      new_proportion = target_proportion
    )
  )
  
  return(report)
}

#' Create a summary matrix of data availability by groups
#'
#' @param df Data frame with species data
#' @param grouping_col Column name for grouping
#' @param data_cols Vector of column names to check for data completeness
#' @return Data frame with summary statistics
create_bias_summary_matrix <- function(df, grouping_col, data_cols) {
  
  # Ensure grouping column is character
  df[[grouping_col]] <- as.character(df[[grouping_col]])
  
  # Create complete data indicator
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  # Create summary matrix
  summary_matrix <- df %>%
    filter(!is.na(!!sym(grouping_col))) %>%
    group_by(!!sym(grouping_col)) %>%
    summarise(
      n_with_data = sum(has_complete_data),
      n_without_data = sum(!has_complete_data),
      n_total = n(),
      prop_with_data = n_with_data / n_total,
      .groups = "drop"
    )
  
  # Add totals row
  totals <- df %>%
    filter(!is.na(!!sym(grouping_col))) %>%
    summarise(
      grouping_col = "Total",
      n_with_data = sum(has_complete_data),
      n_without_data = sum(!has_complete_data),
      n_total = n(),
      prop_with_data = n_with_data / n_total
    )
  names(totals)[1] <- grouping_col
  
  summary_matrix <- bind_rows(summary_matrix, totals)
  
  return(summary_matrix)
}

#' Run comprehensive bias tests
#'
#' @param df Data frame with species data
#' @param tree Phylogenetic tree
#' @param output_dir Directory for saving results
#' @param geographic_col Column name for geographic regions (default: "GeographicRegion_Jetz")
#' @param territoriality_cols Vector of column names for territoriality measures
#' @param dimorphism_cols Vector of column names for sexual dimorphism measures
#' @param dimorphism_thresholds List of thresholds for dimorphism measures (default: median split)
#' @param data_cols Vector of column names to check for data completeness
#' @param save_plots Whether to save diagnostic plots
#' @return List of bias test results
run_bias_tests <- function(df,
                          tree,
                          output_dir = "Outputs/BiasTests",
                          geographic_col = "GeographicRegion_Jetz",
                          territoriality_cols = c("TerritorialityWeakVsStrong", 
                                                "TerritorialityWeakVsStrongHighConf",
                                                "Territory_12vs3"),
                          dimorphism_cols = c("logMaleFemalePlumageDiffAbs", 
                                            "PercentAbsLogWingDimorphism"),
                          dimorphism_thresholds = list(
                            logMaleFemalePlumageDiffAbs = 0.5,  # median split by default
                            PercentLogWingDimorphism_AVONET = 0.5  # median split by default
                          ),
                          data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop"),
                          save_plots = TRUE) {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Prepare data
  df_work <- df
  
  # Add geographic regions if needed
  if ("GeographicRegion_Jetz" %in% geographic_col && !"GeographicRegion_Jetz" %in% colnames(df_work)) {
    df_work$GeographicRegion_Jetz <- NA
    df_work$GeographicRegion_Jetz[which(df_work$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
    df_work$GeographicRegion_Jetz[which(df_work$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
  }
  
  # Filter to species in tree
  df_work <- df_work[df_work$species %in% tree$tip.label, ]
  
  results <- list()
  
  # Test geographic bias
  if (!is.null(geographic_col) && geographic_col %in% colnames(df_work)) {
    cat("\nTesting geographic bias using", geographic_col, "\n")
    
    geo_matrix <- create_bias_summary_matrix(df_work, geographic_col, data_cols)
    print(geo_matrix)
    
    geo_downsample <- calculate_downsampling(
      df = df_work,
      bias_type = "geographic",
      grouping_col = geographic_col,
      group_value_high = "Holarctic",
      group_value_low = "Tropical",
      data_cols = data_cols
    )
    
    results$geographic <- list(
      summary_matrix = geo_matrix,
      downsample_report = geo_downsample
    )
    
    # Save results
    write.csv(geo_matrix, 
              file.path(output_dir, "geographic_bias_summary.csv"),
              row.names = FALSE)
  }
  
  # Test territoriality biases
  for (terr_col in territoriality_cols) {
    if (terr_col %in% colnames(df_work)) {
      cat("\nTesting territoriality bias using", terr_col, "\n")
      
      # Convert to character for consistent handling
      df_work[[terr_col]] <- as.character(df_work[[terr_col]])
      
      terr_matrix <- create_bias_summary_matrix(df_work, terr_col, data_cols)
      print(terr_matrix)
      
      # Determine high/low values based on column
      if (terr_col == "Territory_12vs3") {
        high_val <- "1"  # Territory = 3 (year-round)
        low_val <- "0"   # Territory = 1 or 2 (seasonal)
      } else {
        high_val <- "1"  # Strong territoriality
        low_val <- "0"   # Weak territoriality
      }
      
      terr_downsample <- calculate_downsampling(
        df = df_work,
        bias_type = "territoriality",
        grouping_col = terr_col,
        group_value_high = high_val,
        group_value_low = low_val,
        data_cols = data_cols
      )
      
      results[[terr_col]] <- list(
        summary_matrix = terr_matrix,
        downsample_report = terr_downsample
      )
      
      # Save results
      write.csv(terr_matrix,
                file.path(output_dir, paste0(terr_col, "_bias_summary.csv")),
                row.names = FALSE)
    }
  }
  
  # Test sexual dimorphism biases
  for (dim_col in dimorphism_cols) {
    if (dim_col %in% colnames(df_work)) {
      cat("\nTesting sexual dimorphism bias using", dim_col, "\n")
      
      # Get non-NA values for this dimorphism measure
      dim_values <- df_work[[dim_col]][!is.na(df_work[[dim_col]])]
      
      # Determine threshold (use provided or calculate median)
      if (dim_col %in% names(dimorphism_thresholds) && 
          dimorphism_thresholds[[dim_col]] != 0.5) {
        threshold <- dimorphism_thresholds[[dim_col]]
      } else {
        threshold <- median(dim_values, na.rm = TRUE)
      }
      
      # Create binary variable for high/low dimorphism
      df_work$dimorphism_binary <- NA
      df_work$dimorphism_binary[!is.na(df_work[[dim_col]])] <- 
        ifelse(df_work[[dim_col]][!is.na(df_work[[dim_col]])] > threshold, "High", "Low")
      
      # Create summary matrix
      dim_matrix <- create_bias_summary_matrix(df_work, "dimorphism_binary", data_cols)
      # Add threshold info to matrix
      attr(dim_matrix, "threshold") <- threshold
      attr(dim_matrix, "original_column") <- dim_col
      
      cat("  Threshold used:", round(threshold, 3), "\n")
      print(dim_matrix)
      
      # Calculate downsampling (use "territoriality" as bias_type for compatibility)
      dim_downsample <- calculate_downsampling(
        df = df_work,
        bias_type = "territoriality",
        grouping_col = "dimorphism_binary",
        group_value_high = "High",
        group_value_low = "Low",
        data_cols = data_cols
      )
      
      results[[dim_col]] <- list(
        summary_matrix = dim_matrix,
        downsample_report = dim_downsample,
        threshold = threshold
      )
      
      # Save results
      write.csv(dim_matrix,
                file.path(output_dir, paste0(dim_col, "_bias_summary.csv")),
                row.names = FALSE)
    }
  }
  
  # Create summary plots if requested
  if (save_plots) {
    create_bias_diagnostic_plots(results, output_dir)
  }
  
  # Save full report
  saveRDS(results, file.path(output_dir, "bias_test_results.rds"))
  
  # Print summary
  cat("\n=== BIAS TEST SUMMARY ===\n")
  for (test_name in names(results)) {
    report <- results[[test_name]]$downsample_report
    if (report$downsample_info$n_to_remove > 0) {
      cat("\n", test_name, ":\n", sep = "")
      cat("  Group to downsample:", report$downsample_info$group_to_downsample, "\n")
      cat("  Species to remove:", report$downsample_info$n_to_remove, "\n")
      cat("  Current proportion:", round(report$downsample_info$original_proportion, 3), "\n")
      cat("  Target proportion:", round(report$target_proportion, 3), "\n")
    } else {
      cat("\n", test_name, ": No bias detected (proportions are equal)\n", sep = "")
    }
  }
  
  return(results)
}

#' Create diagnostic plots for bias tests
#'
#' @param results Results from run_bias_tests
#' @param output_dir Directory for saving plots
create_bias_diagnostic_plots <- function(results, output_dir) {
  
  # Create a combined plot showing all biases
  plot_data <- data.frame()
  
  for (test_name in names(results)) {
    test_results <- results[[test_name]]
    props <- test_results$summary_matrix
    props <- props[props[[1]] != "Total", ]  # Remove total row
    
    temp_df <- data.frame(
      test = test_name,
      group = props[[1]],
      proportion = props$prop_with_data,
      n_total = props$n_total,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, temp_df)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = group, y = proportion, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(proportion * 100, 1), "%\n(n=", n_total, ")")),
              vjust = -0.5) +
    facet_wrap(~ test, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Data Availability Biases",
         subtitle = "Proportion of species with complete data for FemaleSong and CooperativeBreeding",
         y = "Proportion with complete data",
         x = "") +
    ylim(0, max(plot_data$proportion) * 1.2) +
    guides(fill = "none")
  
  ggsave(file.path(output_dir, "bias_diagnostic_plots.pdf"), p, 
         width = 10, height = 6)
  ggsave(file.path(output_dir, "bias_diagnostic_plots.png"), p, 
         width = 10, height = 6, dpi = 150)
}

#' Get species to exclude for downsampling
#'
#' @param df Data frame with species data
#' @param downsample_report Report from calculate_downsampling
#' @param seed Random seed for reproducibility
#' @return Vector of species names to exclude
get_species_to_exclude <- function(df, downsample_report, seed = NULL) {
  
  if (downsample_report$downsample_info$n_to_remove == 0) {
    return(character(0))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Get the grouping column and value to downsample
  grouping_col <- downsample_report$grouping_col
  group_value <- downsample_report$downsample_info$group_to_downsample
  n_to_remove <- downsample_report$downsample_info$n_to_remove
  
  # Create complete data indicator
  data_cols <- c("FemaleSong_Agg01", "HighConfidence_Coop")  # Should be passed as parameter
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  # Get species in the group to downsample that have complete data
  eligible_species <- df$species[df[[grouping_col]] == group_value & df$has_complete_data]
  
  # Randomly select species to exclude
  if (length(eligible_species) < n_to_remove) {
    warning("Not enough eligible species to remove. Removing all eligible species.")
    return(eligible_species)
  }
  
  species_to_exclude <- sample(eligible_species, n_to_remove, replace = FALSE)
  
  return(species_to_exclude)
}