# calculate_stratified_downsampling_with_territoriality.R
# Perform detailed stratified downsampling calculations for bias correction
# Including territoriality bias calculations
# Kate Snyder
# Created: 2025-06-08

library(dplyr)
library(tidyr)

#' Calculate territoriality bias downsampling
#'
#' @param df Data frame with species data
#' @param territoriality_col Column name for territoriality data
#' @param territory_value_high Value representing high territoriality
#' @param territory_value_low Value representing low territoriality
#' @param data_cols Columns that must have data (e.g., c("FemaleSong_Agg01", "HighConfidence_Coop"))
#' @return List with downsampling information
calculate_territoriality_downsampling <- function(df,
                                                territoriality_col = "TerritorialityWeakVsStrong",
                                                territory_value_high = "1",
                                                territory_value_low = "0",
                                                data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")) {
  
  # Check if required columns exist
  missing_cols <- setdiff(c(territoriality_col, data_cols), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add HaveData indicator (all specified data columns must be non-NA)
  df$HaveData <- apply(df[data_cols], 1, function(x) all(!is.na(x)))
  
  # Convert territoriality to character for comparison
  df[[territoriality_col]] <- as.character(df[[territoriality_col]])
  
  # Calculate proportions for high and low territoriality
  high_terr <- df %>%
    filter(!!sym(territoriality_col) == territory_value_high) %>%
    summarise(
      species_with_data = sum(HaveData),
      total_species = n(),
      proportion = ifelse(n() > 0, species_with_data / total_species, 0)
    )
  
  low_terr <- df %>%
    filter(!!sym(territoriality_col) == territory_value_low) %>%
    summarise(
      species_with_data = sum(HaveData),
      total_species = n(),
      proportion = ifelse(n() > 0, species_with_data / total_species, 0)
    )
  
  # Check if we have data for both groups
  if (nrow(high_terr) == 0 || high_terr$total_species == 0) {
    warning(paste("No species found with", territoriality_col, "=", territory_value_high))
    high_terr <- data.frame(species_with_data = 0, total_species = 0, proportion = 0)
  }
  
  if (nrow(low_terr) == 0 || low_terr$total_species == 0) {
    warning(paste("No species found with", territoriality_col, "=", territory_value_low))
    low_terr <- data.frame(species_with_data = 0, total_species = 0, proportion = 0)
  }
  
  # Create proportions data frame
  proportions <- data.frame(
    territoriality = c("High", "Low"),
    species_with_data = c(high_terr$species_with_data, low_terr$species_with_data),
    total_species = c(high_terr$total_species, low_terr$total_species),
    proportion = c(high_terr$proportion, low_terr$proportion)
  )
  
  # Check if we can perform downsampling
  if (high_terr$total_species == 0 && low_terr$total_species == 0) {
    return(list(
      proportions = proportions,
      target_proportion = 0,
      n_to_remove = 0,
      downsample_info = list(
        group_to_downsample = "None - no data available",
        territory_value = NA
      ),
      calculation = "No species found with specified territoriality values"
    ))
  }
  
  # Determine which group has higher proportion and needs downsampling
  if (high_terr$proportion > low_terr$proportion) {
    # Need to downsample high territoriality group
    target_proportion <- low_terr$proportion
    n_to_remove <- high_terr$species_with_data - round(high_terr$total_species * target_proportion)
    group_to_downsample <- "High territoriality"
    territory_value <- territory_value_high
  } else {
    # Need to downsample low territoriality group
    target_proportion <- high_terr$proportion
    n_to_remove <- low_terr$species_with_data - round(low_terr$total_species * target_proportion)
    group_to_downsample <- "Low territoriality"
    territory_value <- territory_value_low
  }
  
  # Make sure n_to_remove is not negative
  n_to_remove <- max(0, n_to_remove)
  
  # Create report
  report <- list(
    proportions = proportions,
    target_proportion = target_proportion,
    n_to_remove = n_to_remove,
    downsample_info = list(
      group_to_downsample = group_to_downsample,
      territory_value = territory_value
    ),
    calculation = paste0(
      "High territoriality: ", high_terr$species_with_data, "/", high_terr$total_species,
      " = ", round(high_terr$proportion, 3), "\n",
      "Low territoriality: ", low_terr$species_with_data, "/", low_terr$total_species,
      " = ", round(low_terr$proportion, 3), "\n",
      "Target proportion: ", round(target_proportion, 3), "\n",
      "Remove ", n_to_remove, " species from ", group_to_downsample
    )
  )
  
  return(report)
}

#' Calculate stratified downsampling for specific trait combinations
#'
#' @param df Data frame with species data
#' @param stratify_vars Vector of column names to stratify by
#' @param data_col Column to check for data presence (e.g., "FemaleSong_Agg01")
#' @param output_file Optional markdown file to save results
#' @return List with downsampling calculations
calculate_stratified_downsampling <- function(df,
                                            stratify_vars = c("GeographicRegion_Jetz", "HighConfidence_Coop"),
                                            data_col = "FemaleSong_Agg01",
                                            output_file = NULL) {
  
  # Add data availability indicator
  # For cooperative breeding analyses, need BOTH FS and CB data
  if ("HighConfidence_Coop" %in% colnames(df) && "HighConfidence_Coop" %in% stratify_vars) {
    df$HaveData <- !is.na(df[[data_col]]) & !is.na(df$HighConfidence_Coop)
  } else {
    df$HaveData <- !is.na(df[[data_col]])
  }
  
  # Create summary by strata
  strata_summary <- df %>%
    filter(!is.na(!!sym(stratify_vars[1])) & !is.na(!!sym(stratify_vars[2]))) %>%
    group_by(across(all_of(stratify_vars))) %>%
    summarise(
      n_with_data = sum(HaveData),
      n_total = n(),
      prop_with_data = n_with_data / n_total,
      .groups = "drop"
    ) %>%
    arrange(desc(prop_with_data))
  
  results <- list()
  results$summary <- strata_summary
  
  # Calculate specific downsampling scenarios
  downsampling_calcs <- list()
  
  # 1. Holarctic Non-cooperative vs Tropical Non-cooperative
  if ("GeographicRegion_Jetz" %in% stratify_vars & "HighConfidence_Coop" %in% stratify_vars) {
    hol_noncoop <- strata_summary %>% 
      filter(GeographicRegion_Jetz == "Holarctic" & HighConfidence_Coop == 0)
    trop_noncoop <- strata_summary %>% 
      filter(GeographicRegion_Jetz == "Tropical" & HighConfidence_Coop == 0)
    
    if (nrow(hol_noncoop) > 0 & nrow(trop_noncoop) > 0) {
      target_prop <- trop_noncoop$prop_with_data[1]
      n_to_keep <- round(hol_noncoop$n_total[1] * target_prop)
      n_to_remove <- hol_noncoop$n_with_data[1] - n_to_keep
      
      downsampling_calcs$holarctic_noncoop <- list(
        current = hol_noncoop,
        target = trop_noncoop,
        n_to_remove = n_to_remove,
        calculation = paste0(
          "Holarctic Non-coop: ", hol_noncoop$n_with_data[1], "/", hol_noncoop$n_total[1],
          " = ", round(hol_noncoop$prop_with_data[1], 3), "\n",
          "Tropical Non-coop: ", trop_noncoop$n_with_data[1], "/", trop_noncoop$n_total[1],
          " = ", round(trop_noncoop$prop_with_data[1], 3), "\n",
          "To balance: x/", hol_noncoop$n_total[1], " = ", round(target_prop, 3),
          " → x = ", n_to_keep, "\n",
          "Remove: ", hol_noncoop$n_with_data[1], " - ", n_to_keep, " = ", n_to_remove, " species"
        )
      )
    }
    
    # 2. Tropical Cooperative overrepresentation
    trop_coop <- strata_summary %>% 
      filter(GeographicRegion_Jetz == "Tropical" & HighConfidence_Coop == 1)
    
    if (nrow(trop_coop) > 0) {
      # Get overall tropical proportion (only among species with CB data)
      trop_all <- df %>%
        filter(GeographicRegion_Jetz == "Tropical" & !is.na(HighConfidence_Coop)) %>%
        summarise(
          n_coop = sum(HighConfidence_Coop == 1, na.rm = TRUE),
          n_total = n(),
          prop_coop = n_coop / n_total
        )
      
      # Among species with data
      trop_with_data <- df %>%
        filter(GeographicRegion_Jetz == "Tropical" & HaveData) %>%
        summarise(
          n_coop = sum(HighConfidence_Coop == 1, na.rm = TRUE),
          n_total = n(),
          prop_coop = n_coop / n_total
        )
      
      # Solve for x: (n_coop - x) / (n_total - x) = target_prop
      n_coop <- trop_with_data$n_coop[1]
      n_total <- trop_with_data$n_total[1]
      target_prop <- trop_all$prop_coop[1]
      
      x <- (n_coop - target_prop * n_total) / (1 - target_prop)
      
      downsampling_calcs$tropical_coop <- list(
        current_prop = trop_with_data$prop_coop[1],
        target_prop = target_prop,
        n_to_remove = round(x),
        calculation = paste0(
          "Tropical+Coop1+HaveData: ", n_coop, " species\n",
          "Tropical+HaveData (total): ", n_total, " species\n",
          "Overall Tropical Coop proportion: ", round(target_prop, 3), "\n",
          "Current proportion among those with data: ", round(trop_with_data$prop_coop[1], 3), "\n",
          "Solve: (", n_coop, " - x) / (", n_total, " - x) = ", round(target_prop, 3), "\n",
          "x = ", round(x, 1), " → Remove ", round(x), " species"
        )
      )
    }
  }
  
  # 3. Global Cooperative overrepresentation
  global_all <- df %>%
    filter(!is.na(HighConfidence_Coop)) %>%
    summarise(
      n_coop = sum(HighConfidence_Coop == 1, na.rm = TRUE),
      n_total = n(),
      prop_coop = n_coop / n_total
    )
  
  global_with_data <- df %>%
    filter(HaveData) %>%
    summarise(
      n_coop = sum(HighConfidence_Coop == 1, na.rm = TRUE),
      n_total = n(),
      prop_coop = n_coop / n_total
    )
  
  if (global_with_data$n_coop[1] > 0) {
    n_coop <- global_with_data$n_coop[1]
    n_total <- global_with_data$n_total[1]
    target_prop <- global_all$prop_coop[1]
    
    x <- (n_coop - target_prop * n_total) / (1 - target_prop)
    
    downsampling_calcs$global_coop <- list(
      current_prop = global_with_data$prop_coop[1],
      target_prop = target_prop,
      n_to_remove = round(x),
      calculation = paste0(
        "Global Coop1 rate: ", round(target_prop, 3), "\n",
        "Rate among those with data: ", round(global_with_data$prop_coop[1], 3), "\n",
        "Coop1 + HaveData: ", n_coop, " species\n",
        "HaveData (total): ", n_total, " species\n",
        "Solve: (", n_coop, " - x) / (", n_total, " - x) = ", round(target_prop, 3), "\n",
        "x = ", round(x, 1), " → Remove ", round(x), " species"
      )
    )
  }
  
  # 4. Territoriality bias calculations
  if ("TerritorialityWeakVsStrong" %in% colnames(df)) {
    terr_result <- calculate_territoriality_downsampling(
      df = df,
      territoriality_col = "TerritorialityWeakVsStrong",
      territory_value_high = "1",
      territory_value_low = "0",
      data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")
    )
    
    downsampling_calcs$territoriality_bias <- list(
      current_prop = terr_result$proportions,
      target_prop = terr_result$target_proportion,
      n_to_remove = terr_result$n_to_remove,
      calculation = terr_result$calculation
    )
  }
  
  # 5. Territory_12vs3 calculations if available
  if ("Territory_12vs3" %in% colnames(df)) {
    # For Territory_12vs3, use values 1 for weak/non-territorial and 3 for strong
    terr123_result <- calculate_territoriality_downsampling(
      df = df,
      territoriality_col = "Territory_12vs3",
      territory_value_high = "1",  # Strong territoriality
      territory_value_low = "0",   # Weak/non-territorial
      data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")
    )
    
    downsampling_calcs$territory_12vs3_bias <- list(
      current_prop = terr123_result$proportions,
      target_prop = terr123_result$target_proportion,
      n_to_remove = terr123_result$n_to_remove,
      calculation = terr123_result$calculation
    )
  }
  
  results$downsampling <- downsampling_calcs
  
  # Generate report if requested
  if (!is.null(output_file)) {
    generate_downsampling_report(results, output_file)
  }
  
  return(results)
}

#' Generate markdown report for stratified downsampling
generate_downsampling_report <- function(results, output_file) {
  
  report <- c(
    "# Stratified Downsampling Calculations",
    paste("\n**Date:**", Sys.Date()),
    "\n## Summary Table\n",
    "```",
    capture.output(print(results$summary, n = 100)),
    "```",
    "\n## Detailed Downsampling Calculations\n"
  )
  
  # Add each calculation
  for (name in names(results$downsampling)) {
    calc <- results$downsampling[[name]]
    
    section_title <- switch(name,
      holarctic_noncoop = "### Holarctic Non-cooperative Birds",
      tropical_coop = "### Tropical Cooperative Birds", 
      global_coop = "### Global Cooperative Birds",
      territoriality_bias = "### Territoriality Bias (TerritorialityWeakVsStrong)",
      territory_12vs3_bias = "### Territoriality Bias (Territory_12vs3)",
      name
    )
    
    report <- c(report,
      paste("\n", section_title, "\n"),
      "```",
      calc$calculation,
      "```"
    )
    
    # Add proportion table for territoriality
    if (name %in% c("territoriality_bias", "territory_12vs3_bias")) {
      report <- c(report,
        "\n**Proportions Table:**\n",
        "```",
        capture.output(print(calc$current_prop)),
        "```"
      )
    }
  }
  
  # Add recommendations
  report <- c(report,
    "\n## Downsampling Summary\n",
    "Based on these calculations, the following species should be removed:"
  )
  
  for (name in names(results$downsampling)) {
    calc <- results$downsampling[[name]]
    if (calc$n_to_remove > 0) {
      section_name <- switch(name,
        holarctic_noncoop = "Holarctic non-cooperative",
        tropical_coop = "Tropical cooperative", 
        global_coop = "Global cooperative",
        territoriality_bias = "TerritorialityWeakVsStrong",
        territory_12vs3_bias = "Territory_12vs3",
        gsub("_", " ", name)
      )
      
      report <- c(report,
        paste("- **", section_name, "**: Remove", 
              calc$n_to_remove, "species")
      )
    }
  }
  
  writeLines(report, output_file)
  cat("Report saved to:", output_file, "\n")
}

# Test the function if called directly
if (!interactive()) {
  library(phytools)
  
  # Load data
  df <- read.csv("Data_R_2025-06-07.csv")
  tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  
  # Filter to tree species
  df <- df[df$species %in% tree$tip.label, ]
  
  # Add geographic regions if needed
  if (!"GeographicRegion_Jetz" %in% colnames(df)) {
    df$GeographicRegion_Jetz <- NA
    df$GeographicRegion_Jetz[which(df$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
    df$GeographicRegion_Jetz[which(df$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
  }
  
  # Calculate stratified downsampling
  results <- calculate_stratified_downsampling(
    df = df,
    stratify_vars = c("GeographicRegion_Jetz", "HighConfidence_Coop"),
    data_col = "FemaleSong_Agg01",
    output_file = "Stratified_Downsampling_Calculations.md"
  )
}