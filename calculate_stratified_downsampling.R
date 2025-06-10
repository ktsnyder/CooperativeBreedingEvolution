# calculate_stratified_downsampling.R
# Perform detailed stratified downsampling calculations for bias correction
# Kate Snyder
# Created: 2025-06-08

library(dplyr)
library(tidyr)

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
      
      # Algebra: (n_coop - x) = target_prop * (n_total - x)
      # n_coop - x = target_prop * n_total - target_prop * x
      # n_coop - target_prop * n_total = x - target_prop * x
      # n_coop - target_prop * n_total = x(1 - target_prop)
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
      name
    )
    
    report <- c(report,
      paste("\n", section_title, "\n"),
      "```",
      calc$calculation,
      "```"
    )
  }
  
  # Add recommendations
  report <- c(report,
    "\n## Downsampling Summary\n",
    "Based on these calculations, the following species should be removed:"
  )
  
  for (name in names(results$downsampling)) {
    calc <- results$downsampling[[name]]
    if (calc$n_to_remove > 0) {
      report <- c(report,
        paste("- **", gsub("_", " ", name), "**: Remove", 
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
    output_file = "claude_code_sessions/Stratified_Downsampling_Calculations.md"
  )
}