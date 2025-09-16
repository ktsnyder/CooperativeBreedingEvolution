## Helper functions for transition count pipeline
## Kate Snyder / Claude
## 2025-06-28

# Source debug control if not already loaded
if (!exists("debugPrint")) {
  if (file.exists("debug_control.R")) {
    source("debug_control.R")
  } else {
    # Fallback debug functions
    debugPrint <- function(...) { }  # No-op by default
    validationPrint <- function(...) { cat("VALIDATION:", ..., "\n") }
  }
}

# Map generic column names to arrow plot format (q12, q13, etc.)
# CRITICAL: This mapping must be exact for transition_plot.R to work correctly
mapGenericToArrowNames <- function(generic_names) {
  # Create mapping from generic names to q codes
  # Based on the original mapping in excerpt_from_RunAnalyses_simmap overlap execution.R
  mapping <- c(
    "trait2_0to1_in_trait1_0" = "q12",  # FS0to1inCoop0
    "trait1_0to1_in_trait2_0" = "q13",  # Coop0to1inFS0
    "trait2_1to0_in_trait1_0" = "q21",  # FS1to0inCoop0
    "trait1_0to1_in_trait2_1" = "q24",  # Coop0to1inFS1
    "trait1_1to0_in_trait2_0" = "q31",  # Coop1to0inFS0
    "trait2_0to1_in_trait1_1" = "q34",  # FS0to1inCoop1
    "trait1_1to0_in_trait2_1" = "q42",  # Coop1to0inFS1
    "trait2_1to0_in_trait1_1" = "q43"   # FS1to0inCoop1
  )
  
  # Return the q code for each generic name
  result <- mapping[generic_names]
  names(result) <- generic_names
  return(result)
}

# Reverse mapping from q codes to generic names
mapArrowToGenericNames <- function(q_codes) {
  # Create reverse mapping
  mapping <- c(
    "q12" = "trait2_0to1_in_trait1_0",
    "q13" = "trait1_0to1_in_trait2_0",
    "q21" = "trait2_1to0_in_trait1_0",
    "q24" = "trait1_0to1_in_trait2_1",
    "q31" = "trait1_1to0_in_trait2_0",
    "q34" = "trait2_0to1_in_trait1_1",
    "q42" = "trait1_1to0_in_trait2_1",
    "q43" = "trait2_1to0_in_trait1_1"
  )
  
  result <- mapping[q_codes]
  names(result) <- q_codes
  return(result)
}

# Create a dataframe mapping generic names to legacy names for a specific trait pair
createTransitionNameMapping <- function(trait1_name, trait2_name) {
  # This function helps with debugging and validation
  mapping_df <- data.frame(
    generic_name = c(
      "trait1_0to1_in_trait2_0", "trait1_1to0_in_trait2_0",
      "trait1_0to1_in_trait2_1", "trait1_1to0_in_trait2_1",
      "trait2_0to1_in_trait1_0", "trait2_1to0_in_trait1_0",
      "trait2_0to1_in_trait1_1", "trait2_1to0_in_trait1_1"
    ),
    arrow_code = c(
      "q13", "q31",
      "q24", "q42",
      "q12", "q21",
      "q34", "q43"
    ),
    stringsAsFactors = FALSE
  )
  
  # Add example legacy names if traits are Coop and FS
  if (grepl("Coop", trait1_name, ignore.case = TRUE) && grepl("FS|FemaleSong", trait2_name, ignore.case = TRUE)) {
    mapping_df$legacy_name <- c(
      "Coop0to1inFS0", "Coop1to0inFS0",
      "Coop0to1inFS1", "Coop1to0inFS1",
      "FS0to1inCoop0", "FS1to0inCoop0",
      "FS0to1inCoop1", "FS1to0inCoop1"
    )
  }
  
  return(mapping_df)
}

# Prepare transition data for plotting
prepareTransitionDataForPlot <- function(df, trait1_name, trait2_name) {
  # First, check if we need to add ObsProp columns for transition_plot.R compatibility
  if (!any(grepl("^ObsProp", names(df)))) {
    # Add ObsProp columns that transition_plot.R expects
    # These are the overlap proportions with legacy naming
    if ("Overlap_0_0" %in% names(df)) {
      df$ObsProp0Absent <- df$Overlap_0_0  # trait1=0, trait2=0
      df$ObsProp1Absent <- df$Overlap_0_1  # trait1=0, trait2=1
      df$ObsProp0Present <- df$Overlap_1_0  # trait1=1, trait2=0
      df$ObsProp1Present <- df$Overlap_1_1  # trait1=1, trait2=1
    }
  }
  
  # Get the transition columns (excluding Expected columns)
  # Handle both regular names and names with .x suffix (from merge)
  transition_cols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01](\\.x)?$", names(df), value = TRUE)
  
  # If we have .x columns, rename them to remove the suffix
  if (any(grepl("\\.x$", transition_cols))) {
    for (col in transition_cols) {
      if (grepl("\\.x$", col)) {
        new_name <- sub("\\.x$", "", col)
        names(df)[names(df) == col] <- new_name
      }
    }
    # Update the list after renaming
    transition_cols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01]$", names(df), value = TRUE)
  }
  
  # Removed debug output - issue resolved
  
  # Handle Expected columns that might have .y suffix
  expected_cols <- grep("Expected(\\.y)?$", names(df), value = TRUE)
  if (any(grepl("\\.y$", expected_cols))) {
    for (col in expected_cols) {
      if (grepl("\\.y$", col)) {
        new_name <- sub("\\.y$", "", col)
        names(df)[names(df) == col] <- new_name
      }
    }
  }
  
  # Calculate differences from expected for each transition
  diff_created <- 0
  for (col in transition_cols) {
    expected_col <- paste0(col, "Expected")
    diff_col <- paste0(col, "DifferenceFromExpected")
    if (expected_col %in% names(df)) {
      df[[diff_col]] <- df[[col]] - df[[expected_col]]
      diff_created <- diff_created + 1
    }
  }
  
  # Removed debug output - issue resolved
  
  # Rename difference columns to q codes for arrow plot
  diff_cols <- grep("DifferenceFromExpected$", names(df), value = TRUE)
  q_created <- 0
  for (col in diff_cols) {
    # Extract the base transition name
    base_name <- sub("DifferenceFromExpected$", "", col)
    # Get the corresponding q code
    q_code <- mapGenericToArrowNames(base_name)
    if (!is.na(q_code)) {
      # Rename the column to the q code
      names(df)[names(df) == col] <- q_code
      q_created <- q_created + 1
    }
  }
  
  # Removed debug output - issue resolved
  
  return(df)
}

# Local copy of getLabels function for independence
getLabels_local <- function(trait) {
  require(stringr)
  if (trait == "Final.polygyny") {
    return(c("Monogamy", "Polygyny"))
  } else if (grepl("coop", trait, ignore.case = TRUE)) {
    return(c("Non-Cooperative", "Cooperative"))
  } else if (str_detect(trait, "Kin")) {
    return(c("Non-kin", "Kin"))
  } else if (str_detect(trait, "Familial")) {
    return(c("Non-Familial Living", "Familial Living"))
  } else if (str_detect(trait, "Colonial")) {
    return(c("Non-Colonial", "Colonial"))
  } else if (str_detect(trait, "GroupsLargerThanPair")) {
    return(c("Asocial or pair", "Small or large groups"))
  } else if (str_detect(trait, "LongSocialBonds")) {
    return(c("Bonds last one season or less", "Multi-year bonds"))
  } else if (str_detect(trait, "MoreThanTwoCaretakers")) {
    return(c("Two or fewer caretakers", "More than two caretakers"))
  } else if (str_detect(trait, "TwoOrMoreCaretakers")) {
    return(c("Fewer than two caretakers", "Two or more caretakers"))
  } else if (str_detect(trait, "Asocial")) {
    return(c("Asocial", "Pair or group sociality"))
  } else if (str_detect(trait, "SeasonOrLonger")) {
    return(c("Bonds last less than one season", "Season or longer social bonds"))
  } else if (str_detect(trait, "LargestGroupSizes")) {
    return(c("Asocial, pair, or small groups", "Large groups"))
  } else if (str_detect(trait, "FemaleSong")) {
    return(c("Female Song Absent", "Female Song Present"))
  } else if (str_detect(trait, "TerritorialityWeakVsStrong")) {
    return(c("Weak Territoriality", "Strong Territoriality"))
  } else if (str_detect(trait, "Territory_12vs3")) {
    return(c("Territory 1 or 2", "Territory 3"))
  } else {
    return(c(paste(trait, "0"), paste(trait, "1")))
  }
}

# Validate transition counts
validateTransitionCounts <- function(countsdf) {
  # First check data types
  count_cols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01]$", names(countsdf), value = TRUE)
  
  # Check if columns are numeric
  non_numeric_cols <- character()
  for (col in count_cols) {
    if (!is.numeric(countsdf[[col]])) {
      non_numeric_cols <- c(non_numeric_cols, col)
    }
  }
  
  if (length(non_numeric_cols) > 0) {
    validationPrint("Non-numeric columns found in transition count data")
    message("Non-numeric transition columns: ", paste(non_numeric_cols, collapse=", "))
    for (col in non_numeric_cols[1:min(2, length(non_numeric_cols))]) {
      message("Column ", col, " class: ", class(countsdf[[col]]))
      message("First few values: ", paste(head(countsdf[[col]], 3), collapse=", "))
    }
    warning("Non-numeric data found in transition count columns")
    return(FALSE)
  }
  
  # Check that all counts are non-negative
  count_cols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01]$", names(countsdf), value = TRUE)
  
  for (col in count_cols) {
    if (any(countsdf[[col]] < 0, na.rm = TRUE)) {
      warning(paste("Negative counts found in column:", col))
      return(FALSE)
    }
  }
  
  # Check that totals match
  # If Total columns don't exist, calculate them for validation
  if ("Total_trait1_0to1" %in% names(countsdf)) {
    total_trait1_0to1 <- countsdf$Total_trait1_0to1
  } else {
    # Calculate total as sum of all trait1 0to1 transitions
    total_trait1_0to1 <- countsdf$trait1_0to1_in_trait2_0 + countsdf$trait1_0to1_in_trait2_1
  }
  
  sum_trait1_0to1 <- countsdf$trait1_0to1_in_trait2_0 + countsdf$trait1_0to1_in_trait2_1
  
  # This check is only meaningful if we have a separate Total column
  if ("Total_trait1_0to1" %in% names(countsdf)) {
    total_check1 <- all.equal(total_trait1_0to1, sum_trait1_0to1)
    
    if (!isTRUE(total_check1)) {
      validationPrint("trait1 0to1 transitions don't match sum of state-specific counts")
      message("First 5 Total_trait1_0to1 values: ", paste(head(total_trait1_0to1, 5), collapse=", "))
      message("First 5 sum of states: ", paste(head(sum_trait1_0to1, 5), collapse=", "))
      message("First 5 differences: ", paste(head(total_trait1_0to1 - sum_trait1_0to1, 5), collapse=", "))
      message("Maximum absolute difference: ", max(abs(total_trait1_0to1 - sum_trait1_0to1)))
      message("all.equal result: ", total_check1)
      warning("Total trait1 0to1 transitions don't match sum of state-specific counts")
      return(FALSE)
    }
  }
  
  # Check other total columns if they exist
  if ("Total_trait1_1to0" %in% names(countsdf)) {
    total_trait1_1to0 <- countsdf$Total_trait1_1to0
    sum_trait1_1to0 <- countsdf$trait1_1to0_in_trait2_0 + countsdf$trait1_1to0_in_trait2_1
    
    total_check2 <- all.equal(total_trait1_1to0, sum_trait1_1to0)
    
    if (!isTRUE(total_check2)) {
      validationPrint("trait1 1to0 transitions don't match sum of state-specific counts")
      debugPrint("- Total_trait1_1to0:", paste(head(total_trait1_1to0, 5), collapse=", "), "...")
      debugPrint("- Sum of states:", paste(head(sum_trait1_1to0, 5), collapse=", "), "...")
      debugPrint("- Differences:", paste(head(total_trait1_1to0 - sum_trait1_1to0, 5), collapse=", "), "...")
      debugPrint("- Max difference:", max(abs(total_trait1_1to0 - sum_trait1_1to0)))
      warning("Total trait1 1to0 transitions don't match sum of state-specific counts")
      return(FALSE)
    }
  }
  
  # Check proportion times sum to 1
  prop_check1 <- all(abs(countsdf$trait1_PropTime0 + countsdf$trait1_PropTime1 - 1) < 1e-6)
  prop_check2 <- all(abs(countsdf$trait2_PropTime0 + countsdf$trait2_PropTime1 - 1) < 1e-6)
  
  if (!prop_check1 || !prop_check2) {
    if (!prop_check1) {
      validationPrint("trait1 proportion times don't sum to 1")
      prop_sums <- countsdf$trait1_PropTime0 + countsdf$trait1_PropTime1
      debugPrint("- Prop sums:", paste(head(prop_sums, 5), collapse=", "), "...")
      debugPrint("- Max deviation from 1:", max(abs(prop_sums - 1)))
    }
    if (!prop_check2) {
      validationPrint("trait2 proportion times don't sum to 1")
      prop_sums <- countsdf$trait2_PropTime0 + countsdf$trait2_PropTime1
      debugPrint("- Prop sums:", paste(head(prop_sums, 5), collapse=", "), "...")
      debugPrint("- Max deviation from 1:", max(abs(prop_sums - 1)))
    }
    warning("Proportion times don't sum to 1")
    return(FALSE)
  }
  
  return(TRUE)
}

# Create RateRef dataframe for transition plot
createRateRef <- function(trait1_name, trait2_name) {
  # This creates the mapping between transition names and q rates
  # Using generic names internally
  RateRef <- data.frame(
    Transitions = c(
      "trait2_0to1_in_trait1_0", "trait1_0to1_in_trait2_0",
      "trait2_1to0_in_trait1_0", "trait1_0to1_in_trait2_1",
      "trait1_1to0_in_trait2_0", "trait2_0to1_in_trait1_1",
      "trait1_1to0_in_trait2_1", "trait2_1to0_in_trait1_1"
    ),
    qRate = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"),
    stringsAsFactors = FALSE
  )
  
  return(RateRef)
}