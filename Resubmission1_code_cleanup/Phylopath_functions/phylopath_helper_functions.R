# Helper functions for phylopath bias robustness analysis

#' Convert trait set string to a label for filenames
#' Keeps the full trait set string for clarity
#'
#' @param trait_set Character string with full trait set name
#' @return Character string with label (same as input but with spaces replaced)
get_trait_set_label <- function(trait_set) {
  if (is.null(trait_set)) {
    return(NULL)
  }
  
  # Replace spaces with underscores for filename compatibility
  label <- gsub(" ", "_", trait_set)
  return(label)
}

#' Find or create full dataset result for a specific trait set
#'
#' @param trait_set Character string with full trait set name
#' @param run_if_missing Logical - whether to run phylopath if result doesn't exist
#' @return List with phylopath result or NULL if not found
find_or_create_full_dataset_result <- function(trait_set, run_if_missing = FALSE) {
  
  # Get label for filename (full trait set with underscores)
  trait_label <- get_trait_set_label(trait_set)
  
  # Define possible RDS filenames and locations
  rds_files <- c(
    # New format with full trait set in filename in Outputs/PhylopathFigures
    file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_result_", trait_label, ".rds")),
    # Also check without "_result" suffix
    file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_", trait_label, ".rds")),
    # New format with full trait set in filename in Outputs/PhylopathFigures, trait_set possibly with spaces
    file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_result_", trait_set, ".rds")),
    # Also check without "_result" suffix, trait_set possibly with spaces
    file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_", trait_set, ".rds")),
    # Legacy format in current directory
    "phylopath_full_dataset_result.rds",
    # Legacy format with trait set in current directory
    paste0("phylopath_full_dataset_", trait_label, "_result.rds")
  )
  
  # Check for existing RDS files
  for (rds_file in rds_files) {
    if (file.exists(rds_file)) {
      cat("Loading full dataset result from:", rds_file, "\n")
      result <- readRDS(rds_file)
      
      # Verify this is for the correct trait set if possible
      if (!is.null(result$trait_set)) {
        if (result$trait_set == trait_set) {
          return(result)
        } else {
          cat("Warning: RDS file is for different trait set:", result$trait_set, "\n")
        }
      } else if (rds_file == "phylopath_full_dataset_result.rds" && 
                 grepl("TerritorialityWeakVsStrong", trait_set)) {
        # Legacy file is likely for the default trait set
        return(result)
      }
    }
  }
  
  return(NULL)
}

#' Extract coefficient value from phylopath result for a specific path
#'
#' @param phylopath_result Result object from phylopath analysis
#' @param from_var Source variable name
#' @param to_var Target variable name
#' @param avg_method Method for averaging ("conditional" or "full")
#' @return Numeric coefficient value or NA if not found
extract_path_coefficient <- function(phylopath_result, from_var, to_var, avg_method = "conditional") {
  
  if (is.null(phylopath_result) || is.null(phylopath_result$result)) {
    return(NA)
  }
  
  # Get averaged result
  if (avg_method == "conditional") {
    avg_result <- phylopath::average(phylopath_result$result, avg_method = "conditional")
  } else {
    avg_result <- phylopath::average(phylopath_result$result, avg_method = "full")
  }
  
  # Extract coefficient
  if (!is.null(avg_result$coef) && from_var %in% rownames(avg_result$coef) && 
      to_var %in% colnames(avg_result$coef)) {
    return(avg_result$coef[from_var, to_var])
  } else {
    return(NA)
  }
}

#' Get variable name mappings for different trait sets
#'
#' @param trait_set Character string with full trait set name
#' @return List with variable name mappings
get_variable_mappings <- function(trait_set) {
  
  # Default mappings
  mappings <- list(
    CB = "HighConfidence_Coop",
    FS = "FemaleSong_Agg01",
    MASS = "logMass_AVONET"
  )
  
  # Add territory mapping based on trait set
  if (grepl("Territory_12vs3", trait_set)) {
    mappings$TERR <- "Territory_12vs3"
  } else if (grepl("TerritorialityWeakVsStrong", trait_set)) {
    mappings$TERR <- "TerritorialityWeakVsStrong"
  }
  
  return(mappings)
}