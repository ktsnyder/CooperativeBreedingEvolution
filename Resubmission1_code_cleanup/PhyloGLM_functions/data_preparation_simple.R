# Data Preparation for PhyloGLM Framework
# Minimal functions for matching data and tree

library(ape)

#' Prepare data for phyloglm analysis (simplified version)
#' 
#' @param data Original data frame
#' @param tree Original phylogenetic tree  
#' @param config Analysis configuration with $response, $predictors, $controls
#' @return List containing prepared data, tree, and metadata
prepare_analysis_data <- function(data, tree, config) {
  
  # Extract variables from config
  all_vars <- unique(c(config$response, config$predictors, config$controls))
  
  # Check variables exist
  missing_vars <- setdiff(all_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Ensure species column exists
  if (!"species" %in% colnames(data)) {
    if ("Species" %in% colnames(data)) {
      data$species <- data$Species
    } else if (rownames(data)[1] != "1") {
      data$species <- rownames(data)
    } else {
      stop("No species identifier found in data")
    }
  }
  
  # Remove duplicate species
  if (any(duplicated(data$species))) {
    data <- data[!duplicated(data$species), ]
  }
  
  # Filter to complete cases
  data_complete <- data[complete.cases(data[, all_vars]), ]
  
  # Match tree and data
  tree_species <- tree$tip.label
  data_species <- data_complete$species
  common_species <- intersect(tree_species, data_species)
  
  if (length(common_species) < 10) {
    stop(paste("Too few species for analysis:", length(common_species)))
  }
  
  # Prune tree and filter data
  tree_pruned <- drop.tip(tree, setdiff(tree_species, common_species))
  data_final <- data_complete[data_complete$species %in% common_species, ]
  
  # Ensure row names match tree tips
  rownames(data_final) <- data_final$species
  data_final <- data_final[tree_pruned$tip.label, ]
  
  return(list(
    data = data_final,
    tree = tree_pruned,
    n_species = nrow(data_final)
  ))
}

#' Check data quality (simplified version)
#' 
#' @param prepared_data Output from prepare_analysis_data
#' @return List with basic quality info
check_data_quality <- function(prepared_data) {
  list(
    n_species = prepared_data$n_species,
    is_acceptable = prepared_data$n_species >= 10
  )
}

#' Create data summary (simplified version) 
#' 
#' @param prepared_data Output from prepare_analysis_data
#' @return Basic data frame summary
create_data_summary <- function(prepared_data) {
  data.frame(
    n_species = prepared_data$n_species,
    stringsAsFactors = FALSE
  )
}