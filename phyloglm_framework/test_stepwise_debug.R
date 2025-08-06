# Debug script for stepwise expansion

library(phylolm)
library(ape)

# Load data
cat("Loading data...\n")
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Get the base result
fs_result <- all_results[["FS_vs_CB_TerrWS_Mass"]]
cat("Base model AIC:", fs_result$comparison$comparison$AIC[1], "\n")

# Check the prepared data from the original analysis
cat("\nOriginal analysis had", nrow(fs_result$prepared_data$data), "species\n")

# Get the BEST model, not the first model
best_model_name <- fs_result$comparison$comparison$Model[1]
cat("Best model name:", best_model_name, "\n")
base_formula <- formula(fs_result$models$models[[best_model_name]])
cat("\nBase formula:", deparse(base_formula), "\n")

# Get variables needed
vars_needed <- all.vars(base_formula)
cat("\nVariables needed:", paste(vars_needed, collapse=", "), "\n")

# Check if all variables exist
missing <- setdiff(vars_needed, names(data))
if (length(missing) > 0) {
  cat("Missing variables:", paste(missing, collapse=", "), "\n")
}

# Prepare data subset
data_subset <- data[, c("species", vars_needed)]
complete_rows <- complete.cases(data_subset)
data_subset <- data_subset[complete_rows, ]
cat("\nData after removing NAs:", nrow(data_subset), "species\n")

# Match with tree
shared_species <- intersect(data_subset$species, tree$tip.label)
cat("Shared species between data and tree:", length(shared_species), "\n")

data_subset <- data_subset[data_subset$species %in% shared_species, ]
tree_subset <- keep.tip(tree, shared_species)

# Order data to match tree
data_subset <- data_subset[match(tree_subset$tip.label, data_subset$species), ]
rownames(data_subset) <- data_subset$species

# Check alignment
cat("\nData rows match tree tips:", all(rownames(data_subset) == tree_subset$tip.label), "\n")

# Try fitting base model
cat("\nTrying to fit base model...\n")
tryCatch({
  base_model <- phyloglm(
    formula = base_formula,
    data = data_subset,
    phy = tree_subset,
    method = "logistic_MPLE",
    btol = 30,
    log.alpha.bound = 4
  )
  cat("Success! AIC =", -2 * base_model$logLik + 2 * base_model$d, "\n")
  
  # Now try adding Territory as numeric
  cat("\nTrying to add Territory as numeric...\n")
  
  # First check if Territory exists and what it contains
  if ("Territory" %in% names(data)) {
    # Get Territory data for matched species
    territory_data <- data[data$species %in% rownames(data_subset), c("species", "Territory")]
    territory_data <- territory_data[match(rownames(data_subset), territory_data$species), ]
    
    cat("Territory values:\n")
    print(table(territory_data$Territory, useNA = "always"))
    
    data_subset$Territory_num <- as.numeric(territory_data$Territory)
  } else {
    cat("Territory column not found in data!\n")
    cat("Available columns:", paste(head(names(data), 20), collapse=", "), "...\n")
  }
  
  new_formula <- update(base_formula, ~ . + Territory_num)
  cat("New formula:", deparse(new_formula), "\n")
  
  new_model <- phyloglm(
    formula = new_formula,
    data = data_subset,
    phy = tree_subset,
    method = "logistic_MPLE",
    btol = 30,
    log.alpha.bound = 4
  )
  
  new_aic <- -2 * new_model$logLik + 2 * new_model$d
  cat("New model AIC:", new_aic, "\n")
  cat("AIC improvement:", (-2 * base_model$logLik + 2 * base_model$d) - new_aic, "\n")
  
}, error = function(e) {
  cat("Error:", e$message, "\n")
  
  # Additional debugging
  cat("\nDebugging info:\n")
  cat("Tree tips (first 5):", head(tree_subset$tip.label), "\n")
  cat("Data species (first 5):", head(rownames(data_subset)), "\n")
  cat("Tree is ultrametric:", is.ultrametric(tree_subset), "\n")
  cat("Tree has branch lengths:", !is.null(tree_subset$edge.length), "\n")
})