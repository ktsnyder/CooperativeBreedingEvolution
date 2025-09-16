# Debug script to diagnose file search issues

# Source the updated functions
source("create_phylopath_bias_robustness_figure_updated.R")

# Define trait sets
trait_sets <- list(
  weak_vs_strong = "FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET",
  terr_12vs3 = "FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET"
)

# Check what files are actually in the directories
cat("=== Checking actual directory contents ===\n\n")

# Check PhylopathJackknife
jackknife_dir <- "Outputs/PhylopathJackknife"
if (dir.exists(jackknife_dir)) {
  cat("Contents of", jackknife_dir, ":\n")
  subdirs <- list.dirs(jackknife_dir, full.names = FALSE, recursive = FALSE)
  for (subdir in subdirs) {
    if (subdir != "") {
      cat("  Subdirectory:", subdir, "\n")
      csv_files <- list.files(file.path(jackknife_dir, subdir), 
                              pattern = "\\.csv$", 
                              recursive = FALSE)
      if (length(csv_files) > 0) {
        cat("    CSV files:", paste(csv_files, collapse = ", "), "\n")
      }
    }
  }
} else {
  cat(jackknife_dir, "does not exist\n")
}

cat("\n")

# Check PhylopathDownsampled
downsampled_dir <- "Outputs/PhylopathDownsampled"
if (dir.exists(downsampled_dir)) {
  cat("Contents of", downsampled_dir, ":\n")
  subdirs <- list.dirs(downsampled_dir, full.names = FALSE, recursive = FALSE)
  for (subdir in subdirs) {
    if (subdir != "" && grepl("models", subdir)) {
      cat("  Subdirectory:", subdir, "\n")
      csv_files <- list.files(file.path(downsampled_dir, subdir), 
                              pattern = "detailed_models.*\\.csv$", 
                              recursive = FALSE)
      if (length(csv_files) > 0) {
        cat("    Detailed model files:", paste(csv_files[1:min(3, length(csv_files))], collapse = ", "))
        if (length(csv_files) > 3) cat(", ...")
        cat("\n")
      }
    }
  }
} else {
  cat(downsampled_dir, "does not exist\n")
}

cat("\n=== Testing file extraction with verbose mode ===\n")

# Test extraction for each trait set
for (trait_name in names(trait_sets)) {
  trait_set <- trait_sets[[trait_name]]
  
  cat("\n--- Testing", trait_name, "---\n")
  cat("Trait set:", trait_set, "\n\n")
  
  # Extract with verbose output
  results <- extract_phylopath_results(
    results_dir = "Outputs",
    trait_set = trait_set,
    verbose = TRUE  # This will show what directories and patterns it's searching
  )
  
  cat("\nFound", length(results), "bias corrections total\n")
  if (length(results) > 0) {
    cat("Bias corrections found:", paste(names(results), collapse = ", "), "\n")
  }
  
  # Check specifically for JackknifeSpecies
  if ("Jackknife by species" %in% names(results)) {
    cat("✓ JackknifeSpecies found\n")
  } else {
    cat("✗ JackknifeSpecies NOT found\n")
  }
}

cat("\n=== Checking for full dataset RDS files ===\n")

# Check for RDS files
rds_dir <- "Outputs/PhylopathPlots"
if (dir.exists(rds_dir)) {
  rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = FALSE)
  if (length(rds_files) > 0) {
    cat("RDS files in", rds_dir, ":\n")
    for (f in rds_files) {
      cat("  -", f, "\n")
    }
  } else {
    cat("No RDS files found in", rds_dir, "\n")
  }
} else {
  cat(rds_dir, "does not exist\n")
}

cat("\n=== Done ===\n")