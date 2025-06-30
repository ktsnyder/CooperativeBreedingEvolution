## Script to run CharacterSimmaps_modified and calcHuel_corrected on all trait pairs
## Traits: HighConfidence_Coop, FemaleSong_Agg01, TerritorialityWeakVsStrong
## 150 real and dummy sims for each pair
## Kate Snyder / Claude
## 2025-06-27

# Load required libraries
library(phytools)
library(tidyverse)

# Source the functions
source("CharacterSimmaps_modified.R")
source("calcHuel_corrected.R")

# Set parameters
nsims <- 10
traits <- c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong")
traits <- c("FemaleSong_Agg01", "Territory_12vs3")
traits <- c("TerritorialityWeakVsStrong", "FemaleSong_Agg01")

# Load tree and data (adjust paths as needed)
tree_file <- "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
data_file <- "Data_R_2025-06-09.csv"

# Read tree and data
tree <- read.nexus(tree_file)
df <- read.csv(data_file)

# Create all pairwise combinations
trait_pairs <- combn(traits, 2, simplify = FALSE)

# Create summary dataframe to store results
all_results <- data.frame()

# Process each trait pair
for (pair_idx in 1:length(trait_pairs)) {
  trait_pair <- trait_pairs[[pair_idx]]
  trait1 <- trait_pair[1]
  trait2 <- trait_pair[2]
  
  cat("\n=====================================\n")
  cat("Processing trait pair:", trait1, "vs", trait2, "\n")
  cat("=====================================\n\n")
  
  # Generate timestamp for this analysis
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create output directory for this pair
  pair_output_dir <- paste0("Simmap_Overlap_Outputs_", trait1, "_vs_", trait2, "_", nsims, "_", timestamp)
  if (!dir.exists(pair_output_dir)) {
    dir.create(pair_output_dir)
  }
  
  ditree = multi2di(tree)
  ditree$edge.length[which(ditree$edge.length == 0)] <- 0.0000000000001
  # Run CharacterSimtree_info# Run CharacterSimmaps_modified for REAL data
  cat("Generating REAL data simmaps...\n")
  dfout_real <- CharacterSimmaps_modified(
    columns = c(trait1, trait2),
    df = df,
    tree = ditree,
    dummy = FALSE,
    nsims = nsims,
    treelabel = "Re-refactored",
    datalabel = paste(trait1, trait2),
    plotSampleSimmaps = TRUE,
    calculate_transitions = FALSE,
    save_simmaps = TRUE,
    output_dir = pair_output_dir
  )
  
  # Run CharacterSimmaps_modified for DUMMY data
  cat("\nGenerating DUMMY data simmaps...\n")
  dfout_dummy <- CharacterSimmaps_modified(
    columns = c(trait1, trait2),
    df = df,
    tree = ditree,
    dummy = TRUE,
    nsims = nsims,
    treelabel = "Re-refactored",
    datalabel = paste(trait1, trait2),
    dummyMethod = "makeSimmap",
    plotSampleSimmaps = TRUE,
    calculate_transitions = FALSE,
    save_simmaps = TRUE,
    output_dir = pair_output_dir
  )
  
  # Run calcHuel_corrected
  cat("\nRunning calcHuel_corrected analysis...\n")
  huel_results <- calcHuel_corrected(
    dfout = dfout_real,
    dfDummy = dfout_dummy,
    nsims_real = nsims,
    nsims_dummy = nsims,
    otherlabel = "Re-refactored",
    newplot = TRUE,
    plot_ggplots_pdf = TRUE,
    trait1_name = trait1,
    trait2_name = trait2
  )
  
  # Save the plots
  if (!is.null(huel_results$plots)) {
    # Save combined plot
    library(gridExtra)
    pdf_filename <- file.path(pair_output_dir, 
                             paste0("Simmap_Overlap_Analysis_", trait1, "_vs_", trait2, "_", nsims, "_", timestamp, ".pdf"))
    pdf(pdf_filename, width = 8, height = 14)
    
    if (length(huel_results$plots) >= 3) {
      grid.arrange(huel_results$plots$p1, huel_results$plots$p2, huel_results$plots$p3, 
                   ncol = 1, heights = c(1, 1, 1.2))
    } else {
      grid.arrange(huel_results$plots$p1, huel_results$plots$p2, ncol = 1)
    }
    dev.off()
    cat("Saved analysis plots to:", pdf_filename, "\n")
  }
  
  # Store summary results
  if (!is.null(huel_results$mediansRow)) {
    all_results <- rbind(all_results, huel_results$mediansRow)
  }
  
  # Save intermediate results
  saveRDS(list(
    dfout_real = dfout_real,
    dfout_dummy = dfout_dummy,
    huel_results = huel_results,
    trait1 = trait1,
    trait2 = trait2,
    timestamp = timestamp
  ), file.path(pair_output_dir, paste0("analysis_results_", trait1, "_vs_", trait2, "_", timestamp, ".rds")))
  
  cat("\nCompleted analysis for", trait1, "vs", trait2, "\n")
}

# Save summary results
summary_filename <- paste0("trait_pairs_summary_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
write.csv(all_results, summary_filename, row.names = FALSE)
cat("\n\nAll analyses complete. Summary results saved to:", summary_filename, "\n")

# Print summary table
cat("\nSummary of p-values:\n")
print(all_results[, c("trait1", "trait2", "pval", "D_real")])


#### Trait pairs with pregenerated simmaps ----

source("claude_code_sessions/CharacterSimmaps_modified.R")

trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
timestamp = ""
Realresult <- CharacterSimmaps_modified(
  columns = c("HighConfidence_Coop", "FemaleSong_Agg01"),
  df = df,
  tree = tree,
  dummy = FALSE,
  nsims = 100,  # Will be overridden by actual number in RDS
  treelabel = "PregeneratedSimmaps",
  trait1_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps HighConfidence_Coop CONSISTENT_RealData 2025-06-25_2007 .rds",
  trait2_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps FemaleSong_Agg01 CONSISTENT_RealData 2025-06-25_1959 .rds",
  calculate_transitions = TRUE,
  save_simmaps = FALSE,
  output_dir = "output"
)

Dummyresult <- CharacterSimmaps_modified(
  columns = c("HighConfidence_Coop", "FemaleSong_Agg01"),
  df = df,
  tree = tree,
  dummy = FALSE,
  nsims = 100,  # Will be overridden by actual number in RDS
  treelabel = "PregeneratedSimmaps",
  trait1_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps HighConfidence_Coop CONSISTENT_DummyData 2025-06-25 .rds",
  trait2_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps FemaleSong_Agg01 CONSISTENT_DummyData 2025-06-25 .rds",
  calculate_transitions = TRUE,
  save_simmaps = FALSE,
  output_dir = "output"
)

# Run calcHuel_corrected
cat("\nRunning calcHuel_corrected analysis...\n")
huel_results <- calcHuel_corrected(
  dfout = Realresult,
  dfDummy = Dummyresult,
  nsims_real = 100,
  nsims_dummy = 100,
  otherlabel = "Pre-generatedSimmaps",
  newplot = TRUE,
  plot_ggplots_pdf = TRUE,
  trait1_name = trait1,
  trait2_name = trait2
)

# Save the plots
if (!is.null(huel_results$plots)) {
  # Save combined plot
  library(gridExtra)
  
  # Create output directory for pre-generated simmaps
  output_dir_for_plots <- paste0("Simmap_Overlap_Outputs_", trait1, "_vs_", trait2, "_pregenerated_100")
  if (!dir.exists(output_dir_for_plots)) {
    dir.create(output_dir_for_plots)
  }
  
  pdf_filename <- file.path(output_dir_for_plots, 
                            paste0("Simmap_Overlap_Analysis_", trait1, "_vs_", trait2, "_", 100, "sims_", timestamp, ".pdf"))
  pdf(pdf_filename, width = 8, height = 14)
  
  if (length(huel_results$plots) >= 3) {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, huel_results$plots$p3, 
                 ncol = 1, heights = c(1, 1, 1.2))
  } else {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, ncol = 1)
  }
  dev.off()
  cat("Saved analysis plots to:", pdf_filename, "\n")
}

# Store summary results
if (!is.null(huel_results$mediansRow)) {
  all_results <- rbind(all_results, huel_results$mediansRow)
}

# Save intermediate results
saveRDS(list(
  dfout_real = Realresults,
  dfout_dummy = Dummyresults,
  huel_results = huel_results,
  trait1 = trait1,
  trait2 = trait2,
  timestamp = timestamp
), file.path(pair_output_dir, paste0("analysis_results_", trait1, "_vs_", trait2, "_", timestamp, ".rds")))

cat("\nCompleted analysis for", trait1, "vs", trait2, "\n")



## using simmap generated from full dataset 
trait1 = "TerritorialityWeakVsStrong"
trait2 = "FemaleSong_Agg01"
timestamp = ""
Realresult <- CharacterSimmaps_modified(
  columns = c("TerritorialityWeakVsStrong", "FemaleSong_Agg01"),
  df = df,
  tree = tree,
  dummy = FALSE,
  nsims = 100,  # Will be overridden by actual number in RDS
  treelabel = "PregeneratedSimmaps",
  trait1_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps TerritorialityWeakVsStrong CONSISTENT_RealData 2025-06-25_1957 .rds",
  trait2_real_multisimmapRDS = "Simmap Overlap Outputs/100 multisimmaps FemaleSong_Agg01 CONSISTENT_RealData 2025-06-25_1959 .rds",
  calculate_transitions = FALSE,
  save_simmaps = FALSE,
  output_dir = "output_reorderedTrait1Trait2"
)

Dummyresult <- CharacterSimmaps_modified(
  columns = c("TerritorialityWeakVsStrong", "FemaleSong_Agg01"),
  df = df,
  tree = tree,
  dummy = TRUE,
  nsims = 100,  # Will be overridden by actual number in RDS
  treelabel = "PregeneratedDummySimmaps",
  trait1_dummy_multisimmapRDS = "",
  trait2_dummy_multisimmapRDS = "",
  calculate_transitions = FALSE,
  save_simmaps = FALSE,
  output_dir = "output_reorderedTrait1Trait2"
)

# Run calcHuel_corrected
cat("\nRunning calcHuel_corrected analysis...\n")
huel_results <- calcHuel_corrected(
  dfout = Realresult,
  dfDummy = Dummyresult,
  nsims_real = 100,
  nsims_dummy = 100,
  otherlabel = "Pre-generatedSimmaps",
  newplot = TRUE,
  plot_ggplots_pdf = TRUE,
  trait1_name = trait1,
  trait2_name = trait2
)

# Save the plots
if (!is.null(huel_results$plots)) {
  # Save combined plot
  library(gridExtra)
  
  # Create output directory for pre-generated simmaps
  output_dir_for_plots <- paste0("Simmap_Overlap_Outputs_", trait1, "_vs_", trait2, "_pregenerated_100")
  if (!dir.exists(output_dir_for_plots)) {
    dir.create(output_dir_for_plots)
  }
  
  pdf_filename <- file.path(output_dir_for_plots, 
                            paste0("Simmap_Overlap_Analysis_", trait1, "_vs_", trait2, "_", 100, "sims_", timestamp, ".pdf"))
  pdf(pdf_filename, width = 8, height = 14)
  
  if (length(huel_results$plots) >= 3) {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, huel_results$plots$p3, 
                 ncol = 1, heights = c(1, 1, 1.2))
  } else {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, ncol = 1)
  }
  dev.off()
  cat("Saved analysis plots to:", pdf_filename, "\n")
}



## using simmaps generated from subsetted trees --
ditree = multi2di(tree)
ditree$edge.length[which(ditree$edge.length == 0)] <- 0.000000000000001
trait1 = "TerritorialityWeakVsStrong"
trait2 = "FemaleSong_Agg01"
timestamp = "FSdummyOnlyMaybeFixed"
Realresult <- CharacterSimmaps_modified(
  columns = c("TerritorialityWeakVsStrong", "FemaleSong_Agg01"),
  df = df,
  tree = ditree,
  dummy = FALSE,
  nsims = 500,  # Will be overridden by actual number in RDS
  treelabel = "SimmapsFromPrevRun",
  trait1_real_multisimmapRDS = "Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/SIMMAPS_FOR_TerritorialityWeakVsStrong_whenPairedWith_FemaleSong_Agg01_500sims_REAL_drop_tips_then_simmap_2025-06-27_1845.rds",
  trait2_real_multisimmapRDS = "Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/SIMMAPS_FOR_FemaleSong_Agg01_whenPairedWith_TerritorialityWeakVsStrong_500sims_REAL_drop_tips_then_simmap_2025-06-27_1845.rds",
  calculate_transitions = FALSE,
  save_simmaps = FALSE,
  output_dir = "output"
)



Dummyresult <- CharacterSimmaps_modified(
  columns = c("TerritorialityWeakVsStrong", "FemaleSong_Agg01"),
  df = df,
  tree = tree,
  dummy = TRUE,
  nsims = 500,  # Will be overridden by actual number in RDS
  treelabel = "OnlyFSDummy",
  trait1_dummy_multisimmapRDS = "Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/SIMMAPS_FOR_TerritorialityWeakVsStrong_whenPairedWith_FemaleSong_Agg01_500sims_REAL_drop_tips_then_simmap_2025-06-27_1845.rds",
  trait2_dummy_multisimmapRDS = "Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/SIMMAPS_FOR_FemaleSong_Agg01_whenPairedWith_TerritorialityWeakVsStrong_500sims_DUMMY_makeSimmap_drop_tips_then_simmap_2025-06-27_1917.rds",
  calculate_transitions = FALSE,
  save_simmaps = FALSE,
  output_dir = "output"
)

# Run calcHuel_corrected
cat("\nRunning calcHuel_corrected analysis...\n")
huel_results <- calcHuel_corrected(
  dfout = Realresult,
  dfDummy = Dummyresult,
  nsims_real = 500,
  nsims_dummy = 500,
  otherlabel = "OnlyFSDummySimmaps",
  newplot = TRUE,
  plot_ggplots_pdf = TRUE,
  trait1_name = trait1,
  trait2_name = trait2
)

# Save the plots
if (!is.null(huel_results$plots)) {
  # Save combined plot
  library(gridExtra)
  
  # Create output directory for this analysis
  output_dir_for_plots <- paste0("Simmap_Overlap_Outputs_", trait1, "_vs_", trait2, "_pregenerated_500")
  if (!dir.exists(output_dir_for_plots)) {
    dir.create(output_dir_for_plots)
  }
  
  pdf_filename <- file.path(output_dir_for_plots, 
                            paste0("Simmap_Overlap_Analysis_", trait1, "_vs_", trait2, "_", 500, "sims_", timestamp, ".pdf"))
  pdf(pdf_filename, width = 8, height = 14)
  
  if (length(huel_results$plots) >= 3) {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, huel_results$plots$p3, 
                 ncol = 1, heights = c(1, 1, 1.2))
  } else {
    grid.arrange(huel_results$plots$p1, huel_results$plots$p2, ncol = 1)
  }
  dev.off()
  cat("Saved analysis plots to:", pdf_filename, "\n")
}


