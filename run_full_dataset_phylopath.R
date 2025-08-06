# run_full_dataset_phylopath.R
# Script to run phylopath analysis on the full dataset

# Load required libraries and functions
library(phylopath)
library(phytools)
library(dplyr)
library(ggplot2)

# Source the core phylopath functions
source("run_phylopath_fxns.R")

female_song_var = "FemaleSong_Agg01"
coop_breeding_var = "HighConfidence_Coop"
territoriality_var = "Territory_12vs3"
mass_var = "logMass_AVONET"

all_traits_phylopath_label <- paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, sep = "_")

# Load data
cat("Loading data...\n")
data <- read.csv("Data_R_2025-06-09.csv", stringsAsFactors = FALSE)
tree <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Run the full dataset phylopath analysis
cat("Running phylopath analysis on full dataset...\n")
result <- run_CB_FS_Terr_phylopath(
  dfIn = data,
  tree = tree,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  plots2pdf = TRUE,
  plots2png = TRUE,
  output_dir = file.path("Outputs","PhylopathPlots", all_traits_phylopath_label)
)

# Save the result object
cat("Saving phylopath result object...\n")
saveRDS(result, file.path("Outputs","PhylopathPlots",paste0("phylopath_full_dataset_result_", all_traits_phylopath_label, ".rds")))

# Create plots
cat("Creating plots...\n")

# Get conditional averaged model
cond_avg <- phylopath::average(result$result, cut_off = 2, avg_method = "conditional")

# Create DAG plot
phylopath_map_positions <- data.frame(
  name = c(female_song_var, coop_breeding_var, territoriality_var, mass_var),
  x = c(8, 2, 5, 5),
  y = c(9, 9, 5, 1),
  stringsAsFactors = FALSE
)

dag_plot <- plot(cond_avg,
                manual_layout = phylopath_map_positions,
                text_size = 4,
                box_x = 24, box_y = 18,
                edge_width = 2) +
  ggtitle("Full Dataset Phylogenetic Path Analysis") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.margin = margin(30, 30, 30, 30)) +
  coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)

ggsave(file.path("Outputs","PhylopathPlots",paste0("phylopath_full_dataset_dag_", all_traits_phylopath_label, ".png")), dag_plot, width = 8, height = 8, dpi = 300)

# Extract path coefficients for reference
cat("\nConditional averaged path coefficients:\n")
coef_matrix <- cond_avg$coef
print(coef_matrix)

# Save coefficients to CSV for reference
coef_df <- as.data.frame(as.table(coef_matrix))
names(coef_df) <- c("From", "To", "Coefficient")
coef_df <- coef_df[!is.na(coef_df$Coefficient), ]
write.csv(coef_df, file.path("Outputs","PhylopathPlots",paste0("phylopath_full_dataset_coefficients_", all_traits_phylopath_label,".csv")), row.names = FALSE)

cat("\nFull dataset phylopath analysis complete!\n")
cat("Results saved to:\n")
cat("- phylopath_full_dataset_result.rds (R object)\n")
cat("- phylopath_full_dataset_dag.png (DAG plot)\n")
cat("- phylopath_full_dataset_coefficients.csv (coefficients)\n")

# Return the result
result