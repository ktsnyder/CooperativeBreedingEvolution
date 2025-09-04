## Run several non-downsampled, non-resampled Phylopath analyses
# Generates Figure 4A; Extended Data Figures 7A

source("Phylopath_functions/run_phylopath_fxns.R")
source("Phylopath_functions/create_enhanced_DAG.R")

data <- read.csv("Data_R.csv")
tree <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")


#### Extended Data Figure 7A - FemaleSong_Agg01_HighConfidence_Coop_Territory_12vs3_logMass_AVONET DAG ----
phylopath_output_dir <- file.path("Outputs", "PhylopathFigures", "FemaleSong_Agg01_HighConfidence_Coop_Territory_12vs3_logMass_AVONET")

female_song_var = "FemaleSong_Agg01"
coop_breeding_var = "HighConfidence_Coop"
territoriality_var = "Territory_12vs3"
mass_var = "logMass_AVONET"

result_phylopath <- run_CB_FS_Terr_phylopath(
  dfIn = data,
  tree = tree,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  plots2pdf = TRUE, 
  output_dir = phylopath_output_dir
)

trait_label = paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, sep = "_")
RDSpath <- file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_result_", trait_label, ".rds"))
saveRDS(result_phylopath, RDSpath)

dag_file_path = file.path(phylopath_output_dir, paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, "enhanced_dag.pdf"))
print(dag_file_path)
create_enhanced_dag(phylopath_result = result_phylopath, output_file = dag_file_path, 
                                     title = paste(coop_breeding_var, territoriality_var), 
                                     female_song_var = female_song_var,
                                     coop_breeding_var = coop_breeding_var,
                                     territoriality_var = territoriality_var,
                                     mass_var = mass_var)


#### Figure 4A - FemaleSong_Agg01_HighConfidence_Coop_TerritorialityWeakVsStrong_logMass_AVONET DAG ----
phylopath_output_dir <- file.path("Outputs", "PhylopathFigures", "FemaleSong_Agg01_HighConfidence_Coop_TerritorialityWeakVsStrong_logMass_AVONET")

female_song_var = "FemaleSong_Agg01"
coop_breeding_var = "HighConfidence_Coop"
territoriality_var = "TerritorialityWeakVsStrong"
mass_var = "logMass_AVONET"

result_phylopath <- run_CB_FS_Terr_phylopath(
  dfIn = data,
  tree = tree,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  plots2pdf = TRUE, 
  output_dir = phylopath_output_dir
)

trait_label = paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, sep = "_")
RDSpath <- file.path("Outputs", "PhylopathFigures", paste0("phylopath_full_dataset_result_", trait_label, ".rds"))
saveRDS(result_phylopath, RDSpath)

dag_file_path = file.path(phylopath_output_dir, paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, "enhanced_dag.pdf"))
print(dag_file_path)
create_enhanced_dag(phylopath_result = result_phylopath, output_file = dag_file_path, 
                    title = paste(coop_breeding_var, territoriality_var), 
                    female_song_var = female_song_var,
                    coop_breeding_var = coop_breeding_var,
                    territoriality_var = territoriality_var,
                    mass_var = mass_var)


