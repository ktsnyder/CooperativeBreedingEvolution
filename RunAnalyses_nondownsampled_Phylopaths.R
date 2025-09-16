## Run several non-downsampled, non-resampled Phylopath analyses

source("run_phylopath_fxns.R")
source("create_enhanced_DAG.R")

data <- read.csv("Data_R_2025-07-23.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

phylopath_output_dir <- "Outputs/PhylopathPlots/AltCoops"
if (!dir.exists(phylopath_output_dir)) dir.create(phylopath_output_dir, recursive = TRUE)

AltCoops = c("MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnCoop", "CockburnInferred", "BiagoliniCoop", "DaleCoop", "DowningCoop", "JetzCoopInclCockburn", "Griesser2017Coop", "CornwallisCoop")

plotlist = list()

for (i in 1:length(AltCoops)) {
  female_song_var = "FemaleSong_Agg01"
  coop_breeding_var = AltCoops[i]
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
  
  panel_letter = LETTERS[i]
  
  dag_file_path = file.path(phylopath_output_dir, paste(panel_letter, female_song_var, coop_breeding_var, territoriality_var, mass_var, "enhanced_dag.pdf"))
  plotlist[[i]] <- create_enhanced_dag(phylopath_result = result_phylopath, output_file = dag_file_path, 
                      title = paste(panel_letter, ")", coop_breeding_var, territoriality_var), 
                      female_song_var = female_song_var,
                      coop_breeding_var = coop_breeding_var,
                      territoriality_var = territoriality_var,
                      mass_var = mass_var)
}

library(gridExtra)
library(ggplot2)

# Save as PDF with all plots in a grid
pdf("all_AltCoop_enhanced_DAG_plots7.pdf", width = 42, height = 40)
do.call(grid.arrange, c(plotlist, ncol = 4))  # 3 columns, 5 rows
dev.off()



#### FemaleSong_Agg01_HighConfidence_Coop_Territory_12vs3_logMass_AVONET DAG ----
phylopath_output_dir <- "Outputs/PhylopathFigures/FemaleSong_Agg01_HighConfidence_Coop_Territory_12vs3_logMass_AVONET"

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

dag_file_path = file.path(phylopath_output_dir, paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, "enhanced_dag.pdf"))
print(dag_file_path)
create_enhanced_dag(phylopath_result = result_phylopath, output_file = dag_file_path, 
                                     title = paste(coop_breeding_var, territoriality_var), 
                                     female_song_var = female_song_var,
                                     coop_breeding_var = coop_breeding_var,
                                     territoriality_var = territoriality_var,
                                     mass_var = mass_var)


#### FemaleSong_Agg01_HighConfidence_Coop_TerritorialityWeakVsStrong_logMass_AVONET DAG ----
phylopath_output_dir <- "Outputs/PhylopathFigures/FemaleSong_Agg01_HighConfidence_Coop_TerritorialityWeakVsStrong_logMass_AVONET"

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

dag_file_path = file.path(phylopath_output_dir, paste(female_song_var, coop_breeding_var, territoriality_var, mass_var, "enhanced_dag.pdf"))
print(dag_file_path)
create_enhanced_dag(phylopath_result = result_phylopath, output_file = dag_file_path, 
                    title = paste(coop_breeding_var, territoriality_var), 
                    female_song_var = female_song_var,
                    coop_breeding_var = coop_breeding_var,
                    territoriality_var = territoriality_var,
                    mass_var = mass_var)





