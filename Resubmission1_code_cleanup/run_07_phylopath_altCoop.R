### Run phylopath using alternative cooperative breeding classification methods/sources ----

# Supplemental Table 15; Extended Data Figure 6 ----
phylopath_output_dir <- file.path("Outputs","PhylopathFigures","AltCoops")
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
pdf(file.path(phylopath_output_dir,"all_AltCoop_enhanced_DAG_plots.pdf"), width = 42, height = 40)
do.call(grid.arrange, c(plotlist, ncol = 4))  # 3 columns, 5 rows
dev.off()
