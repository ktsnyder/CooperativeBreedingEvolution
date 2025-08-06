
jackpath = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathJackknife/detailed_models_JackknifeSpecies_n875_2025-07-08.csv")

data_file = "Data_R_2025-06-09.csv"
df = read.csv(data_file)

jackpath = merge(jackpath, df[,c("species", "HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong", "Territory_12vs3", "logMass_AVONET", "GeographicRegion_Jetz", "logMaleFemalePlumageDiffAbs", "PercentAbsLogWingDimorphism")], by.x = "RemovedSpecies", by.y = "species", all.x = T)

write.csv(jackpath, "/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathJackknife/JackknifedSpeciesData_detailed_models_JackknifeSpecies_n875_2025-07-08.csv", row.names = F)


modelfreqs = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathJackknife/model_frequencies_JackknifeSpecies_n875_2025-07-08.csv")



batch_results <- run_phyloglm_batch(
  analysis_configs = configs,
  data = read.csv("Data_R_2025-06-09.csv"),
  tree = ape::read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"),
  output_dir = "Outputs/PhyloglmResults",
  bootstrap_n = 10,
  save_intermediate = TRUE  # This will skip the plotting step
)