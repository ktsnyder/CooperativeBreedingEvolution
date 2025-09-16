### Simmap overlap runner
### Kate Snyder
### 7/1/2025

source("simmap_overlap_runner_helpers.R")

#### TerrWeakStrong/FemaleSong ----
starttime = Sys.time()
result <- runSimmapOverlapAnalysis(
  trait1 = "TerritorialityWeakVsStrong",
  trait2 = "FemaleSong_Agg01",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = TRUE,
  plot_transitions = TRUE,
  save_outputs = TRUE,
  other_label = "BothDummy"
)
endtime = Sys.time()
endtime-starttime

result <- runSimmapOverlapAnalysis(
  trait1 = "TerritorialityWeakVsStrong",
  trait2 = "FemaleSong_Agg01",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = FALSE,  # Changed to FALSE
  plot_transitions = FALSE,       # Changed to FALSE
  save_outputs = TRUE,
  use_pregenerated_simmaps = TRUE,  # Added - tells it to use existingfiles
  trait1_real_simmaps_file = "Simmap Overlap Outputs/TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_500_BothDummy_20250701_015401/simmaps/TerritorialityWeakVsStrong_FemaleSong_Agg01_simmaps_REAL_500_BothDummy.rds",
  trait2_real_simmaps_file = "Simmap Overlap Outputs/TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_500_BothDummy_20250701_015401/simmaps/FemaleSong_Agg01_TerritorialityWeakVsStrong_simmaps_REAL_500_BothDummy.rds",
  trait1_dummy_simmaps_file = "Simmap Overlap Outputs/TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_500_BothDummy_20250701_015401/simmaps/TerritorialityWeakVsStrong_FemaleSong_Agg01_simmaps_REAL_500_BothDummy.rds",
  trait2_dummy_simmaps_file = "Simmap Overlap Outputs/TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_500_BothDummy_20250701_015401/simmaps/FemaleSong_Agg01_TerritorialityWeakVsStrong_simmaps_DUMMY_500_BothDummy.rds",
  other_label = "OnlyFSDummy_pregeneratedSimmaps_NoTransitions"  # Changed label to distinguish
)




#### TerrWeakStrong/Coop ----
result <- runSimmapOverlapAnalysis(
  trait1 = "TerritorialityWeakVsStrong",
  trait2 = "HighConfidence_Coop",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = TRUE,
  plot_transitions = TRUE,
  save_outputs = FALSE,
  other_label = "BothDummy"
)


result <- runSimmapOverlapAnalysis(
  trait1 = "TerritorialityWeakVsStrong",
  trait2 = "HighConfidence_Coop",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = FALSE,  # Changed to FALSE
  plot_transitions = FALSE,       # Changed to FALSE
  save_outputs = TRUE,
  use_pregenerated_simmaps = TRUE,  # Added - tells it to use existingfiles
  trait1_real_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/TerritorialityWeakVsStrong_HighConfidence_Coop_simmaps_REAL_500_BothDummy.rds",
  trait2_real_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/HighConfidence_Coop_TerritorialityWeakVsStrong_simmaps_REAL_500_BothDummy.rds",
  trait1_dummy_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/TerritorialityWeakVsStrong_HighConfidence_Coop_simmaps_REAL_500_BothDummy.rds",
  trait2_dummy_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/HighConfidence_Coop_TerritorialityWeakVsStrong_simmaps_DUMMY_500_BothDummy.rds",
  other_label = "OnlyCBDummy_pregeneratedSimmaps_NoTransitions"  # Changed label to distinguish
)
Sys.time()

# Because save_outputs was accidentally set to FALSE before, redo BothDummy with TRUE using pre-generated simmaps - has not been run yet 7/1 4:42pm
result <- runSimmapOverlapAnalysis(
  trait1 = "TerritorialityWeakVsStrong",
  trait2 = "HighConfidence_Coop",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = TRUE,
  plot_transitions = TRUE,
  save_outputs = TRUE,
  use_pregenerated_simmaps = TRUE,  # Added - tells it to use existingfiles
  trait1_real_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/TerritorialityWeakVsStrong_HighConfidence_Coop_simmaps_REAL_500_BothDummy.rds",
  trait2_real_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/HighConfidence_Coop_TerritorialityWeakVsStrong_simmaps_REAL_500_BothDummy.rds",
  trait1_dummy_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/TerritorialityWeakVsStrong_HighConfidence_Coop_simmaps_DUMMY_500_BothDummy.rds",
  trait2_dummy_simmaps_file = "Simmap_Overlap_Outputs/TerritorialityWeakVsStrong_vs_HighConfidence_Coop_500_500_BothDummy_20250701_030906/simmaps/HighConfidence_Coop_TerritorialityWeakVsStrong_simmaps_DUMMY_500_BothDummy.rds",
  other_label = "BothDummy"
)


#### Terr12vs3/Coop ----
result <- runSimmapOverlapAnalysis(
  trait1 = "Territory_12vs3",
  trait2 = "HighConfidence_Coop",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = TRUE,
  plot_transitions = TRUE,
  save_outputs = TRUE,
  other_label = ""
)

result <- runSimmapOverlapAnalysis(
  trait1 = "Territory_12vs3",
  trait2 = "FemaleSong_Agg01",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-06-09.csv",
  calculate_transitions = TRUE,
  plot_transitions = TRUE,
  save_outputs = TRUE,
  other_label = ""
)


#### Terr12vs3 with coop & FS - no randomization of Terr simmaps 8/5/2025 ----
result <- runSimmapOverlapAnalysis(
  trait1 = "Territory_12vs3",
  trait2 = "HighConfidence_Coop",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-07-23.csv",
  calculate_transitions = FALSE,
  plot_transitions = FALSE,
  save_outputs = TRUE,
  use_pregenerated_simmaps = TRUE,  # Added - tells it to use existing files
  trait1_real_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_HighConfidence_Coop_500_500_20250701_121109/simmaps/Territory_12vs3_HighConfidence_Coop_simmaps_REAL_500.rds",
  trait2_real_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_HighConfidence_Coop_500_500_20250701_121109/simmaps/HighConfidence_Coop_Territory_12vs3_simmaps_REAL_500.rds",
  trait1_dummy_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_HighConfidence_Coop_500_500_20250701_121109/simmaps/Territory_12vs3_HighConfidence_Coop_simmaps_REAL_500.rds",
  trait2_dummy_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_HighConfidence_Coop_500_500_20250701_121109/simmaps/HighConfidence_Coop_Territory_12vs3_simmaps_DUMMY_500.rds",
  other_label = "OnlyCBDummy_pregeneratedSimmaps_NoTransitions" 
)

result <- runSimmapOverlapAnalysis(
  trait1 = "Territory_12vs3",
  trait2 = "FemaleSong_Agg01",
  nsims_real = 500,
  nsims_dummy = 500,
  tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  data_file = "Data_R_2025-07-23.csv",
  calculate_transitions = FALSE,
  plot_transitions = FALSE,
  save_outputs = TRUE,
  use_pregenerated_simmaps = TRUE,  # Added - tells it to use existing files
  trait1_real_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_FemaleSong_Agg01_500_500_20250701_203938/simmaps/Territory_12vs3_FemaleSong_Agg01_simmaps_REAL_500.rds",
  trait2_real_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_FemaleSong_Agg01_500_500_20250701_203938/simmaps/FemaleSong_Agg01_Territory_12vs3_simmaps_REAL_500.rds",
  trait1_dummy_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_FemaleSong_Agg01_500_500_20250701_203938/simmaps/Territory_12vs3_FemaleSong_Agg01_simmaps_REAL_500.rds",
  trait2_dummy_simmaps_file = "Simmap_Overlap_Outputs/Territory_12vs3_vs_FemaleSong_Agg01_500_500_20250701_203938/simmaps/FemaleSong_Agg01_Territory_12vs3_simmaps_DUMMY_500.rds",
  other_label = "OnlyFSDummy_pregeneratedSimmaps_NoTransitions" 
)
