
# Remove83HolarcticNoncoop
prefix_geo <- "Remove83HolarcticNoncoop_n500"
plots_geo <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = NULL, #result_geo_holarctic,
  downsampling_info = NULL, #geo_holarctic_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  detailed_models_input = "Outputs/PhylopathDownsampled/detailed_models_Remove83HolarcticNoncoop_500_2025-06-11.csv",
  model_frequencies_input = "Outputs/PhylopathDownsampled/model_frequencies_Remove83HolarcticNoncoop_500_2025-06-11.csv",
  output_prefix = prefix_geo,
  output_dir = phylopath_output_dir,
  save_png = TRUE
)

# Remove24TropicalCoop
prefix_trop <- "Remove24TropicalCoop_n500"
plots_trop <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = NULL, #result_geo_tropical,
  downsampling_info = NULL, #trop_coop_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  detailed_models_input = "Outputs/PhylopathDownsampled/detailed_models_Remove24TropicalCoop_500_2025-06-11.csv",
  model_frequencies_input = "Outputs/PhylopathDownsampled/model_frequencies_Remove24TropicalCoop_500_2025-06-11.csv",
  output_prefix = prefix_trop,
  output_dir = phylopath_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

# Remove15GlobalCoop
prefix_global <- "Remove15GlobalCoop_n500"
plots_global <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = NULL, #result_global_coop,
  downsampling_info = NULL, #global_coop_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  detailed_models_input = "Outputs/PhylopathDownsampled/detailed_models_Remove15GlobalCoop_500_2025-06-11.csv",
  model_frequencies_input = "Outputs/PhylopathDownsampled/model_frequencies_Remove15GlobalCoop_500_2025-06-11.csv",
  output_prefix = prefix_global,
  output_dir = phylopath_output_dir,
  save_png = TRUE
)


# Remove266StrongTerr
prefix_terr <- paste0("Remove266StrongTerr_n500")
plots_terr <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = NULL,
  downsampling_info = NULL, #266,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  detailed_models_input = "Outputs/PhylopathDownsampled/detailed_models_Remove266StrongTerr_500_2025-06-11.csv",
  model_frequencies_input = "Outputs/PhylopathDownsampled/model_frequencies_Remove266StrongTerr_500_2025-06-11.csv",
  output_prefix = prefix_terr,
  output_dir = phylopath_output_dir,
  save_png = TRUE, 
)


# Remove155Terr3
prefix_terr123 <- paste0("FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET Remove155Terr3_n500")
plots_terr123 <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = NULL, #result_terr,
  downsampling_info = NULL, #terr_downsample_info$report,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  detailed_models_input = "Outputs/PhylopathDownsampled/detailed_models_Remove155Terr3_500_2025-06-11.csv",
  model_frequencies_input = "Outputs/PhylopathDownsampled/model_frequencies_Remove155Terr3_500_2025-06-11.csv",
  output_prefix = prefix_terr123,
  output_dir = phylopath_output_dir,
  save_png = TRUE
)
