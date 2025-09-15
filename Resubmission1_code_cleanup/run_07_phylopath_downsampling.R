#### Run phylopath analyses with downsampling - binary traits ----
# Performs analyses used in RunAnalyses_generate_bias_figures_bothtraits.R to make Figure 4B-C, Extended Data Figure 7B-D

# Must first run "run_05_bias_analyses.R"

female_song_var = "FemaleSong_Agg01"
coop_breeding_var = "HighConfidence_Coop"
if (!exists("territoriality_var")) {
  territoriality_var =  "TerritorialityWeakVsStrong" 
  print("Defaulting to TerritorialityWeakVsStrong as the territoriality variable in phylopath models. To run using year-round territoriality, set territoriality_var = Territory_12vs3.")
}
#territoriality_var =  "Territory_12vs3" # uncomment to test using year-round territoriality as the territoriality feature
mass_var = "logMass_AVONET"

if(!exists("n_iterations")) {n_iterations = 5}

if (!exists("include_species_jackknife")) {
  include_species_jackknife = FALSE
  print("Defaulting to only running data-availability bias downsampling, and not running phylopath jackknife by species. To also run phylopath jackknife by species, define include_species_jackknife = TRUE")
}

dfIn_phylo = read.csv("Data_R.csv")
tree_phylo = read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

all_traits_phylopath_label = paste(female_song_var, coop_breeding_var, territoriality_var, mass_var)

trait_set_output_dir = file.path("Outputs", "PhylopathDownsampled", paste0(all_traits_phylopath_label, " models"))

if (!dir.exists(trait_set_output_dir)) {
  dir.create(trait_set_output_dir, recursive = T)
}

# Add geographic regions if needed
if (!"GeographicRegion_Jetz" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$GeographicRegion_Jetz <- NA
  dfIn_phylo$GeographicRegion_Jetz[which(dfIn_phylo$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
  dfIn_phylo$GeographicRegion_Jetz[which(dfIn_phylo$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
}

# Add data availability indicators
if (!"HaveFSData" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$HaveFSData <- !is.na(dfIn_phylo$FemaleSong_Agg01)
}
if (!"HaveFSCBData" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$HaveFSCBdata <- !is.na(dfIn_phylo$HighConfidence_Coop) & !is.na(dfIn_phylo$FemaleSong_Agg01)
}

################################# 1. Geographic bias correction - HOLARCTIC NONCOOPERATIVE ---
cat("\n\nRunning geographic bias correction - Holarctic non-cooperative...\n")
cat("Removing", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "species\n")
result_geo_holarctic <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
  downsample_values = c("Holarctic", 0, TRUE),
  numToRemove = downsampling_results$downsampling$holarctic_noncoop$n_to_remove,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "HolarcticNoncoop")
)

# Create plots with proper file naming
geo_holarctic_info <- list(
  proportions = data.frame(
    territoriality = c("Holarctic Non-coop", "Tropical Non-coop"),
    species_with_data = c(209, 428),
    total_species = c(369, 977),
    proportion = c(0.566, 0.438)
  ),
  target_proportion = 0.438,
  n_to_remove = downsampling_results$downsampling$holarctic_noncoop$n_to_remove, #83,
  downsample_info = list(
    group_to_downsample = "Holarctic non-cooperative",
    territory_value = "0"
  )
)

prefix_geo <- paste(all_traits_phylopath_label, "Remove83HolarcticNoncoop")
saveRDS(result_geo_holarctic,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_geo, "_n", n_iterations, "_", Sys.Date(), ".rds")))

plots_geo <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_geo_holarctic,
  downsampling_info = geo_holarctic_info,
  full_dataset = newdata,
  tree = treefile,
  output_prefix = prefix_geo,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 2. Geographic bias correction - TROPICAL COOPERATIVE ---
cat("\n\nRunning geographic bias correction - Tropical cooperative...\n")
result_geo_tropical <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
  downsample_values = c("Tropical", 1, TRUE),
  numToRemove = downsampling_results$downsampling$tropical_coop$n_to_remove, #24,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove24TropicalCoop")
)

# Create plots
trop_coop_info <- list(
  proportions = data.frame(
    group = c("Tropical Coop", "Other"),
    species_with_data = c(68, 807),
    total_species = c(159, 2005),
    proportion = c(0.428, 0.402)
  ),
  target_proportion = 0.402,
  n_to_remove = downsampling_results$downsampling$tropical_coop$n_to_remove,
  downsample_info = list(
    group_to_downsample = "Tropical cooperative",
    territory_value = "1"
  )
)

prefix_trop <- paste(all_traits_phylopath_label, "Remove24TropicalCoop") 
saveRDS(result_geo_tropical,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_trop, "_n", n_iterations, "_", Sys.Date(), ".rds")))

plots_trop <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_geo_tropical,
  downsampling_info = trop_coop_info,
  full_dataset = newdata,
  tree = treefile,
  output_prefix = prefix_trop,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 3. GLOBAL COOPERATIVE bias correction ---
cat("\n\nRunning global cooperative bias correction...\n")
result_global_coop <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("HighConfidence_Coop", "HaveFSData"),
  downsample_values = c(1, TRUE),
  numToRemove = downsampling_results$downsampling$global_coop$n_to_remove, #15,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove15GlobalCoop")
)

# Create plots - info
global_coop_info <- list(
  proportions = data.frame(
    group = c("Cooperative", "Non-cooperative"),
    species_with_data = c(80, 795),
    total_species = c(226, 1938),
    proportion = c(0.354, 0.410)
  ),
  target_proportion = 0.069,
  n_to_remove = downsampling_results$downsampling$global_coop$n_to_remove, #15,
  downsample_info = list(
    group_to_downsample = "Cooperative",
    territory_value = "1"
  )
)

# actual plot creation
prefix_global <- paste(all_traits_phylopath_label, "Remove15GlobalCoop")
saveRDS(result_global_coop,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_global, "_n", n_iterations, "_", Sys.Date(), ".rds")))
plots_global <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_global_coop,
  downsampling_info = global_coop_info,
  full_dataset = newdata,
  tree = treefile,
  output_prefix = prefix_global,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 4. TERRITORIALITY WEAK/STRONG bias correction ---
cat("\n\nRunning territoriality bias correction...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
dfIn_phylo_char$TerritorialityWeakVsStrong <- as.character(dfIn_phylo_char$TerritorialityWeakVsStrong)
dfIn_phylo_char$Territory_12vs3 <- as.character(dfIn_phylo_char$Territory_12vs3)

if (downsampling_results$downsampling$territoriality_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territoriality_bias$n_to_remove, "species from StrongTerr\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("TerritorialityWeakVsStrong", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territoriality_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territoriality_bias$n_to_remove, "StrongTerr")
  )
  
  # Create plots
  prefix_terr <- paste0(all_traits_phylopath_label, "Remove",
                        downsampling_results$downsampling$territoriality_bias$n_to_remove, "StrongTerr")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = downsampling_results$downsampling$territoriality_bias$n_to_remove,
    full_dataset = newdata,
    tree = treefile,
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

################################# 5. TERRITORIALITY 12 VS 3 bias correction using Territory_12vs3 as TERR variable ---
cat("\n\nRunning Territory_12vs3 bias correction...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
dfIn_phylo_char$Territory_12vs3 <- as.character(dfIn_phylo_char$Territory_12vs3)

if (downsampling_results$downsampling$territory_12vs3_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "species from", "Terr3\n")
  # terr_downsample_info$report$downsample_info$group_to_downsample, "\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("Territory_12vs3", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territory_12vs3_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  )
  
  # Create plots
  prefix_terr <- paste0("FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = NULL, #terr_downsample_info$report,
    full_dataset = newdata,
    tree = treefile,
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}


################################# 6. TERRITORIALITY 12 VS 3 bias correction using TerritorialityWeakVsStrong as TERR variable ---
cat("\n\nRunning Territory_12vs3 bias correction - uses TerritorialityWeakVsStrong as TERR var...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
#dfIn_phylo_char$Terr <- as.character(dfIn_phylo_char$Territory_12vs3)
dfIn_phylo_char$TerritorialityWeakVsStrong <- as.character(dfIn_phylo_char$TerritorialityWeakVsStrong)

if (downsampling_results$downsampling$territory_12vs3_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "species from Terr3", "\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("Territory_12vs3", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territory_12vs3_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3-uses-TerritorialityWeakVsStrong-as-TERR")
  )
  
  # Create plots
  prefix_terr <- paste0("FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET Remove",
                        downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = NULL, #terr_downsample_info$report,
    full_dataset = newdata,
    tree = treefile,
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

# 5. Process all detailed model files with enhanced plots
cat("\n\nProcessing detailed model files for enhanced plots...\n")
detailed_files <- list.files(path = trait_set_output_dir, pattern = paste0("detailed_models_.*", ".*\\.csv$"),  full.names = TRUE)


#### Run phylopath with downsampling for dimorphism ----

# Define dimorphism variables to test
dimorphism_vars <- list(
  plumage = list(
    col = "logMaleFemalePlumageDiffAbs",
    label = "PlumageDimorphism",
    description = "log Plumage Dimorphism (Absolute Value)"
  ),
  wing = list(
    col = "PercentAbsLogWingDimorphism",
    label = "WingDimorphism", 
    description = "Percent Absolute-Value Log Wing Dimorphism"
  )
)

# Loop through each dimorphism variable
for (dim_type in names(dimorphism_vars)) {
  dim_info <- dimorphism_vars[[dim_type]]
  
  # Run the improved phylopath dimorphism correction
  result_dimorphism <- run_phylopath_dimorphism_correction(
    dfIn_phylo = dfIn_phylo,
    tree = tree_phylo,
    dim_info = dim_info,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    phylopath_output_dir = phylopath_output_dir,
    save_outputs = TRUE
  )
  
  # If results were obtained, create plots
  if (!is.null(result_dimorphism)) {
    # Create output prefix with dimorphism type
    prefix_dimorphism <- paste0("Remove", 
                                result_dimorphism$downsampling_info$n_to_remove, 
                                "High", 
                                dim_info$label,"ByPropensity")
    saveRDS(result_dimorphism,
            file.path(trait_set_output_dir,
                      paste0("result_", prefix_dimorphism, "_n", n_iterations, "_", Sys.Date(), ".rds")))
    
    # Use the flexible plotting function
    plots_dimorphism <- create_downsampled_plots(
      downsampling_results = result_dimorphism,
      detailed_models_input = result_dimorphism$detailed_models,  # Pass dataframe directly
      downsampling_info = NULL,
      full_dataset = dfIn_phylo,
      tree = tree_phylo,
      full_data_phylopath_input = NULL,
      output_prefix = prefix_dimorphism,
      output_dir = file.path(phylopath_output_dir, dim_info$label),
      save_png = TRUE,
      save_pdf = TRUE
    )
  }
}

cat("\n\nAll phylopath bias downsampling complete, beginning jackknifing by species \n")


#### Jackknife by species ----

if (include_species_jackknife) {
  source(file.path("Phylopath_functions", "run_jackknife_species_phylopath.R"))
  
  jackknife_results <- run_jackknife_species_phylopath(
    dfIn = dfIn_phylo,
    tree = tree,
    female_song_var = "FemaleSong_Agg01",
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = territoriality_var,
    mass_var = "logMass_AVONET",
    save_conditional_plots = FALSE,  # Set to TRUE if you want plots (will create large PDF)
    save_path_coefficients = TRUE,
    output_dir = file.path("Outputs", "PhylopathJackknife")
  )
}

#### Generate summary figures - forest plots, heat map ----
cat("\n\ngenerating plots using RunAnalyses_generate_bias_figures_both_traits.R \n")
n_iterations_to_use = n_iterations
source(file.path("Phylopath_functions", "RunAnalyses_generate_bias_figures_both_traits.R"))


