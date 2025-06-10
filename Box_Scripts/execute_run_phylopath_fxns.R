## Code to execute repeated downsampling phylopath 

source("run_phylopath_fxns.R")

boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')
allcols = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/2024-8-4AllCol_Duet_Female_AllTraits_fromJiaying2025-03-04.csv")
tree = read.nexus(treefile)

df = read.csv(newdata)
dfIn = merge(df, allcols[,c("BirdtreeSpecies", "Female_plumage_score_Dale2015", "Male_plumage_score_Dale2015", "Realm_Jetz2011", "FemalePC1sum_Dunn2015", "FemalePC2sum_Dunn2015", "MalePC1sum_Dunn2015", "MalePC2sum_Dunn2015", "sumDiffPC1_Dunn2015", "sumDiffPC2_Dunn2015", "Tropical_life_history_ppca_Dale2015", "Region_Cockburn2006", "location_Downing2015")], by.x = "species", by.y = "BirdtreeSpecies")
rownames(dfIn) = dfIn$species
dfIn$logMass_AVONET = log(dfIn$Mass_AVONET)
dfIn$absCentroid.Latitude.AVONET = abs(dfIn$Centroid.Latitude_AVONET)
dfIn$HaveFSCBdata = !is.na(dfIn$HighConfidence_Coop) & !is.na(dfIn$FemaleSong_Agg01)


## Run one downsample
downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData")
downsample_values = c("Holarctic", 0, TRUE)
numToRemove = 83
seed = 123

output <- downsample_run_phylopath(dfIn = dfIn, tree = tree, downsample_columns = downsample_columns, downsample_values = downsample_values, numToRemove = numToRemove, seed = seed, female_song_var = "FemaleSong_Agg01", coop_breeding_var = "HighConfidence_Coop", territoriality_var = "TerritorialityWeakVsStrong", mass_var = "logMass_AVONET")


#### Run many downsamples - logMass continuous var  - downsampling holarctic non-cooperative species ----

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive", "TerritorialityPermissiveExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
    downsample_values = c("Holarctic", 0, TRUE),
    numToRemove = 83,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "GeographicRegion_Jetz ")
  )
}

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive", "TerritorialityPermissiveExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
    downsample_values = c("Holarctic", 0, TRUE),
    numToRemove = 88,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "GeographicRegion_Jetz ")
  )
}

# To examine the results:
head(results$results_df)
head(results$model_frequencies)

# If you saved path coefficients:
head(results$path_summary)

# For a quick summary:
summary_stats <- results$results_df %>%
  summarize(
    n_iterations = n(),
    mean_species = mean(nSpecies),
    mean_sub2_models = mean(nSub2dCIC_Models),
    most_common_best_model = names(sort(table(best_model), decreasing = TRUE)[1]),
    proportion_most_common = max(table(best_model))/n()
  )

print(summary_stats)



#### Run many downsamples - Centroid.Latitude continuous var ----

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive", "TerritorialityPermissiveExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
    downsample_values = c("Holarctic", 0, TRUE),
    numToRemove = 83,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "",
    plotlabel = paste(terrcol, "GeographicRegion_Jetz ")
  )
}


#### Run many downsamples - logMass continuous var  - downsampling TROPICAL COOPERATIVE species - Jetz regions ----

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
    downsample_values = c("Tropical", 1, TRUE),
    numToRemove = 24,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "GeographicRegion_Jetz Remove24TropicalCoop ")
  )
}

#### Run many downsamples - logMass continuous var  - downsampling tropical cooperative species - Cockburn regions ----
# Because tropical species with cooperative breeding may be preferentially studied and therefore more likely to have FS data

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("GeographicRegion_Cockburn", "HighConfidence_Coop", "HaveFSData"),
    downsample_values = c("Tropical", 1, TRUE),
    numToRemove = 26,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "GeographicRegion_Cockburn ")
  )
}

#### Run many downsamples - logMass continuous var  - downsampling GLOBAL COOPERATIVE species ----

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive") #, "TerritorialityPermissiveExclusive")

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("HighConfidence_Coop", "HaveFSData"),
    downsample_values = c(1, TRUE),
    numToRemove = 15,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "Remove15GlobalCoop ")
  )
}


#### Run many downsamples - logMass continuous var  - downsampling GLOBAL TERRIRORY3 species ----

TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "Territory_12vs3", "TerritorialityPermissiveColonialCoopVsExclusive") #, "TerritorialityPermissiveExclusive")
TerritorialityCols = "Territory_12vs3"

for (terrcol in TerritorialityCols) {
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c("Territory_12vs3", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = 155,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01", 
    coop_breeding_var = "HighConfidence_Coop", 
    territoriality_var = terrcol, 
    mass_var = "logMass_AVONET",
    plotlabel = paste(terrcol, "Remove155Territory3 ")
  )
}


#### Analyze outputs of downsamples ----
# Example usage:
detailed_models_df <- read.csv("detailed_models_TerritorialityWeakVsStrong GeographicRegion_Jetz Remove24TropicalCoop _500_2025-05-16.csv")
results <- create_all_plots(detailed_models_df, save_plots = TRUE, file_prefix = "TerritorialityWeakVsStrong GeographicRegion_Jetz Remove24TropicalCoop _500_2025-05-16")

# View individual plots:
print(results$plots$coef_violin)
print(results$plots$sig_heatmap)

# View summaries:
head(results$path_summary)
head(results$model_summary)


detailed_model_list = list.files(pattern = "detailed_models_")
for (i in 3: length(detailed_model_list)) {
  tempfile = detailed_model_list[i]
  basename = gsub(".csv", "", tempfile)
  basename = gsub("detailed_models_", "", basename)
  output_basename = paste("phylopath_analyses", basename)
  detailed_models_df = read.csv(tempfile)
  results <- create_all_plots(detailed_models_df, save_plots = TRUE, file_prefix = output_basename)
}
  
  
  