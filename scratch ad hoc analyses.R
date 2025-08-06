# ad-hoc data analyses


browniedf <- read.csv('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/OutputFiles/2025-02-08 Territory Song.rep.final multistate aceARD Brownie test 500 sims.csv')

browniedf <- read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/2025-06-13 TerrWeakStrongXHighConfCoop Song.rep.final multistate aceARD Brownie IntersectionTrait 1500 sims.csv')
table(browniedf$convergence)
browniedf = browniedf[which(browniedf$convergence == "Optimization has converged."),]


#calculate overall mean pval
ERloglikmean <- mean(browniedf$ERloglik)
ARDloglikmean <- mean(browniedf$ARDloglik)
ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 3)



'/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Simmap Overlap Outputs/transition-counts_Coop_FS_TerritorialityWeakVsStrong_500sims.csv'

'/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Simmap Overlap Outputs/2025-05-09_transition-counts HighConfidence_Coop FemaleSong_Agg01 TerritorialityWeakVsStrong_500sims.csv'

data = read.csv("Data_R_2025-06-09.csv")
tree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
source("run_phyloglm_fxns_updated.R")
source("phyloglm_plot_fixes.R")

results <- run_phyloglm_model_set(
  response_var = "FemaleSong_Agg01",
  pred_var = "HighConfidence_Coop",  # Your custom predictor
  terr_var = "TerritorialityWeakVsStrong",
  original_data = data,
  original_tree = tree,
  boot = 1000
)

# 1. AIC Ranking Plot
create_aic_ranking_plot_fixed(
  model_results = results,  # Pass the full results object
  threshold = 2,
  save_plot = TRUE,
  output_dir = "Outputs/PhyloglmResults"
)

# 2. Interaction Plot for each model
# For the best model:
best_model_name <- results$comparison$Model[1]
best_model <- results$models[[best_model_name]]

create_phyloglm_interaction_plot_fixed(
  model = best_model,
  data = results$data,
  response_var = results$response,
  pred_var = results$predictor,
  terr_var = results$territoriality,
  save_plot = TRUE,
  output_dir = "Outputs/PhyloglmResults"
)

# 3. Coefficient Plot (after extracting coefficients)
coef_summary <- extract_phyloglm_coefficients(
  models = list(best_model_name = best_model),
  response_var = results$response,
  pred_var = results$predictor,
  terr_var = results$territoriality
)

create_coefficient_plot(
  coef_summary = coef_summary,
  model_name = best_model_name,
  save_plot = TRUE,
  output_dir = "Outputs/PhyloglmResults"
)

# 4. For all top models (deltaAIC < 2)
top_models <- results$comparison$Model[results$comparison$deltaAIC < 2]
for (model_name in top_models) {
  model <- results$models[[model_name]]
  
  create_phyloglm_interaction_plot_fixed(
    model = model,
    data = results$data,
    response_var = results$response,
    pred_var = results$predictor,
    terr_var = results$territoriality,
    save_plot = TRUE,
    output_dir = "Outputs/PhyloglmResults"
  )
}


# Hacky FS/TerrWeakStrong simmap overlap

simmapOverlapdf= read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/2025-06-13_transition-counts HighConfidence_Coop FemaleSong_Agg01 TerritorialityWeakVsStrong_1000sims.csv')

# compare AVONET latitude to Odom et al 2025 latitude
newdata = "Data_R_2025-07-21.csv"
data = read.csv(newdata)

passertree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex")

OdomData = read.csv("Odom et al 2025 final published data.csv")
sum(OdomData$TipLabel.for.tree %in% data$species)
sum(!OdomData$TipLabel.for.tree %in% data$species)

sum(data$species %in% OdomData$TipLabel.for.tree)
sum(!OdomData$TipLabel.for.tree %in% data$species)
OdomData$TipLabel.for.tree[which(!OdomData$TipLabel.for.tree %in% data$species)]
OdomData$TipLabel.for.tree[which(!OdomData$TipLabel.for.tree %in% data$species)] %in% passertree$tip.label
OdomData$TipLabel.for.tree[which(!OdomData$TipLabel.for.tree %in% passertree$tip.label)]

inconsistent_names <- read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/inconsistent_species_names_InclTobias.csv')
OdomData$TipLabel.for.tree[which(!OdomData$TipLabel.for.tree %in% data$species)] %in% inconsistent_names$in_database

bothdata <- merge(data, OdomData, by.x = "species", by.y = "TipLabel.for.tree", all=T)

plot(bothdata$Centroid.Latitude_AVONET, bothdata$degrees_from_equator)
plot(abs(bothdata$Centroid.Latitude_AVONET), bothdata$degrees_from_equator)
abline(b=1)

bothdata %>% group_by(FemaleSong_Agg01, FemSongFinal_PrsAbs) %>% count


