# adding more factors to phylopath models to test

boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

tree = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

complete_vars <- c("Territory", "FemaleSong_Agg01", "HighConfidence_Coop", "Social.bond", "Chorus", "Duet")
binary_vars <- c("FemaleSong_Agg01", "HighConfidence_Coop", "Chorus", "Duet")
dfIn_clean <- dfIn[complete.cases(dfIn[,complete_vars]),]
# dfIn_clean[["FemaleSong_Agg01"]] <- as.factor(dfIn_clean[["FemaleSong_Agg01"]])
# dfIn_clean[["HighConfidence_Coop"]] <- as.factor(dfIn_clean[["HighConfidence_Coop"]])
# dfIn_clean[["Chorus"]] <- as.factor(dfIn_clean[["Chorus"]])
# dfIn_clean[["Duet"]] <- as.factor(dfIn_clean[["Duet"]])
dfIn_clean[["Territory"]] <- as.numeric(dfIn_clean[["Territory"]])
dfIn_clean[["Social.bond"]] <- as.numeric(dfIn_clean[["Social.bond"]])

for(var in binary_vars) {
  print(dfIn_clean[[var]])
  # Convert to 0/1 numeric first
  data_clean[[var]] <- as.numeric(as.character(data_clean[[var]]))
  
  # Then explicitly convert to binary factor with labels
  dfIn_clean[[var]] <- factor(dfIn_clean[[var]], levels = c(0, 1), labels = c("absent", "present"))
  
  # Print confirmation
  cat("Converted", var, "to binary factor with levels:", paste(levels(dfIn_clean[[var]]), collapse=", "), "\n")
}

tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfIn_clean$species))

song_models <- define_model_set(
  direct = c(FemaleSong_Agg01 ~ Territory + Social.bond + Duet),
  indirect = c(FemaleSong_Agg01 ~ Territory + Social.bond, Duet ~ Territory + Social.bond),
 # reciprocal = c(FemaleSong_Agg01 ~ Duet, Duet ~ FemaleSong_Agg01),
  .common = c(Territory ~ Social.bond)
)

models <- define_model_set(
    # Basic models with single predictors
    m01 = c(FemaleSong_Agg01 ~ Territory),
    m02 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
    m03 = c(FemaleSong_Agg01 ~ Social.bond),
    m04 = c(FemaleSong_Agg01 ~ Chorus),
    m05 = c(FemaleSong_Agg01 ~ Duet),
    
    # Models with pair-wise combinations
    m06 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop),
    m07 = c(FemaleSong_Agg01 ~ Territory + Social.bond),
    m08 = c(FemaleSong_Agg01 ~ Territory + Chorus),
    m09 = c(FemaleSong_Agg01 ~ Territory + Duet),
    m10 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond),
    m11 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Chorus),
    m12 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Duet),
    m13 = c(FemaleSong_Agg01 ~ Social.bond + Chorus),
    m14 = c(FemaleSong_Agg01 ~ Social.bond + Duet),
    m15 = c(FemaleSong_Agg01 ~ Chorus + Duet),
    
    # Models with interactions between pairs
    m16 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop),
    m17 = c(FemaleSong_Agg01 ~ Territory * Social.bond),
    m18 = c(FemaleSong_Agg01 ~ Territory * Chorus),
    m19 = c(FemaleSong_Agg01 ~ Territory * Duet),
    m20 = c(FemaleSong_Agg01 ~ HighConfidence_Coop * Social.bond),
    m21 = c(FemaleSong_Agg01 ~ HighConfidence_Coop * Chorus),
    m22 = c(FemaleSong_Agg01 ~ HighConfidence_Coop * Duet),
    m23 = c(FemaleSong_Agg01 ~ Social.bond * Chorus),
    m24 = c(FemaleSong_Agg01 ~ Social.bond * Duet),
    m25 = c(FemaleSong_Agg01 ~ Chorus * Duet),
    
    # Models with three predictors
    m26 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond),
    m27 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Chorus),
    m28 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Duet),
    m29 = c(FemaleSong_Agg01 ~ Territory + Social.bond + Chorus),
    m30 = c(FemaleSong_Agg01 ~ Territory + Social.bond + Duet),
    m31 = c(FemaleSong_Agg01 ~ Territory + Chorus + Duet),
    m32 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond + Chorus),
    m33 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond + Duet),
    m34 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Chorus + Duet),
    m35 = c(FemaleSong_Agg01 ~ Social.bond + Chorus + Duet),
    
    # Selected interaction models with three predictors
    m36 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond),
    m37 = c(FemaleSong_Agg01 ~ Territory * Social.bond + HighConfidence_Coop),
    m38 = c(FemaleSong_Agg01 ~ HighConfidence_Coop * Social.bond + Territory),
    m39 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Chorus),
    m40 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Duet),
    
    # Models with four predictors
    m41 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond + Chorus),
    m42 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond + Duet),
    m43 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Chorus + Duet),
    m44 = c(FemaleSong_Agg01 ~ Territory + Social.bond + Chorus + Duet),
    m45 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond + Chorus + Duet),
    
    # # Selected models with key interactions and four predictors
    m46 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond + Chorus),
    m47 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond + Duet),

    # # Full model with all five predictors
    m48 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond + Chorus + Duet),

    # Full model with important interactions
    m49 = c(FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond + Chorus + Duet),
    m50 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop * Social.bond + Chorus + Duet),

    # Path analysis models considering HighConfidence_Coop as dependent on Territory
    m51 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory,
            HighConfidence_Coop ~ Territory),

    # Path analysis models considering Social.bond as dependent on Territory
    m52 = c(FemaleSong_Agg01 ~ Social.bond + Territory,
            Social.bond ~ Territory),

    # Path analysis with Duet and Chorus
    m53 = c(FemaleSong_Agg01 ~ Duet + Territory,
            Duet ~ Territory),

    m54 = c(FemaleSong_Agg01 ~ Chorus + Territory,
            Chorus ~ Territory),

    # More complex path models
    m55 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond,
            HighConfidence_Coop ~ Territory,
            Social.bond ~ Territory),

    m56 = c(FemaleSong_Agg01 ~ Duet + HighConfidence_Coop,
            Duet ~ Territory,
            HighConfidence_Coop ~ Territory),

    # Path model with more complex relationships
    m57 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Duet,
            HighConfidence_Coop ~ Territory + Social.bond,
            Social.bond ~ Territory),

    # Models with Territory and HighConfidence_Coop as correlated predictors
    m58 = c(FemaleSong_Agg01 ~ Territory + Duet,
            Duet ~ HighConfidence_Coop,
            HighConfidence_Coop ~ Territory),

    # Models testing evolutionary relationships
    m59 = c(FemaleSong_Agg01 ~ Territory,
            HighConfidence_Coop ~ FemaleSong_Agg01),

    m60 = c(FemaleSong_Agg01 ~ Territory,
            Duet ~ FemaleSong_Agg01,
            Chorus ~ FemaleSong_Agg01),
    .common = NULL
)
  

# Create a subset of models you want to test
selected_models <- models#[c("m01", "m02", "m16", "m51", "m59")]  # Example selection

# Use the models in your phylopath analysis
result <- phylo_path(models, data = dfIn_clean, tree = tree_clean, model = 'lambda')



# Define a more comprehensive set of models
test_models <- define_model_set(
  # Basic single-predictor models
  m01 = FemaleSong_Agg01 ~ Territory,
  m02 = FemaleSong_Agg01 ~ HighConfidence_Coop,
  m03 = FemaleSong_Agg01 ~ Social.bond,
  m04 = FemaleSong_Agg01 ~ Chorus,
  m05 = FemaleSong_Agg01 ~ Duet,
  
  # Two-predictor models
  m06 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop,
  m07 = FemaleSong_Agg01 ~ Territory + Social.bond,
  m08 = FemaleSong_Agg01 ~ Territory + Chorus, 
  m09 = FemaleSong_Agg01 ~ Territory + Duet,
  m10 = FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond,
  m11 = FemaleSong_Agg01 ~ HighConfidence_Coop + Chorus,
  m12 = FemaleSong_Agg01 ~ HighConfidence_Coop + Duet,
  
  # Key interaction models
  m13 = FemaleSong_Agg01 ~ Territory * HighConfidence_Coop,
  m14 = FemaleSong_Agg01 ~ Territory * Social.bond,
  m15 = FemaleSong_Agg01 ~ HighConfidence_Coop * Social.bond,
  
  # Three-predictor models
  m16 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond,
  m17 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Chorus,
  m18 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Duet,
  
  # Three-predictor models with interactions
  m19 = FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond,
  m20 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop * Social.bond,
  
  # Full models
  m21 = FemaleSong_Agg01 ~ Territory + HighConfidence_Coop + Social.bond + Chorus + Duet,
  m22 = FemaleSong_Agg01 ~ Territory * HighConfidence_Coop + Social.bond + Chorus + Duet,
  
  .common = NULL
)

# Define path models separately
path_models <- define_model_set(
  # Sequential relationships
  p01 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ Territory),
  
  p02 = c(FemaleSong_Agg01 ~ Social.bond,
         Social.bond ~ Territory),
  
  p03 = c(FemaleSong_Agg01 ~ Territory,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  
  p04 = c(FemaleSong_Agg01 ~ Territory,
         Social.bond ~ FemaleSong_Agg01),
  
  p05 = c(FemaleSong_Agg01 ~ Duet,
         Duet ~ Territory),
  
  p06 = c(FemaleSong_Agg01 ~ Chorus,
         Chorus ~ Territory),

  # Models with mediation
  p07 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory,
         HighConfidence_Coop ~ Territory),

  p08 = c(FemaleSong_Agg01 ~ Social.bond + Territory,
         Social.bond ~ Territory),
  
  # Models with multiple mediators
  p09 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond,
         HighConfidence_Coop ~ Territory,
         Social.bond ~ Territory),

  p10 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Duet,
         HighConfidence_Coop ~ Territory,
         Duet ~ Territory),
  # 
  # # Models with feedback loops
  # p11 = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop,
  #        HighConfidence_Coop ~ FemaleSong_Agg01),

  # p12 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Social.bond,
  #        HighConfidence_Coop ~ Territory,
  #        Social.bond ~ FemaleSong_Agg01),
         
  .common = NULL
)

# Run phylopath with the regression models first
result_reg <- phylo_path(test_models, data = dfIn_complete, tree = tree_complete, model = 'lambda')

# If successful, view and save results
if(!is.null(result_reg)) {
  s_reg <- summary(result_reg)
  print(s_reg)
  
  # Save regression model plots
  model_plot <- plot_model_set(test_models)
  ggsave("phylopath_regression_models.png", model_plot, width = 15, height = 12)
  
  summary_plot <- plot(s_reg)
  ggsave("phylopath_regression_summary.png", summary_plot, width = 15, height = 12)
  
  if(nrow(s_reg) > 0) {
    best_model <- best(result_reg)
    best_plot <- plot(best_model)
    ggsave("phylopath_regression_best_model.png", best_plot, width = 12, height = 10)
  }
}

# Run phylopath with the path models 
result_path <- phylo_path(path_models, data = dfIn_complete, tree = tree_complete, model = 'lambda')

# If successful, view and save results
if(!is.null(result_path)) {
  s_path <- summary(result_path)
  print(s_path)
  
  # Save path model plots
  path_model_plot <- plot_model_set(path_models)
  ggsave("phylopath_path_models.png", path_model_plot, width = 15, height = 12)
  
  path_summary_plot <- plot(s_path)
  ggsave("phylopath_path_summary.png", path_summary_plot, width = 15, height = 12)
  
  if(nrow(s_path) > 0) {
    path_best_model <- best(result_path)
    path_best_plot <- plot(path_best_model)
    ggsave("phylopath_path_best_model.png", path_best_plot, width = 12, height = 10)
  }
}

# Run phylopath directly
result <- phylo_path(test_models, data = dfIn_complete, tree = tree_complete, model = 'lambda')

# If successful, view the results
if(!is.null(result)) {
  s <- summary(result)
  print(s)
  
  # Save plots using ggsave
  model_plot <- plot_model_set(test_models)
  ggsave("phylopath_models.png", model_plot, width = 12, height = 8)
  
  summary_plot <- plot(s)
  ggsave("phylopath_summary.png", summary_plot, width = 12, height = 8)
  
  if(nrow(s) > 0) {
    best_model <- best(result)
    best_plot <- plot(best_model)
    ggsave("phylopath_best_model.png", best_plot, width = 10, height = 8)
  }
}

