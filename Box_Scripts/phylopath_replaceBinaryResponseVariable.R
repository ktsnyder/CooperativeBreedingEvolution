## phylopath with additional factors from Tobias

library(phylopath)
library(phytools)


boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

#newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

consensus_phylo = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

# 3-level Territory
dfIn$Territory <- as.character(dfIn$Territory)
dfIn$Territory <- factor(dfIn$Territory, levels = c("1", "2", "3"))
unique(dfIn$HighConfidence_Coop)
unique(dfIn$FemaleSong_Agg01)
unique(dfIn$Territory)

dfIn = dfIn[which(!is.na(dfIn$Territory)),]
# Create split binary variables for Territory

Territory.matrix <- model.matrix(~ Territory - 1, data = dfIn)
dfIn$Territory1 <- Territory.matrix[,1]
dfIn$Territory2 <- Territory.matrix[,2]
dfIn$Territory3 <- Territory.matrix[,3]



# 2-level territory with logMass_AVONET
dfIn_clean <- dfIn[complete.cases(dfIn[, c("HighConfidence_Coop", "TerritorialityWeakVsStrong", "FemaleSong_Agg01", "Mass_AVONET")]), ]
dfIn_clean <- dfIn[complete.cases(dfIn[, c("HighConfidence_Coop", "Territory", "FemaleSong_Agg01", "Mass_AVONET")]), ]
#dfIn_clean <- dfIn[complete.cases(dfIn[, c("HighConfidence_Coop", "TerritorialityWeakVsStrong", "FemaleSong_Agg01")]), ]
dfIn_clean$logMass_AVONET = log(dfIn_clean$Mass_AVONET)

# Prune tree to match dfIn_clean species
tree_clean <- drop.tip(consensus_phylo, setdiff(consensus_phylo$tip.label, dfIn_clean$species))


# incomplete model set
run_phylopath_analysis <- function(data, tree, response_var, predictor_vars, 
                                   territory_vars = c("Territory2", "Territory3")) {
  
  library(phylopath)
  library(ggplot2)
  
  # Make a copy of the data to avoid modifying the original
  data_clean <- data
  
  # Ensure data is properly prepared
  complete_vars <- c(territory_vars, response_var, predictor_vars)
  data_clean <- data_clean[complete.cases(data_clean[,complete_vars]),]
  
  # Force binary variables to be recognized as binary factors
  binary_vars <- c(response_var, predictor_vars, territory_vars)
  for(var in binary_vars) {
    # Convert to 0/1 numeric first
    data_clean[[var]] <- as.numeric(as.character(data_clean[[var]]))
    
    # Then explicitly convert to binary factor with labels
    data_clean[[var]] <- factor(data_clean[[var]], levels = c(0, 1), labels = c("absent", "present"))
    
    # Print confirmation
    cat("Converted", var, "to binary factor with levels:", paste(levels(data_clean[[var]]), collapse=", "), "\n")
  }
  
  # Subset tree
  tree_clean <- drop.tip(tree, setdiff(tree$tip.label, data_clean$species))
  rownames(data_clean) <- data_clean$species
  
  # Create formula strings
  f0 <- as.formula(paste(response_var, "~", paste(territory_vars, collapse = " + ")))
  
  f1_1 <- as.formula(paste(response_var, "~", paste(c(territory_vars, predictor_vars), collapse = " + ")))
  f1_2 <- as.formula(paste(predictor_vars, "~", paste(territory_vars, collapse = " + ")))
  
  f2_1 <- as.formula(paste(response_var, "~", predictor_vars))
  f2_2 <- as.formula(paste(predictor_vars, "~", paste(territory_vars, collapse = " + ")))
  
  f3_1 <- as.formula(paste(response_var, "~", paste(territory_vars, collapse = " + ")))
  f3_2 <- as.formula(paste(predictor_vars, "~", paste(territory_vars, collapse = " + ")))
  
  f4_1 <- as.formula(paste(response_var, "~", paste(territory_vars, collapse = " + ")))
  f4_2 <- as.formula(paste(predictor_vars, "~", response_var))
  
  f5_1 <- as.formula(paste(response_var, "~", paste(territory_vars, collapse = " + ")))
  f5_2 <- as.formula(paste(predictor_vars, "~", paste(c(response_var, territory_vars), collapse = " + ")))
  
  # Now create the models using eval() to evaluate the formulas in the current environment
  models <- eval(substitute(
    define_model_set(
      m0 = c(F0),
      m1 = c(F1_1, F1_2),
      m2 = c(F2_1, F2_2),
      m3 = c(F3_1, F3_2),
      m4 = c(F4_1, F4_2),
      m5 = c(F5_1, F5_2)
    ),
    list(
      F0 = f0,
      F1_1 = f1_1, F1_2 = f1_2,
      F2_1 = f2_1, F2_2 = f2_2,
      F3_1 = f3_1, F3_2 = f3_2,
      F4_1 = f4_1, F4_2 = f4_2,
      F5_1 = f5_1, F5_2 = f5_2
    )
  ))
  
  # Run analysis with error handling
  cat("\nRunning phylo_path...\n")
  result <- tryCatch({
    phylo_path(models, data = data_clean, tree = tree_clean, model = 'lambda', na.rm = TRUE)
  }, 
  warning = function(w) {
    cat("Warning from phylo_path:", conditionMessage(w), "\n")
    # Continue despite warnings
    return(phylo_path(models, data = data_clean, tree = tree_clean, model = 'lambda', na.rm = TRUE))
  },
  error = function(e) {
    cat("Error from phylo_path:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # If successful, save plots and return results
  if (!is.null(result)) {
    # Try to generate and save plots with error handling
    tryCatch({
      cat("\nGenerating plots...\n")
      
      # Create file names
      plot_file <- paste0("phylopath_", response_var, " " , paste(territory_vars, collapse = " "), "_plots_", format(Sys.Date(), "%Y%m%d"), ".png")
      summary_file <- paste0("phylopath_", response_var, " " , paste(territory_vars, collapse = " "), "_summary_", format(Sys.Date(), "%Y%m%d"), ".png")
      best_model_file <- paste0("phylopath_", response_var, " " , paste(territory_vars, collapse = " "), "_best_model_", format(Sys.Date(), "%Y%m%d"), ".png")
      
      # Plot and save model set as ggplot
      p1 <- plot_model_set(models) +  coord_cartesian(clip = "off") +
        theme(
          plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
        )
      ggsave(plot_file, p1, width = 10, height = 10, dpi = 300)
      cat("Model set plot saved to:", plot_file, "\n")
      
      # Generate and save summary plot
      s <- summary(result)
      p2 <- plot(s)+  coord_cartesian(clip = "off") +
        theme(
          plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
        )
      ggsave(summary_file, p2, width = 10, height = 10, dpi = 300)
      cat("Summary plot saved to:", summary_file, "\n")
      
      # Generate and save best model plot
      b <- best(result)
      p3 <- plot(b)+  coord_cartesian(clip = "off") +
        theme(
          plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
        )
      ggsave(best_model_file, p3, width = 10, height = 10, dpi = 300)
      cat("Best model plot saved to:", best_model_file, "\n")
    }, 
    warning = function(w) {
      cat("Warning when generating plots:", conditionMessage(w), "\n")
    },
    error = function(e) {
      cat("Error when generating plots:", conditionMessage(e), "\n")
    })
    
    return(list(
      result = result,
      summary = summary(result),
      best_model = best(result),
      data = data_clean,
      tree = tree_clean,
      models = models
    ))
  } else {
    # Return just the data and models if analysis failed
    return(list(
      data = data_clean,
      tree = tree_clean,
      models = models
    ))
  }
}

# Example usage:
# To run with female song as response:
fs_analysis <- run_phylopath_analysis(
  data = dfIn_clean,
  tree = tree,
  response_var = "FemaleSong_Agg01",
  predictor_vars = "HighConfidence_Coop",
  territory_vars = c("Territory1", "Territory3")
)

# To run with duetting as response:
duet_analysis <- run_phylopath_analysis(
  data = dfIn_clean,
  tree = tree,
  response_var = "Duet",
  predictor_vars = "HighConfidence_Coop"
)




#### Define your variables as strings - this allows easy switching of variables ----
female_song_var <- "FemaleSong_Agg01"
coop_breeding_var <- "HighConfidence_Coop"  
territoriality_var <- "Territory"
mass_var <- "logMass_AVONET"

labels <- setNames(
  c("Female Song", "Cooperative Breeding", territoriality_var, "log(Body Mass)"),
  c(female_song_var, coop_breeding_var, territoriality_var, mass_var)
)

dfIn_clean <- dfIn[complete.cases(dfIn[, c(coop_breeding_var, territoriality_var, female_song_var, "Mass_AVONET")]), ]
dfIn_clean$logMass_AVONET = log(dfIn_clean$Mass_AVONET)
tree_clean <- drop.tip(consensus_phylo, setdiff(consensus_phylo$tip.label, dfIn_clean$species))

dfIn_clean[,female_song_var] <- as.factor(dfIn_clean[,female_song_var])
dfIn_clean[,coop_breeding_var] <- as.factor(dfIn_clean[,coop_breeding_var])
dfIn_clean[,territoriality_var] <- as.factor(dfIn_clean[,territoriality_var])

# Helper function to build formulas from patterns
build_formula <- function(pattern, var_map) {
  for (var_name in names(var_map)) {
    pattern <- gsub(var_name, var_map[[var_name]], pattern, fixed = TRUE)
  }
  return(as.formula(pattern))
}

# Helper function to create a set of formulas
create_formula_set <- function(patterns, var_map) {
  formulas <- lapply(patterns, function(pattern) {
    build_formula(pattern, var_map)
  })
  return(do.call(c, formulas))
}

# Create a variable mapping dictionary
var_map <- list(
  "FS" = female_song_var,
  "COOP" = coop_breeding_var,
  "TERR" = territoriality_var,
  "MASS" = mass_var
)

##### Renamed models, more complete; separated mass-containing models into model_patterns_with_MASS and models that have FS or COOP influencing TERR into model_patterns_with_TERR_response  ----
model_patterns <- list(
  ## SERIES A: Basic models without mass, no direct FS-COOP interaction
  "A0_TERR→FS" = c("FS ~ TERR"),
  "A0_TERR→COOP" = c("COOP ~ TERR"),
  "A0_TERR→FS_TERR→COOP" = c("FS ~ TERR", "COOP ~ TERR"),
  
  ## SERIES N: Null models
  "N0_INDEP" = c("FS ~ 1", "COOP ~ 1", "TERR ~ 1"),
  
  ## SERIES X: COOP influences FS (direction 1)
  "X1_COOP→FS" = c("FS ~ COOP"),
  "X1_COOP→FS_TERR→FS" = c("FS ~ TERR + COOP"),
  "X1_COOP→FS_TERR→COOP" = c("FS ~ COOP", "COOP ~ TERR"),

  ## SERIES Y: FS influences COOP (direction 2)
  "Y2_FS→COOP" = c("COOP ~ FS"),
  "Y2_FS→COOP_TERR→COOP" = c("COOP ~ TERR + FS"),
  "Y2_TERR→FS_FS→COOP" = c("FS ~ TERR", "COOP ~ FS")
)

model_patterns_with_MASS <- list(
  ## SERIES B: Direct mass effects, no direct FS-COOP interaction
  "B0_TERR→FS_MASS→FS" = c("FS ~ TERR + MASS"),
  "B0_TERR→COOP_MASS→COOP" = c("COOP ~ TERR + MASS"),
  "B0_TERR→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ TERR + MASS"),
  "B0_MASS→FS" = c("FS ~ MASS"),
  "B0_MASS→COOP" = c("COOP ~ MASS"),
  "B0_MASS→FS_MASS→COOP" = c("FS ~ MASS", "COOP ~ MASS"),
  
  ## SERIES C: Mass affects only COOP, no direct FS-COOP interaction
  "C0_TERR→FS_TERR→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ TERR + MASS"),
  "C0_TERR→FS_MASS→COOP" = c("FS ~ TERR", "COOP ~ MASS"),
  
  ## SERIES D: Mass affects only FS, no direct FS-COOP interaction
  "D0_TERR→FS_TERR→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ TERR"),
  "D0_MASS→FS_TERR→COOP" = c("FS ~ MASS", "COOP ~ TERR"),
  
  ## SERIES E: Mass affects TERR which affects others, no direct FS-COOP interaction
  "E0_TERR→FS_TERR→COOP_MASS→TERR" = c("FS ~ TERR", "COOP ~ TERR", "TERR ~ MASS"),
  "E0_TERR→FS_TERR→COOP_MASS→TERR_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ TERR", "TERR ~ MASS"),
  "E0_TERR→FS_MASS→TERR" = c("FS ~ TERR", "TERR ~ MASS"),
  "E0_TERR→COOP_MASS→TERR_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ MASS"),
  "E0_TERR→COOP_MASS→TERR" = c("COOP ~ TERR", "TERR ~ MASS"),
  
  ## SERIES N: Null models
  "N0_MASS→FS_MASS→COOP_MASS→TERR" = c("FS ~ MASS", "COOP ~ MASS", "TERR ~ MASS"),
  
  ## SERIES X: COOP influences FS (direction 1)
  "X1_COOP→FS_MASS→FS" = c("FS ~ COOP + MASS"),
  "X1_COOP→FS_TERR→FS_MASS→FS" = c("FS ~ COOP + TERR + MASS"),
  "X1_COOP→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ COOP + MASS", "COOP ~ TERR + MASS"),
  "X1_COOP→FS_TERR→COOP_MASS→COOP" = c("FS ~ COOP", "COOP ~ TERR + MASS"),
  "X1_COOP→FS_TERR→COOP_MASS→FS" = c("FS ~ COOP + MASS", "COOP ~ TERR"),
  "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR"),
  "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→COOP" = c("FS ~ TERR + COOP", "COOP ~ TERR + MASS"),
  "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR + MASS"),
  "X1_COOP→FS_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ COOP + MASS", "COOP ~ TERR", "TERR ~ MASS"),
  "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR", "TERR ~ MASS"),
  "X1_COOP→FS_TERR→FS_TERR→COOP" = c("FS ~ TERR + COOP", "COOP ~ TERR"), # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
  
  ## SERIES Y: FS influences COOP (direction 2)
  "Y2_FS→COOP_MASS→COOP" = c("COOP ~ FS + MASS"),
  "Y2_FS→COOP_TERR→COOP_MASS→COOP" = c("COOP ~ FS + TERR + MASS"),
  "Y2_TERR→FS_FS→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ FS + MASS"),
  "Y2_TERR→FS_FS→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ FS + MASS"),
  "Y2_TERR→FS_FS→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ FS"),
  "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR"),
  "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ FS + TERR + MASS"),
  "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR + MASS"),
  "Y2_TERR→FS_FS→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + MASS", "COOP ~ FS", "TERR ~ MASS"),
  "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR", "TERR ~ MASS"),
  "Y2_TERR→FS_FS→COOP_TERR→COOP" = c("FS ~ TERR", "COOP ~ FS + TERR") # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
)

model_patterns_with_TERR_response <- list(
  ## SERIES F: Basic models with TERR as response or an intermediate factor (no direct FS-COOP interaction)
  "F0_COOP→TERR" = c("TERR ~ COOP"),
  "F0_FS→TERR" = c("TERR ~ FS"),
  "F0_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP"),
  "F0_FS→TERR_TERR→COOP" = c("COOP ~ TERR", "TERR ~ FS"),
  "F0_COOP→TERR_TERR→FS" = c("FS ~ TERR", "TERR ~ COOP"),
  # Series F containing MASS:
  "F0_FS→TERR_TERR→COOP_MASS→FS" = c("COOP ~ TERR", "TERR ~ FS + MASS"),
  "F0_COOP→TERR_TERR→FS_MASS→COOP" = c("FS ~ TERR", "TERR ~ COOP + MASS"),
  "F0_FS→TERR_TERR→COOP_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ FS"),
  "F0_COOP→TERR_TERR→FS_MASS→FS" = c("FS ~ TERR + MASS", "TERR ~ COOP"),
  "F0_FS→TERR_TERR→COOP_MASS→FS_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ FS + MASS"),
  "F0_COOP→TERR_TERR→FS_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "TERR ~ COOP + MASS"),
  "F0_MASS→FS_FS→TERR_TERR→COOP" = c("COOP ~ TERR", "TERR ~ FS", "FS ~ MASS"),
  "F0_MASS→COOP_COOP→TERR_TERR→FS" = c("FS ~ TERR", "TERR ~ COOP", "COOP ~ MASS"),
  "F0_FS→TERR_MASS→TERR" = c("TERR ~ FS + MASS"),
  "F0_COOP→TERR_MASS→TERR" = c("TERR ~ COOP + MASS"),
  "F0_FS→TERR_COOP→TERR_MASS→TERR" = c("TERR ~ FS + COOP + MASS"),
  
  ## SERIES G: FS influences COOP and one/both influence TERR (FS → COOP)
  "G2_FS→COOP_COOP→TERR" = c("TERR ~ COOP", "COOP ~ FS"),
  "G2_FS→COOP_FS→TERR" = c("TERR ~ FS", "COOP ~ FS"),
  ## SERIES G containing MASS
  "G2_FS→COOP_COOP→TERR_MASS→FS" = c("TERR ~ COOP", "COOP ~ FS + MASS"),
  "G2_FS→COOP_COOP→TERR_MASS→TERR" = c("TERR ~ COOP + MASS", "COOP ~ FS"),
  "G2_FS→COOP_COOP→TERR_MASS→FS_MASS→TERR" = c("TERR ~ COOP + MASS", "COOP ~ FS + MASS"),
  "G2_MASS→FS_FS→COOP_COOP→TERR" = c("TERR ~ COOP", "COOP ~ FS", "FS ~ MASS"),
  "G2_FS→COOP_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP", "COOP ~ FS"), # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
  
  ## SERIES H: COOP influences FS and one/both influence TERR (COOP → FS)
  "H1_COOP→FS_FS→TERR" = c("TERR ~ FS", "FS ~ COOP"),
  "H1_COOP→FS_COOP→TERR" = c("TERR ~ COOP", "FS ~ COOP"),
  ## SERIES H containing MASS
  "H1_COOP→FS_FS→TERR_MASS→COOP" = c("TERR ~ FS", "FS ~ COOP + MASS"),
  "H1_COOP→FS_FS→TERR_MASS→TERR" = c("TERR ~ FS + MASS", "FS ~ COOP"),
  "H1_COOP→FS_FS→TERR_MASS→COOP_MASS→TERR" = c("TERR ~ FS + MASS", "FS ~ COOP + MASS"),
  "H1_MASS→COOP_COOP→FS_FS→TERR" = c("TERR ~ FS", "FS ~ COOP", "COOP ~ MASS"),
  "H1_COOP→FS_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP", "FS ~ COOP")
)

model_patterns = c(model_patterns, model_patterns_with_MASS)
#model_patterns = c(model_patterns, model_patterns_with_MASS, model_patterns_with_TERR_response)

# Generate all the models
models_list <- lapply(names(model_patterns), function(model_name) {
  create_formula_set(model_patterns[[model_name]], var_map)
})
names(models_list) <- names(model_patterns)

# Create the final model set
models <- do.call(define_model_set, models_list)

# Verify the models are correctly defined
print(models)


result <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'lambda',
  na.rm = TRUE
)

(s <- summary(result))

plot(s)

# prep to plot maps
phylopath_map_positions = data.frame("name" = c(female_song_var, coop_breeding_var, territoriality_var, mass_var), "x" = c(5, 9, 5, 1), "y" = c(4, 7, 1, 7) )


(best_model <- best(result))

(best_model_plot <- plot(best_model, text_size = 4, manual_layout = phylopath_map_positions, labels = labels) +  coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
    ))
ave_model_full_plot <- plot(average(result, cut_off = 2, avg_method = "full"), text_size = 4, manual_layout = phylopath_map_positions, labels = labels) +  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
  )
ave_model_conditional_plot <- plot(average(result, cut_off = 2, avg_method = "conditional"), labels = labels, text_size = 4, manual_layout = phylopath_map_positions) +  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 20, r = 70, b = 20, l = 50, unit = "pt")  
  )

# Create file name
pdf_file <- paste0("phylopath_best_average_models ", 
                   female_song_var, " ", coop_breeding_var, " ", 
                   territoriality_var, " ", mass_var, " ", 
                   Sys.Date(), ".pdf")

# Open PDF device
pdf(file = pdf_file, width = 6, height = 6)  # Square pages

# Print each plot to a new page
print(best_model_plot)
print(ave_model_full_plot)
print(ave_model_conditional_plot)

# Close PDF device
dev.off()



CICsub2_models = s$model[which(s$delta_CICc<2)]

for (i in 1:length(CICsub2_models)) {
  tempmodel = CICsub2_models[i]
  tempmodel_title <- gsub("→", "to", tempmodel)
  tempplot <- plot(choice(result, choice = tempmodel), text_size = 3, manual_layout = phylopath_map_positions)
  tempplot_title <- tempplot + ggtitle(tempmodel_title)
  print(tempplot_title)
}


plot(choice(result, choice = "H1_MASS→COOP_COOP→FS_FS→TERR"), text_size = 3)
plot(choice(result, choice = "F0_FS→TERR_TERR→COOP"), text_size = 3)
plot(choice(result, choice = "H1_COOP→FS_FS→TERR"), text_size = 3)
plot(choice(result, choice = "F0_FS→TERR_TERR→COOP_MASS→FS"), text_size = 3)
plot(choice(result, choice = "H1_COOP→FS_FS→TERR_MASS→TERR"), text_size = 3)
plot(choice(result, choice = "H1_COOP→FS_FS→TERR_COOP→TERR"), text_size = 3)
plot(choice(result, choice = "H1_COOP→FS_FS→TERR_MASS→COOP"), text_size = 3)




