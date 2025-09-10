### PhyloGLM - test coop vs each sociality type on same species subset
### 

source("subsettreedata.R")
source(file.path("PhyloGLM_functions", "batch_runner_helpers.R"))

data = read.csv("Data_R.csv")
 
data_CoopFSsub = data[which(!is.na(data$HighConfidence_Coop) & !is.na(data$FemaleSong_Agg01)),]
rownames(data_CoopFSsub) <- data_CoopFSsub$species
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
tree= read.nexus(treefile)

if (!exists("nBoot")) {
  nBoot = 100
  print("Defaulting to nBoot = 100. To run the test as performed in the manuscript, set nBoot = 500.")
}


socialityMetrics = c( "Griesser2023.Colonial01" , "Griesser2017FamilialLiving", "caretakers_normalized", "social_bonds_normalized", "grouping_normalized")
territoriality_vars = c("TerritorialityWeakVsStrong", "Territory_12vs3")

#formula_variation <- "Add_Soc_To_Base_KeepCoop"
formula_variation <- "Replace_Coop_Additive"

# Create output directory
output_dir <- file.path("Outputs","PhyloGLM_outputs","Replace_Coop_With_Soc", formula_variation)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize summary table
summary_results <- data.frame()

for (i in 1:length(socialityMetrics)) {
  tempSoc = socialityMetrics[i]
  
  for (j in 1:2) {
    territoriality_var = territoriality_vars[j]
    subsetout <- subsettreedata(columns = c(tempSoc, territoriality_var, "logMass_normalized"), newdata = data_CoopFSsub, newtree = tree)
    subsetdf = subsetout$subsetdf
    subsettree = subsetout$subsettree
    
    nSpecies = length(subsetdf$species)
    nCoop = sum(subsetdf$HighConfidence_Coop == 1)
    nFSpresent = sum(subsetdf$FemaleSong_Agg01 == 1)
    nBoth = sum(subsetdf$HighConfidence_Coop == 1 & subsetdf$FemaleSong_Agg01 == 1, na.rm = T)
    percentCoop = nCoop/nSpecies
    percentFSpresent = nFSpresent/nSpecies
    percentBoth = nBoth/nSpecies
    
    if (territoriality_var == "TerritorialityWeakVsStrong") {
      old_formula_text = "FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop + logMass_normalized"
    } else if (territoriality_var == "Territory_12vs3") {
      old_formula_text = "FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop + logMass_normalized"
    }
    
    if (formula_variation == "Add_Soc_To_Base_KeepCoop") {
      new_formula_text = paste(old_formula_text, "+", tempSoc)  
    } else {
      new_formula_text = gsub("HighConfidence_Coop", tempSoc, old_formula_text) 
    }
    
    old_formula <- as.formula(old_formula_text)
    new_formula <- as.formula(new_formula_text)
    
    # Run bootstrap for old formula (with HighConfidence_Coop)
    boot_results_oldform <- run_bootstrap_model(
      formula = old_formula,
      data = subsetdf,
      tree = subsettree,
      n_boot = nBoot,
      save_prefix = paste0(territoriality_var, "_", tempSoc, "_oldCoop"),
      save_matrices = FALSE,
      matrix_dir = output_dir,
      save_coefficient_csv = FALSE,  # We'll save manually with AIC in filename
      use_bootstrap_pvalues = TRUE
    )
    
    # Run bootstrap for new formula (with sociality metric)
    boot_results_newform <- run_bootstrap_model(
      formula = new_formula,
      data = subsetdf,
      tree = subsettree,
      n_boot = nBoot,
      save_prefix = paste0(territoriality_var, "_", tempSoc, "_newSoc"),
      save_matrices = FALSE,
      matrix_dir = output_dir,
      save_coefficient_csv = FALSE,  # We'll save manually with AIC in filename
      use_bootstrap_pvalues = TRUE
    )
    
    old_deg_freedom = boot_results_oldform$fit$d
    new_deg_freedom = boot_results_newform$fit$d
    
    # Calculate AIC for both models
    old_aic_standard <- -2 * boot_results_oldform$fit$logLik + 2 * boot_results_oldform$fit$d
    new_aic_standard <- -2 * boot_results_newform$fit$logLik + 2 * boot_results_newform$fit$d
    
    old_aic_penalized <- -2 * boot_results_oldform$fit$penlogLik + 2 * boot_results_oldform$fit$d
    new_aic_penalized <- -2 * boot_results_newform$fit$penlogLik + 2 * boot_results_newform$fit$d
    
    # get phylogenetic signal for both models
    old_alpha = boot_results_oldform$fit$alpha
    new_alpha = boot_results_newform$fit$alpha
    
    # Fit null model (intercept only)
    null_model <- phyloglm(FemaleSong_Agg01 ~ 1, 
                           data = subsetdf, 
                           phy = subsettree,
                           method = "logistic_MPLE")  # or your method
    
    logLik_null <- null_model$logLik
    
    # Then calculate pseudo-R²
    pseudo_R2_old_CB <- 1 - (boot_results_oldform$fit$logLik / logLik_null)
    pseudo_R2_new_soc <- 1 - (boot_results_newform$fit$logLik / logLik_null)
    
    
    # Save coefficient CSVs with AIC in filename
    old_coef_filename <- paste0("coefficients_boot", nBoot, "_", territoriality_var, "_", tempSoc, "_oldCoop_AIC", round(old_aic_standard, 2), ".csv")
    new_coef_filename <- paste0("coefficients_boot", nBoot, "_", territoriality_var, "_", tempSoc, "_newSoc_AIC", round(new_aic_standard, 2), ".csv")
    
    write.csv(boot_results_oldform$coefficients, 
              file.path(output_dir, old_coef_filename), 
              row.names = FALSE)
    write.csv(boot_results_newform$coefficients, 
              file.path(output_dir, new_coef_filename), 
              row.names = FALSE)
    
    # Extract statistics for HighConfidence_Coop from old model
    coop_row <- boot_results_oldform$coefficients[boot_results_oldform$coefficients$Parameter == "HighConfidence_Coop", ]
    
    # Extract statistics for sociality variable from new model
    soc_row <- boot_results_newform$coefficients[boot_results_newform$coefficients$Parameter == tempSoc, ]
    
    # Create summary row
    summary_row <- data.frame(
      territoriality_var = territoriality_var,
      sociality_var = tempSoc,
      nSpecies = nSpecies,
      nCoop = nCoop,
      nFSpresent = nFSpresent,
      nBoth = nBoth,
      percentCoop = percentCoop,
      percentFSpresent = percentFSpresent,
      percentBoth = percentBoth,
      n_boot = nBoot,
      old_model_formula = old_formula_text,
      old_model_AIC_standard = old_aic_standard,
      old_model_AIC_penalized = old_aic_penalized,
      pseudo_R2_old_CB = pseudo_R2_old_CB,
      old_deg_freedom = old_deg_freedom,
      old_phy_signal = old_alpha,
      old_Coop_Odds_Ratio = ifelse(nrow(coop_row) > 0, coop_row$Odds_Ratio, NA),
      old_Coop_OR_CI_Lower = ifelse(nrow(coop_row) > 0, coop_row$OR_CI_Lower, NA),
      old_Coop_OR_CI_Upper = ifelse(nrow(coop_row) > 0, coop_row$OR_CI_Upper, NA),
      old_Coop_p_value = ifelse(nrow(coop_row) > 0, coop_row$p_value, NA),
      new_model_formula = new_formula_text,
      new_model_AIC_standard = new_aic_standard,
      new_model_AIC_penalized = new_aic_penalized,
      pseudo_R2_new_soc = pseudo_R2_new_soc,
      new_deg_freedom = new_deg_freedom,
      new_phy_signal = new_alpha,
      new_Soc_Odds_Ratio = ifelse(nrow(soc_row) > 0, soc_row$Odds_Ratio, NA),
      new_Soc_OR_CI_Lower = ifelse(nrow(soc_row) > 0, soc_row$OR_CI_Lower, NA),
      new_Soc_OR_CI_Upper = ifelse(nrow(soc_row) > 0, soc_row$OR_CI_Upper, NA),
      new_Soc_p_value = ifelse(nrow(soc_row) > 0, soc_row$p_value, NA),
      AIC_standard_difference = old_aic_standard - new_aic_standard,  # Positive means new model is better
      AIC_penalized_difference = old_aic_penalized - new_aic_penalized,
      OR_Ratio = coop_row$Odds_Ratio/soc_row$Odds_Ratio,
      stringsAsFactors = FALSE
    )
    
    # Add to summary table
    summary_results <- rbind(summary_results, summary_row)
    
    # Print progress
    cat("Completed:", territoriality_var, "-", tempSoc, 
        "| nSpecies:", nSpecies, "| nCoop:", nCoop,
        "| Old AIC:", round(old_aic_standard, 2), "| New AIC:", round(new_aic_standard, 2), "\n")
  }
}

# Save summary table
write.csv(summary_results, 
          file.path(output_dir, paste0("Replace_Coop_Summary_Table_boot", nBoot,".csv")), 
          row.names = FALSE)

cat("\nAnalysis complete. Results saved to:", output_dir, "\n")


