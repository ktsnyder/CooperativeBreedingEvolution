### PhyloGLM - test coop vs each sociality type on same species subset
### 

source("subsettreedata.R")
source("phyloglm_framework/phyloglm_unified_forClaude.R")

data = read.csv("Data_R_2025-07-23.csv")
data$caretakers_normalized <- scale(data$caretakers_Griesser2023, center = T, scale = T)
table(data$social_bonds_Griesser2023)
data$social_bonds_ordinal <- NA
data$social_bonds_ordinal[which(data$social_bonds_Griesser2023 == "a-short")] <- 1
data$social_bonds_ordinal[which(data$social_bonds_Griesser2023 == "b-season")] <- 2
data$social_bonds_ordinal[which(data$social_bonds_Griesser2023 == "c-long")] <- 3
data$social_bonds_normalized <- scale(data$social_bonds_ordinal, center = T, scale = T)
#write.csv(data, "Data_R_2025-07-23_wNormCaretakers_SocBonds.csv", row.names = F)

#data = read.csv("Data_R_2025-07-23_wNormCaretakers_SocBonds.csv")
data$grouping_ordinal = NA
data$grouping_ordinal[which(data$grouping_Griesser2023 == "asocial")] <- 1
data$grouping_ordinal[which(data$grouping_Griesser2023 == "pair")] <- 2
data$grouping_ordinal[which(data$grouping_Griesser2023 == "small_groups")] <- 3
data$grouping_ordinal[which(data$grouping_Griesser2023 == "large_groups")] <- 4
data$grouping_normalized = scale(data$grouping_ordinal, center = T, scale = T)
data %>% group_by(grouping_ordinal, grouping_Griesser2023, grouping_normalized) %>% count
#write.csv(data, "Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv", row.names = F)

data = read.csv("Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv")
 
data_CoopFSsub = data[which(!is.na(data$HighConfidence_Coop) & !is.na(data$FemaleSong_Agg01)),]
rownames(data_CoopFSsub) <- data_CoopFSsub$species
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
tree= read.nexus(treefile)
nBoot = 500

socialityMetrics = c( "Griesser2023.Colonial01" , "Griesser2017FamilialLiving", "Final.polygyny", "caretakers_normalized", "social_bonds_normalized", "grouping_normalized", "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LongSocialBonds", "Griesser2023.SeasonOrLongerSocialBonds")
territoriality_vars = c("TerritorialityWeakVsStrong", "Territory_12vs3")

#formula_variation <- "Add_Soc_To_Base_KeepCoop"
formula_variation <- "Replace_Coop_Additive"

# Create output directory
output_dir <- file.path("Outputs/PhyloglmResults/Replace_Coop_With_Soc/", formula_variation)
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
    old_coef_filename <- paste0("coefficients_boot", nBoot, "_", territoriality_var, "_", tempSoc, "_oldCoop_AIC", round(old_aic, 2), ".csv")
    new_coef_filename <- paste0("coefficients_boot", nBoot, "_", territoriality_var, "_", tempSoc, "_newSoc_AIC", round(new_aic, 2), ".csv")
    
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


### test each sociality factor, without cooperative breeding ----

# Source framework components
source("phyloglm_framework/variable_classification.R")
source("phyloglm_framework/formula_builder.R")
source("phyloglm_framework/data_preparation.R")
source("phyloglm_framework/batch_runner.R")
source("phyloglm_framework/batch_runner_helpers.R")
source("phyloglm_framework/config_builder.R")

output_base_dir <- "Outputs/PhyloglmResults"
data_file = "Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv"
tree_file = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
socialityMetrics = c("caretakers_normalized", "social_bonds_normalized", "grouping_normalized", "Griesser2023.Colonial01", "Griesser2017FamilialLiving") # "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes"

run_custom_batch_example <- function(terr_var) {
  
  # Initialize local summary data frame
  function_summary <- data.frame()
  
  # Load data and tree
  data <- read.csv(data_file)
  data$caretakers_normalized <- scale(data$caretakers_Griesser2023, center = T, scale = T)
  tree <- read.nexus(tree_file)
  
  configs <- list()
  
  # Create bidirectional analyses for each CB variable
  for (soc_var in socialityMetrics) {
    # FS -> CB direction
    configs[[paste0("fs_", soc_var)]] <- create_analysis_config(
      response = "FemaleSong_Agg01",
      predictors = c(soc_var, terr_var),
      complexity_levels = c("null", "main", "additive"),
      controls = c("logMass_normalized"),
      name = paste0("FS_vs_", soc_var, "_", terr_var)
    )
  }
  
  # Run batch
  batch_results <- run_phyloglm_batch(
    analysis_configs = configs,
    data = data,
    tree = tree,
    output_dir = output_base_dir,
    parallel = FALSE,
    bootstrap_n = 500
  )
  
  # Extract and process results for comprehensive summary
  # The results are stored under batch_results$results
  if (!is.null(batch_results$results)) {
    for (analysis_name in names(batch_results$results)) {
      analysis <- batch_results$results[[analysis_name]]
      
      if (!is.null(analysis$success) && analysis$success == TRUE) {
        # Extract sociality variable name - need to parse it better
        # The name format is "FS_vs_SOCVAR_TERRVAR"
        soc_var <- gsub(paste0("FS_vs_(.*?)_", terr_var), "\\1", analysis_name)
        
        # Get best model info from comparison
        if (!is.null(analysis$comparison) && !is.null(analysis$comparison$comparison)) {
          best_model_name <- analysis$comparison$comparison$Model[1]
          best_model_aic <- analysis$comparison$comparison$AIC[1]
          best_model_formula <- paste(deparse(analysis$comparison$best_model$formula), collapse = " ")
          nModels_sub2deltaAIC = sum(analysis$comparison$comparison$deltaAIC < 2)
          sub2_models <- paste(analysis$comparison$comparison$Model[which(analysis$comparison$comparison$deltaAIC < 2)], collapse = ", ")
          
          second_best_model_name <- analysis$comparison$comparison$Model[2]
          second_best_model_aic <- analysis$comparison$comparison$AIC[2]
          
          # Get the best model object
          best_model <- analysis$models$models[[best_model_name]]
          
          # Extract coefficient for sociality variable if it exists
          soc_odds_ratio <- NA
          soc_ci_lower <- NA
          soc_ci_upper <- NA
          soc_p_value <- NA
          
          # Look in the coefficients dataframe
          if (!is.null(analysis$effects)) {
            coef_df <- analysis$effects
            # Find rows for the sociality variable (main effect)
            soc_rows <- coef_df[coef_df$Model == best_model_name & coef_df$Parameter == soc_var, ]
            
            if (nrow(soc_rows) > 0) {
              soc_odds_ratio <- soc_rows$OddsRatio[1]
              soc_ci_lower <- soc_rows$OR_CI_lower[1]
              soc_ci_upper <- soc_rows$OR_CI_upper[1]
              soc_p_value <- soc_rows$p_value[1]
            }
          }
          
          # Get number of species
          n_species <- nrow(analysis$prepared_data$data)
          
          # Create summary row
          summary_row <- data.frame(
            territoriality_var = terr_var,
            sociality_var = soc_var,
            n_species = n_species,
            best_model_type = best_model_name,
            best_model_AIC = best_model_aic,
            best_model_formula = best_model_formula,
            soc_odds_ratio = soc_odds_ratio,
            soc_OR_CI_lower = soc_ci_lower,
            soc_OR_CI_upper = soc_ci_upper,
            soc_p_value = soc_p_value,
            n_models_tested = nrow(analysis$comparison$comparison),
            sub2deltaAIC_models = sub2_models,
            second_best_model_name = second_best_model_name,
            second_best_model_aic = second_best_model_aic,
            stringsAsFactors = FALSE
          )
          
          function_summary <- rbind(function_summary, summary_row)
        }
      }
    }
  }
  
  # Return both batch results and summary
  return(list(batch_results = batch_results, summary = function_summary))
}

# Initialize a comprehensive summary data frame
all_results_summary <- data.frame()

# Run analyses for both territoriality variables
cat("\n========== Running analyses with TerritorialityWeakVsStrong ==========\n")
result1 <- run_custom_batch_example(terr_var = "TerritorialityWeakVsStrong")
all_results_summary <- rbind(all_results_summary, result1$summary)

cat("\n========== Running analyses with Territory_12vs3 ==========\n")
result2 <- run_custom_batch_example(terr_var = "Territory_12vs3")
all_results_summary <- rbind(all_results_summary, result2$summary)

# Save comprehensive summary to CSV
summary_output_dir <- "Outputs/PhyloglmResults/Replace_Coop_With_Soc"
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)

# Check if we have any results
if (nrow(all_results_summary) > 0) {
  # Add columns for model selection criteria
  all_results_summary$is_significant <- all_results_summary$soc_p_value < 0.05
  all_results_summary$OR_excludes_1 <- (all_results_summary$soc_OR_CI_lower > 1) | (all_results_summary$soc_OR_CI_upper < 1)
  
  # Sort by territoriality variable and AIC
  all_results_summary <- all_results_summary[order(all_results_summary$territoriality_var, all_results_summary$best_model_AIC), ]
  
  # Save the comprehensive summary
  write.csv(all_results_summary, 
            file.path(summary_output_dir, "Sociality_Models_Best_Summary_wGroupingNorm_Incl2ndBestModelInfo.csv"),
            row.names = FALSE)
  
  # Print summary to console
  cat("\n\n========== ANALYSIS COMPLETE ==========\n")
  cat("Full summary saved to:", file.path(summary_output_dir, "Sociality_Models_Best_Summary.csv"), "\n")
} else {
  cat("\n\nWARNING: No results were collected. Please check if the batch analyses completed successfully.\n")
  cat("The batch results may be saved in individual folders under:", output_base_dir, "\n")
}
