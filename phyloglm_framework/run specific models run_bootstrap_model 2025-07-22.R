## PhyloGLM get effect sizes of best models (individual models)
## All outputs to Outputs/PhyloglmResults/Single_Model_Runs3 on 7/22/2025

setwd("Desktop/CooperativeBreedingEvolution/")

#### Run single analyses for results tables ----
source("phyloglm_framework/phyloglm_unified (2).R")
source("phyloglm_framework/phyloglm_unified_forClaude.R")
source("subsettreedata.R")

##### Best base models -----
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong"), newdata = "Data_R_2025-07-23.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CBxTerrWS_Mass <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized"), data = subsetdf, tree = subsettree, n_boot = 500, 
                    method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_Mass",
                    save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs_MassNorm", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

CB_vs_FSxTerrWS <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * TerritorialityWeakVsStrong"), data = subsetdf, tree = subsettree, n_boot = 500, 
                    method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerrWS",
                    save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs_MassNorm", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3"), newdata = "Data_R_2025-07-23.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr3_Mass2 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_Terr3_MassNorm", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs_MassNorm", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

FS_vs_Terr3_Mass1 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ Territory_12vs3 + logMass_AVONET"), data = subsetdf, tree = subsettree, n_boot = 500,  method = "logistic_MPLE", save_prefix = "FS_vs_Terr3_MassNonNorm1", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

CB_vs_FSxTerr3_Mass <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * Territory_12vs3"), data = subsetdf, tree = subsettree, n_boot = 500, 
                        method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerr3",
                        save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs_MassNorm", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


#### Single model runs, AltCoops ----
AltCoops = c("MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnCoop", "CockburnInferred", "BiagoliniCoop", "DaleCoop", "DowningCoop", "JetzCoopInclCockburn", "Griesser2017Coop", "CornwallisCoop")

for (i in 1:length(AltCoops)) {
  tempcoop = AltCoops[i]
  
  subsetout <- subsettreedata(columns = c(tempcoop, "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong"), newdata = "Data_R_2025-07-23.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  subsetdf = subsetout$subsetdf
  rownames(subsetdf) <- subsetdf$species
  subsettree = subsetout$subsettree
  
  FS_response_formula = paste("FemaleSong_Agg01 ~", tempcoop ,"* TerritorialityWeakVsStrong + logMass_normalized")
  print(FS_response_formula)
  
  FS_vs_CBxTerrWS_Mass <- run_bootstrap_model(formula = as.formula(FS_response_formula), data = subsetdf, tree = subsettree, n_boot = 500, 
                                              method = "logistic_MPLE", save_prefix = paste0("FS_vs_", tempcoop, "xTerrWS_Mass"),
                                              save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3/AltCoops", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)
  print(FS_vs_CBxTerrWS_Mass$n_successful_boots)
  print(FS_vs_CBxTerrWS_Mass$n_converged)
  
  CB_response_formula = paste(tempcoop, "~ FemaleSong_Agg01 * TerritorialityWeakVsStrong")
  print(CB_response_formula)
  
  CB_vs_FSxTerrWS <- run_bootstrap_model(formula = as.formula(CB_response_formula), data = subsetdf, tree = subsettree, n_boot = 500, 
                                         method = "logistic_MPLE", save_prefix = paste0(tempcoop,"_vs_FSxTerrWS"),
                                         save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3/AltCoops", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)
  print(CB_vs_FSxTerrWS$n_successful_boots)
  print(CB_vs_FSxTerrWS$n_converged)
  
}

##### Full models - all possible predictor variables -----
# Terr Weak/Strong
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong", "Migration_num", "abs_Latitude_normalized", "PercentAbsLogWingDimorphism_normalized", "logMaleFemalePlumageDiffAbs_normalized", "Griesser2017FamilialLiving", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CBxTerrWS_All <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + Griesser2017FamilialLiving + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
          method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_AllPredictors",
          save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

CB_vs_FSxTerrWS_All <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * TerritorialityWeakVsStrong + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + Griesser2017FamilialLiving + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                       method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerrWS_AllPredictors",
                                       save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


# Terr 12vs3
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "Migration_num", "abs_Latitude_normalized", "PercentAbsLogWingDimorphism_normalized", "logMaleFemalePlumageDiffAbs_normalized", "Griesser2017FamilialLiving", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr3_All <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + Griesser2017FamilialLiving + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_Terr3_AllPredictors", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

CB_vs_FSxTerr3_All <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * Territory_12vs3 + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + Griesser2017FamilialLiving + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                           method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerr3_AllPredictors",
                                           save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


#### Almost-full models for FS as response variable; remove familial living -----
# Terr Weak/Strong
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong", "Migration_num", "abs_Latitude_normalized", "PercentAbsLogWingDimorphism_normalized", "logMaleFemalePlumageDiffAbs_normalized", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CBxTerrWS_All <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + Griesser2017FamilialLiving + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                           method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_AllExceptFamLiv",
                                           save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

# Terr 12vs3
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "Migration_num", "abs_Latitude_normalized", "PercentAbsLogWingDimorphism_normalized", "logMaleFemalePlumageDiffAbs_normalized", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr3_All <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized + logMass_normalized + Migration_num + abs_Latitude_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_Terr3_AllExceptFamLiv", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

#### just grouping_normalized with best model from each TerrWS and TerrYR ----
subsetout <- subsettreedata(columns = c("FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong", "grouping_normalized"), newdata = "Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree
FS_vs_xTerrWS <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ grouping_normalized + TerritorialityWeakVsStrong + logMass_normalized"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                     method = "logistic_MPLE", save_prefix = "FS_vs_GroupNorm_TerrWS_Mass",
                                     save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

subsetout <- subsettreedata(columns = c("FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "grouping_normalized"), newdata = "Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree
FS_vs_xTerrYR <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ grouping_normalized + Territory_12vs3 + logMass_normalized"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                     method = "logistic_MPLE", save_prefix = "FS_vs_GroupNorm_TerrYR_MassNorm",
                                     save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

#### Previous best models from stepwise ----
# Terr Weak/Strong, FS response
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong", "Territory", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

# This one is based on the model initially summarized in supp table 15
FS_vs_CBxTerrWS_Step <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized + Territory + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                           method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_Step",
                                           save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

# With Terr 1-3 removed
FS_vs_CBxTerrWS_StepNoTerr13 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                            method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_StepNoTerr1-3",
                                            save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

# This one is based on the "iterative" stepwise run from 20250715
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "TerritorialityWeakVsStrong", "PercentAbsLogWingDimorphism_normalized", "logMaleFemalePlumageDiffAbs_normalized"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CBxTerrWS_StepIter <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_normalized + PercentAbsLogWingDimorphism_normalized + logMaleFemalePlumageDiffAbs_normalized"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                                    method = "logistic_MPLE", save_prefix = "FS_vs_CBxTerrWS_StepIter20250715",
                                                    save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)



# Terr 12vs3, FS response
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "PercentAbsLogWingDimorphism_normalized", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr3_Step2 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + logMass_normalized + PercentAbsLogWingDimorphism_normalized + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_Terr3_Step2", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


# Terr Weak/Strong, CB response ( best base model was just FemaleSong*TerrWS, only added familial living )
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong", "Griesser2017FamilialLiving"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

CB_vs_FSxTerrWS_StepIter <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * TerritorialityWeakVsStrong + Griesser2017FamilialLiving"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                                method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerrWS_StepIter20250715",
                                                save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


# Terr 12vs3, CB response ( best base model was FemaleSong*Terr3+Mass, only added familial living )
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "Territory_12vs3", "logMass_normalized", "Griesser2017FamilialLiving"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

CB_vs_FSxTerr3_StepIter <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * Territory_12vs3 + logMass_normalized + Griesser2017FamilialLiving"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                                method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerr3_StepIter20250715",
                                                save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

# Terr 12vs3 step models with Territory instead of Territory_12vs3 - CB response
CB_vs_FSxTerr13_StepIter <- run_bootstrap_model(formula = as.formula("HighConfidence_Coop ~ FemaleSong_Agg01 * Territory + logMass_normalized + Griesser2017FamilialLiving"), data = subsetdf, tree = subsettree, n_boot = 500, 
                                               method = "logistic_MPLE", save_prefix = "CB_vs_FSxTerr13_StepIter20250715",
                                               save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)

# Terr 12vs3 step models with Territory instead of Territory_12vs3 - FS response
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "PercentAbsLogWingDimorphism_normalized", "GeographicRegion_Jetz"), newdata = "Data_R_2025-07-21.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr13_Step2 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop + Territory + logMass_normalized + PercentAbsLogWingDimorphism_normalized + GeographicRegion_Jetz"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_Terr13_Step2", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)


# Test why caretakers_normalized is breaking the model run
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01"), newdata = "Data_R_2025-07-23_wNormCaretakers.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree
run_bootstrap_model(formula = as.formula(FemaleSong_Agg01 ~ HighConfidence_Coop*TerritorialityWeakVsStrong +      logMass_normalized + caretakers_normalized), data = subsetdf, tree =  subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CB_TerrWS_caretakers", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE, bias_threshold_sd = 50)
                    

#### familial living * coop breed interaction ----
subsetout <- subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "logMass_normalized", "Territory_12vs3", "Griesser2017FamilialLiving"), newdata = "Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv", newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
subsetdf = subsetout$subsetdf
rownames(subsetdf) <- subsetdf$species
subsettree = subsetout$subsettree

FS_vs_CB_Terr13_Step2 <- run_bootstrap_model(formula = as.formula("FemaleSong_Agg01 ~ HighConfidence_Coop * Griesser2017FamilialLiving + Territory_12vs3"), data = subsetdf, tree = subsettree, n_boot = 500, method = "logistic_MPLE", save_prefix = "FS_vs_CBxFam_Terr13", save_matrices = TRUE, matrix_dir = "Outputs/PhyloglmResults/Single_Model_Runs3", save_coefficient_csv = TRUE, use_bootstrap_pvalues = TRUE)
