# Comprehensive stepwise expansion for all base models
# Tests all available predictors and creates detailed summaries

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

# Helper function to get model name for labeling
get_model_label <- function(analysis_name) {
  if (analysis_name == "FS_vs_CB_TerrWS_Mass") {
    return("FS~CB*TerrWS+Mass")
  } else if (analysis_name == "FS_vs_CB_Terr3_Mass") {
    return("FS~CB*Terr3+Mass")
  } else if (analysis_name == "CB_vs_FS_TerrWS_Mass") {
    return("CB~FS*TerrWS")
  } else if (analysis_name == "CB_vs_FS_Terr3_Mass") {
    return("CB~FS*Terr3")
  }
  return(analysis_name)
}

run_comprehensive_stepwise_expansion <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_base_dir = "Outputs/PhyloglmResults/stepwise_comprehensive",
  aic_threshold = 2
) {
  
  # Create base directory
  dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data once
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Analyses to run
  base_analyses <- c("FS_vs_CB_TerrWS_Mass", "FS_vs_CB_Terr3_Mass", 
                     "CB_vs_FS_TerrWS_Mass", "CB_vs_FS_Terr3_Mass")
  
  # Master results storage
  all_expansion_results <- list()
  
  # Run for each base model
  for (analysis_name in base_analyses) {
    
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("BASE ANALYSIS:", analysis_name, "\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
    
    # Create output directory for this analysis
    output_dir <- file.path(output_base_dir, analysis_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get base model
    base_result <- all_results[[analysis_name]]
    if (is.null(base_result)) {
      cat("WARNING: Analysis", analysis_name, "not found in results. Skipping.\n")
      next
    }
    
    best_model_name <- base_result$comparison$comparison$Model[1]
    best_model <- base_result$models$models[[best_model_name]]
    base_aic <- base_result$comparison$comparison$AIC[1]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    
    cat("Base model name:", best_model_name, "\n")
    cat("Base AIC:", round(base_aic, 2), "\n")
    cat("Base formula:", deparse(formula(best_model)), "\n")
    cat("N species:", nrow(base_data), "\n\n")
    
    # Store results
    result <- list(
      analysis_name = analysis_name,
      base_model_name = best_model_name,
      base_model = best_model,
      base_aic = base_aic,
      base_formula = deparse(formula(best_model)),
      base_data = base_data,
      base_tree = base_tree,
      base_n_species = nrow(base_data),
      improvements = data.frame(),
      best_single_addition = NULL
    )
    
    # Test each predictor individually
    cat("Testing predictors:\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    
    # 1. Territory as numeric (if not already in as Terr3)
    if (!grepl("Terr3", analysis_name)) {
      cat("\n1. Territory (1-3 numeric):\n")
      territory_vals <- as.numeric(full_data$Territory[match(base_data$species, full_data$species)])
      test_data <- base_data
      test_data$Territory_num <- territory_vals
      test_data <- test_data[!is.na(test_data$Territory_num), ]
      
      if (nrow(test_data) > 100) {
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        new_formula <- update(formula(best_model), ~ . + Territory_num)
        
        tryCatch({
          new_model <- phyloglm(
            formula = new_formula,
            data = test_data,
            phy = test_tree,
            method = "logistic_MPLE",
            btol = 50,
            log.alpha.bound = 4
          )
          
          new_aic <- -2 * new_model$logLik + 2 * new_model$d
          improvement <- base_aic - new_aic
          
          cat("   Variable added: Territory_num\n")
          cat("   N species:", nrow(test_data), "\n")
          cat("   New formula:", deparse(new_formula), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Territory",
            Variable = "Territory_num",
            Formula_Addition = "+ Territory_num",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Territory",
              variable = "Territory_num",
              improvement = improvement,
              model = new_model,
              data = test_data,
              tree = test_tree,
              formula = new_formula,
              aic = new_aic
            )
          }
          
        }, error = function(e) {
          cat("   Error:", e$message, "\n")
        })
      }
    }
    
    # 2. Migration as numeric
    cat("\n2. Migration (1-3 numeric):\n")
    migration_vals <- as.numeric(full_data$Migration_AVONET[match(base_data$species, full_data$species)])
    test_data <- base_data
    test_data$Migration_num <- migration_vals
    test_data <- test_data[!is.na(test_data$Migration_num), ]
    
    if (nrow(test_data) > 100) {
      test_tree <- keep.tip(base_tree, test_data$species)
      test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
      
      new_formula <- update(formula(best_model), ~ . + Migration_num)
      
      tryCatch({
        new_model <- phyloglm(
          formula = new_formula,
          data = test_data,
          phy = test_tree,
          method = "logistic_MPLE",
          btol = 50,
          log.alpha.bound = 4
        )
        
        new_aic <- -2 * new_model$logLik + 2 * new_model$d
        improvement <- base_aic - new_aic
        
        cat("   Variable added: Migration_num\n")
        cat("   N species:", nrow(test_data), "\n")
        cat("   New formula:", deparse(new_formula), "\n")
        cat("   AIC improvement:", round(improvement, 2), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "Migration",
          Variable = "Migration_num",
          Formula_Addition = "+ Migration_num",
          N_Species = nrow(test_data),
          Base_AIC = base_aic,
          New_AIC = new_aic,
          AIC_Improvement = improvement,
          stringsAsFactors = FALSE
        ))
        
        if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- list(
            predictor = "Migration",
            variable = "Migration_num",
            improvement = improvement,
            model = new_model,
            data = test_data,
            tree = test_tree,
            formula = new_formula,
            aic = new_aic
          )
        }
        
      }, error = function(e) {
        cat("   Error:", e$message, "\n")
      })
    }
    
    # 3. Geographic Region (if not already in model)
    if (!grepl("Region", analysis_name)) {
      cat("\n3. Geographic Region (categorical):\n")
      region_vals <- full_data$GeographicRegion_Jetz[match(base_data$species, full_data$species)]
      test_data <- base_data
      test_data$GeographicRegion_Jetz <- region_vals
      test_data <- test_data[!is.na(test_data$GeographicRegion_Jetz) & test_data$GeographicRegion_Jetz != "", ]
      
      if (nrow(test_data) > 100) {
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        new_formula <- update(formula(best_model), ~ . + GeographicRegion_Jetz)
        
        tryCatch({
          new_model <- phyloglm(
            formula = new_formula,
            data = test_data,
            phy = test_tree,
            method = "logistic_MPLE",
            btol = 50,
            log.alpha.bound = 4
          )
          
          new_aic <- -2 * new_model$logLik + 2 * new_model$d
          improvement <- base_aic - new_aic
          
          cat("   Variable added: GeographicRegion_Jetz\n")
          cat("   N species:", nrow(test_data), "\n")
          cat("   New formula:", deparse(new_formula), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Region",
            Variable = "GeographicRegion_Jetz",
            Formula_Addition = "+ GeographicRegion_Jetz",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Region",
              variable = "GeographicRegion_Jetz",
              improvement = improvement,
              model = new_model,
              data = test_data,
              tree = test_tree,
              formula = new_formula,
              aic = new_aic
            )
          }
          
        }, error = function(e) {
          cat("   Error:", e$message, "\n")
        })
      } else {
        cat("   Not enough species with region data (", nrow(test_data), ")\n")
      }
    }
    
    # 4. Absolute Latitude (if not already in model)
    if (!grepl("absLat|Latitude", analysis_name)) {
      cat("\n4. Absolute Latitude (continuous):\n")
      lat_vals <- abs(full_data$Centroid.Latitude_AVONET[match(base_data$species, full_data$species)])
      test_data <- base_data
      test_data$abs_Latitude <- lat_vals
      test_data <- test_data[!is.na(test_data$abs_Latitude), ]
      
      if (nrow(test_data) > 100) {
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        new_formula <- update(formula(best_model), ~ . + abs_Latitude)
        
        tryCatch({
          new_model <- phyloglm(
            formula = new_formula,
            data = test_data,
            phy = test_tree,
            method = "logistic_MPLE",
            btol = 50,
            log.alpha.bound = 4
          )
          
          new_aic <- -2 * new_model$logLik + 2 * new_model$d
          improvement <- base_aic - new_aic
          
          cat("   Variable added: abs_Latitude\n")
          cat("   N species:", nrow(test_data), "\n")
          cat("   New formula:", deparse(new_formula), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Latitude",
            Variable = "abs_Latitude",
            Formula_Addition = "+ abs_Latitude",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Latitude",
              variable = "abs_Latitude",
              improvement = improvement,
              model = new_model,
              data = test_data,
              tree = test_tree,
              formula = new_formula,
              aic = new_aic
            )
          }
          
        }, error = function(e) {
          cat("   Error:", e$message, "\n")
        })
      }
    }
    
    # 5. Plumage Dimorphism (if not already in model)
    if (!grepl("PlumDim|Plumage", analysis_name)) {
      cat("\n5. Plumage Dimorphism (continuous):\n")
      plum_vals <- full_data$logMaleFemalePlumageDiffAbs[match(base_data$species, full_data$species)]
      test_data <- base_data
      test_data$logMaleFemalePlumageDiffAbs <- plum_vals
      test_data <- test_data[!is.na(test_data$logMaleFemalePlumageDiffAbs), ]
      
      if (nrow(test_data) > 100) {
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        new_formula <- update(formula(best_model), ~ . + logMaleFemalePlumageDiffAbs)
        
        tryCatch({
          new_model <- phyloglm(
            formula = new_formula,
            data = test_data,
            phy = test_tree,
            method = "logistic_MPLE",
            btol = 50,
            log.alpha.bound = 4
          )
          
          new_aic <- -2 * new_model$logLik + 2 * new_model$d
          improvement <- base_aic - new_aic
          
          cat("   Variable added: logMaleFemalePlumageDiffAbs\n")
          cat("   N species:", nrow(test_data), "\n")
          cat("   New formula:", deparse(new_formula), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Plumage",
            Variable = "logMaleFemalePlumageDiffAbs",
            Formula_Addition = "+ logMaleFemalePlumageDiffAbs",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Plumage",
              variable = "logMaleFemalePlumageDiffAbs",
              improvement = improvement,
              model = new_model,
              data = test_data,
              tree = test_tree,
              formula = new_formula,
              aic = new_aic
            )
          }
          
        }, error = function(e) {
          cat("   Error:", e$message, "\n")
        })
      } else {
        cat("   Not enough species with plumage data (", nrow(test_data), ")\n")
      }
    }
    
    # 6. Wing Dimorphism (if not already in model)
    if (!grepl("WingDim|Wing", analysis_name)) {
      cat("\n6. Wing Dimorphism (continuous):\n")
      wing_vals <- full_data$PercentAbsLogWingDimorphism[match(base_data$species, full_data$species)]
      test_data <- base_data
      test_data$PercentAbsLogWingDimorphism <- wing_vals
      test_data <- test_data[!is.na(test_data$PercentAbsLogWingDimorphism), ]
      
      if (nrow(test_data) > 100) {
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        new_formula <- update(formula(best_model), ~ . + PercentAbsLogWingDimorphism)
        
        tryCatch({
          new_model <- phyloglm(
            formula = new_formula,
            data = test_data,
            phy = test_tree,
            method = "logistic_MPLE",
            btol = 50,
            log.alpha.bound = 4
          )
          
          new_aic <- -2 * new_model$logLik + 2 * new_model$d
          improvement <- base_aic - new_aic
          
          cat("   Variable added: PercentAbsLogWingDimorphism\n")
          cat("   N species:", nrow(test_data), "\n")
          cat("   New formula:", deparse(new_formula), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Wing",
            Variable = "PercentAbsLogWingDimorphism",
            Formula_Addition = "+ PercentAbsLogWingDimorphism",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Wing",
              variable = "PercentAbsLogWingDimorphism",
              improvement = improvement,
              model = new_model,
              data = test_data,
              tree = test_tree,
              formula = new_formula,
              aic = new_aic
            )
          }
          
        }, error = function(e) {
          cat("   Error:", e$message, "\n")
        })
      } else {
        cat("   Not enough species with wing data (", nrow(test_data), ")\n")
      }
    }
    
    # 7. Familial Living
    cat("\n7. Familial Living (binary):\n")
    fam_vals <- full_data$Griesser2017FamilialLiving[match(base_data$species, full_data$species)]
    test_data <- base_data
    test_data$Griesser2017FamilialLiving <- fam_vals
    test_data <- test_data[!is.na(test_data$Griesser2017FamilialLiving), ]
    
    if (nrow(test_data) > 100) {
      test_tree <- keep.tip(base_tree, test_data$species)
      test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
      
      new_formula <- update(formula(best_model), ~ . + Griesser2017FamilialLiving)
      
      tryCatch({
        new_model <- phyloglm(
          formula = new_formula,
          data = test_data,
          phy = test_tree,
          method = "logistic_MPLE",
          btol = 50,
          log.alpha.bound = 4
        )
        
        new_aic <- -2 * new_model$logLik + 2 * new_model$d
        improvement <- base_aic - new_aic
        
        cat("   Variable added: Griesser2017FamilialLiving\n")
        cat("   N species:", nrow(test_data), "\n")
        cat("   New formula:", deparse(new_formula), "\n")
        cat("   AIC improvement:", round(improvement, 2), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "FamilialLiving",
          Variable = "Griesser2017FamilialLiving",
          Formula_Addition = "+ Griesser2017FamilialLiving",
          N_Species = nrow(test_data),
          Base_AIC = base_aic,
          New_AIC = new_aic,
          AIC_Improvement = improvement,
          stringsAsFactors = FALSE
        ))
        
        if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- list(
            predictor = "FamilialLiving",
            variable = "Griesser2017FamilialLiving",
            improvement = improvement,
            model = new_model,
            data = test_data,
            tree = test_tree,
            formula = new_formula,
            aic = new_aic
          )
        }
        
      }, error = function(e) {
        cat("   Error:", e$message, "\n")
      })
    } else {
      cat("   Not enough species with familial living data (", nrow(test_data), ")\n")
    }
    
    # Show best single predictor
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      cat("\n", paste(rep("-", 50), collapse=""), "\n")
      cat("BEST SINGLE PREDICTOR:", result$best_single_addition$predictor, "\n")
      cat("Variable name:", result$best_single_addition$variable, "\n")
      cat("AIC improvement:", round(result$best_single_addition$improvement, 2), "\n")
      
      # Extract coefficients
      coef_summary <- summary(result$best_single_addition$model)$coefficients
      final_coefs <- data.frame(
        Parameter = rownames(coef_summary),
        Estimate = coef_summary[, "Estimate"],
        StdErr = coef_summary[, "StdErr"],
        z_value = coef_summary[, "z.value"],
        P.Value = coef_summary[, "p.value"],
        stringsAsFactors = FALSE
      )
      
      # Calculate effects
      final_effects <- calculate_effect_sizes(final_coefs)
      
      # Store results
      result$final_model <- result$best_single_addition$model
      result$final_data <- result$best_single_addition$data
      result$final_tree <- result$best_single_addition$tree
      result$final_formula <- deparse(result$best_single_addition$formula)
      result$final_aic <- result$best_single_addition$aic
      result$final_effects <- final_effects
      
      # Show significant effects
      cat("\nSignificant effects (p < 0.05):\n")
      sig_effects <- final_effects[final_effects$P.Value < 0.05, ]
      for (i in 1:nrow(sig_effects)) {
        cat("  ", sig_effects$Parameter[i], ": OR =", 
            round(sig_effects$OddsRatio[i], 2),
            "(", round(sig_effects$OR_CI_lower[i], 2), "-",
            round(sig_effects$OR_CI_upper[i], 2), "),",
            "p =", format(sig_effects$P.Value[i], digits = 3), "\n")
      }
    } else {
      cat("\nNo predictors improved the model by >", aic_threshold, "AIC units\n")
    }
    
    # Save individual results
    saveRDS(result, file.path(output_dir, "expansion_result.rds"))
    
    # Create individual visualizations
    create_individual_plots(result, output_dir, aic_threshold)
    
    # Create detailed summary
    create_detailed_summary(result, output_dir, aic_threshold)
    
    # Store in master results
    all_expansion_results[[analysis_name]] <- result
  }
  
  # Create master comparison plots
  create_master_plots(all_expansion_results, output_base_dir, aic_threshold)
  
  # Create master summary
  create_master_summary(all_expansion_results, output_base_dir, aic_threshold)
  
  cat("\n\nAll results saved to:", output_base_dir, "\n")
  
  return(all_expansion_results)
}

# Function to create individual plots
create_individual_plots <- function(result, output_dir, aic_threshold) {
  
  # 1. Improvement barplot
  if (nrow(result$improvements) > 0) {
    p_improve <- ggplot(result$improvements, 
                       aes(x = reorder(Predictor, AIC_Improvement), 
                           y = AIC_Improvement)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      labs(title = paste("Stepwise Improvements for", result$analysis_name),
           subtitle = paste("Base model:", result$base_model_name, 
                           "| Base AIC:", round(result$base_aic, 1)),
           x = "", y = "AIC Improvement",
           caption = paste("Red line = threshold (", aic_threshold, ")", sep="")) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(output_dir, paste0(result$analysis_name, "_improvements.png")),
           p_improve, width = 8, height = 6, dpi = 300)
  }
  
  # 2. Forest plot if we have a best model
  if (!is.null(result$final_effects)) {
    effects <- result$final_effects
    
    # Clean names
    effects$Variable <- effects$Parameter
    effects$Variable <- gsub("HighConfidence_Coop", "Cooperative Breeding", effects$Variable)
    effects$Variable <- gsub("FemaleSong_Agg01", "Female Song", effects$Variable)
    effects$Variable <- gsub("TerritorialityWeakVsStrong", "Territoriality", effects$Variable)
    effects$Variable <- gsub("Territory12vs3", "Territory (1,2 vs 3)", effects$Variable)
    effects$Variable <- gsub("Strong", " (Strong)", effects$Variable)
    effects$Variable <- gsub("Territory_num", "Territory (1-3)", effects$Variable)
    effects$Variable <- gsub("Migration_num", "Migration (1-3)", effects$Variable)
    effects$Variable <- gsub("abs_Latitude", "Latitude (absolute)", effects$Variable)
    effects$Variable <- gsub("logMass_AVONET", "Body Mass (log)", effects$Variable)
    effects$Variable <- gsub("Griesser2017FamilialLiving", "Familial Living", effects$Variable)
    effects$Variable <- gsub("logMaleFemalePlumageDiffAbs", "Plumage Dimorphism (log)", effects$Variable)
    effects$Variable <- gsub("PercentAbsLogWingDimorphism", "Wing Dimorphism (%)", effects$Variable)
    effects$Variable <- gsub(":", " × ", effects$Variable)
    
    p_forest <- ggplot(effects, aes(x = OddsRatio, y = reorder(Variable, OddsRatio))) +
      geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
      geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), 
                    height = 0.2, color = "gray30") +
      geom_point(aes(color = P.Value < 0.05), size = 4) +
      scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8)) +
      scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray60"),
                        labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
      labs(
        title = paste("Best Model for", result$analysis_name),
        subtitle = paste(result$base_model_name, "+", result$best_single_addition$predictor),
        x = "Odds Ratio (95% CI)",
        y = "",
        color = "Significant"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom"
      )
    
    ggsave(file.path(output_dir, paste0(result$analysis_name, "_forest.png")),
           p_forest, width = 8, height = 6, dpi = 300)
  }
}

# Function to create detailed summary
create_detailed_summary <- function(result, output_dir, aic_threshold) {
  
  summary_text <- paste("STEPWISE EXPANSION SUMMARY FOR", result$analysis_name, "\n")
  summary_text <- paste0(summary_text, paste(rep("=", 70), collapse=""), "\n\n")
  
  summary_text <- paste0(summary_text, "BASE MODEL INFORMATION\n")
  summary_text <- paste0(summary_text, paste(rep("-", 30), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Model name: ", result$base_model_name, "\n")
  summary_text <- paste0(summary_text, "Formula: ", result$base_formula, "\n")
  summary_text <- paste0(summary_text, "Base AIC: ", round(result$base_aic, 2), "\n")
  summary_text <- paste0(summary_text, "N species: ", result$base_n_species, "\n\n")
  
  summary_text <- paste0(summary_text, "PREDICTORS TESTED\n")
  summary_text <- paste0(summary_text, paste(rep("-", 30), collapse=""), "\n")
  
  if (nrow(result$improvements) > 0) {
    for (i in 1:nrow(result$improvements)) {
      imp <- result$improvements[i,]
      summary_text <- paste0(summary_text, "\n", i, ". ", imp$Predictor, "\n")
      summary_text <- paste0(summary_text, "   Variable: ", imp$Variable, "\n")
      summary_text <- paste0(summary_text, "   Formula addition: ", imp$Formula_Addition, "\n")
      summary_text <- paste0(summary_text, "   N species: ", imp$N_Species, "\n")
      summary_text <- paste0(summary_text, "   New AIC: ", round(imp$New_AIC, 2), "\n")
      summary_text <- paste0(summary_text, "   AIC improvement: ", round(imp$AIC_Improvement, 2))
      if (imp$AIC_Improvement > aic_threshold) {
        summary_text <- paste0(summary_text, " **EXCEEDS THRESHOLD**")
      }
      summary_text <- paste0(summary_text, "\n")
    }
  }
  
  if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
    summary_text <- paste0(summary_text, "\n\nBEST MODEL\n")
    summary_text <- paste0(summary_text, paste(rep("-", 30), collapse=""), "\n")
    summary_text <- paste0(summary_text, "Best predictor: ", result$best_single_addition$predictor, "\n")
    summary_text <- paste0(summary_text, "Variable: ", result$best_single_addition$variable, "\n")
    summary_text <- paste0(summary_text, "Final formula: ", result$final_formula, "\n")
    summary_text <- paste0(summary_text, "Final AIC: ", round(result$final_aic, 2), "\n")
    summary_text <- paste0(summary_text, "Total improvement: ", 
                          round(result$best_single_addition$improvement, 2), "\n\n")
    
    summary_text <- paste0(summary_text, "COEFFICIENTS (BEST MODEL)\n")
    summary_text <- paste0(summary_text, paste(rep("-", 30), collapse=""), "\n")
    
    if (!is.null(result$final_effects)) {
      effects <- result$final_effects
      for (i in 1:nrow(effects)) {
        summary_text <- paste0(summary_text, effects$Parameter[i], ": ")
        summary_text <- paste0(summary_text, "OR = ", round(effects$OddsRatio[i], 3))
        summary_text <- paste0(summary_text, " (", round(effects$OR_CI_lower[i], 3), 
                              "-", round(effects$OR_CI_upper[i], 3), ")")
        summary_text <- paste0(summary_text, ", p = ", format(effects$P.Value[i], digits = 3))
        if (effects$P.Value[i] < 0.05) {
          summary_text <- paste0(summary_text, " *")
        }
        summary_text <- paste0(summary_text, "\n")
      }
    }
  } else {
    summary_text <- paste0(summary_text, "\n\nNo predictors improved model by >", 
                          aic_threshold, " AIC units\n")
  }
  
  writeLines(summary_text, file.path(output_dir, paste0(result$analysis_name, "_summary.txt")))
}

# Function to create master plots
create_master_plots <- function(all_results, output_dir, aic_threshold) {
  
  # Combine all improvements
  all_improvements <- data.frame()
  for (name in names(all_results)) {
    imp <- all_results[[name]]$improvements
    if (nrow(imp) > 0) {
      imp$Analysis <- name
      imp$Model_Label <- get_model_label(name)
      all_improvements <- rbind(all_improvements, imp)
    }
  }
  
  if (nrow(all_improvements) > 0) {
    # Master improvement plot
    p_master <- ggplot(all_improvements, 
                      aes(x = reorder(Predictor, AIC_Improvement), 
                          y = AIC_Improvement, 
                          fill = Model_Label)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Comprehensive Stepwise Model Improvements",
           subtitle = "All base models and predictors tested",
           x = "", y = "AIC Improvement",
           fill = "Base Model",
           caption = paste("Red line = threshold (", aic_threshold, ")", sep="")) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 14))
    
    ggsave(file.path(output_dir, "master_improvements_comparison.png"),
           p_master, width = 12, height = 8, dpi = 300)
  }
}

# Function to create master summary
create_master_summary <- function(all_results, output_dir, aic_threshold) {
  
  summary_text <- "COMPREHENSIVE STEPWISE EXPANSION SUMMARY\n"
  summary_text <- paste0(summary_text, paste(rep("=", 70), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Generated: ", Sys.Date(), "\n")
  summary_text <- paste0(summary_text, "AIC improvement threshold: ", aic_threshold, "\n\n")
  
  # Summary table
  summary_text <- paste0(summary_text, "OVERVIEW OF ALL ANALYSES\n")
  summary_text <- paste0(summary_text, paste(rep("-", 70), collapse=""), "\n\n")
  
  for (name in names(all_results)) {
    result <- all_results[[name]]
    summary_text <- paste0(summary_text, name, "\n")
    summary_text <- paste0(summary_text, "  Base model: ", result$base_model_name, "\n")
    summary_text <- paste0(summary_text, "  Base AIC: ", round(result$base_aic, 2), "\n")
    
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      summary_text <- paste0(summary_text, "  Best addition: ", 
                            result$best_single_addition$predictor,
                            " (AIC improvement = ", 
                            round(result$best_single_addition$improvement, 2), ")\n")
    } else {
      summary_text <- paste0(summary_text, "  No significant improvements found\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  # Detailed results by predictor
  summary_text <- paste0(summary_text, "\nDETAILED RESULTS BY PREDICTOR\n")
  summary_text <- paste0(summary_text, paste(rep("-", 70), collapse=""), "\n")
  
  predictors <- c("Territory", "Migration", "Region", "Latitude", "Plumage", "Wing", "FamilialLiving")
  
  for (pred in predictors) {
    summary_text <- paste0(summary_text, "\n", pred, ":\n")
    
    for (name in names(all_results)) {
      result <- all_results[[name]]
      imp_row <- result$improvements[result$improvements$Predictor == pred, ]
      
      if (nrow(imp_row) > 0) {
        summary_text <- paste0(summary_text, "  ", name, ": ")
        summary_text <- paste0(summary_text, "AIC improvement = ", 
                              round(imp_row$AIC_Improvement[1], 2))
        if (imp_row$AIC_Improvement[1] > aic_threshold) {
          summary_text <- paste0(summary_text, " **")
        }
        summary_text <- paste0(summary_text, " (", imp_row$N_Species[1], " species)\n")
      }
    }
  }
  
  summary_text <- paste0(summary_text, "\n** = Exceeds threshold of ", aic_threshold, " AIC units\n")
  
  writeLines(summary_text, file.path(output_dir, "master_summary.txt"))
}

# Run the comprehensive analysis
results <- run_comprehensive_stepwise_expansion()