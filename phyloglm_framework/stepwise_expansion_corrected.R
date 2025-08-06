# Corrected stepwise expansion that properly handles different sample sizes
# This version refits the base model to the same subset of species for each predictor

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

run_stepwise_expansion_proper <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_corrected_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Results storage
  expansion_results <- list()
  
  # Analyze each direction
  for (analysis_name in c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass")) {
    
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("Analyzing:", analysis_name, "\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n")
    
    # Get base model
    base_result <- all_results[[analysis_name]]
    best_model_name <- base_result$comparison$comparison$Model[1]
    best_model <- base_result$models$models[[best_model_name]]
    original_base_aic <- base_result$comparison$comparison$AIC[1]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    base_formula <- formula(best_model)
    
    cat("Base model:", best_model_name, "\n")
    cat("Original AIC (", nrow(base_data), " species):", round(original_base_aic, 2), "\n")
    cat("Base formula:", deparse(base_formula), "\n\n")
    
    # Store results
    result <- list(
      base_model = best_model,
      original_base_aic = original_base_aic,
      original_n_species = nrow(base_data),
      base_data = base_data,
      base_tree = base_tree,
      improvements = data.frame(),
      best_single_addition = NULL
    )
    
    # Test each predictor individually
    cat("Testing predictors (with proper subsetting):\n")
    cat(paste(rep("-", 40), collapse=""), "\n")
    
    # Helper function to test a predictor
    test_predictor <- function(predictor_name, predictor_values, var_name) {
      cat("\n", predictor_name, ":\n", sep="")
      
      # Add predictor to base data
      test_data <- base_data
      test_data[[var_name]] <- predictor_values[match(test_data$species, full_data$species)]
      
      # Remove NAs
      test_data <- test_data[!is.na(test_data[[var_name]]), ]
      
      if (nrow(test_data) < 100) {
        cat("   Not enough species with data (", nrow(test_data), ")\n")
        return(NULL)
      }
      
      # Match tree
      test_tree <- keep.tip(base_tree, test_data$species)
      test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
      
      cat("   N species with", predictor_name, "data:", nrow(test_data), "\n")
      
      # CRITICAL: Refit base model to this subset
      cat("   Refitting base model to same species subset...\n")
      base_subset_model <- phyloglm(
        formula = base_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      base_subset_aic <- -2 * base_subset_model$logLik + 2 * base_subset_model$d
      cat("   Base model AIC on subset:", round(base_subset_aic, 2), "\n")
      
      # Now fit model with predictor
      new_formula <- update(base_formula, paste("~ . +", var_name))
      new_model <- phyloglm(
        formula = new_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      new_aic <- -2 * new_model$logLik + 2 * new_model$d
      
      # Calculate PROPER improvement
      improvement <- base_subset_aic - new_aic
      
      cat("   Model with", predictor_name, "AIC:", round(new_aic, 2), "\n")
      cat("   TRUE AIC improvement:", round(improvement, 2), "\n")
      
      return(list(
        predictor = predictor_name,
        variable = var_name,
        n_species = nrow(test_data),
        base_subset_aic = base_subset_aic,
        new_aic = new_aic,
        improvement = improvement,
        model = new_model,
        data = test_data,
        tree = test_tree,
        formula = new_formula
      ))
    }
    
    # 1. Territory as numeric
    if (!grepl("Terr3", analysis_name)) {
      territory_result <- test_predictor(
        "Territory (1-3)", 
        as.numeric(full_data$Territory),
        "Territory_num"
      )
      
      if (!is.null(territory_result)) {
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = territory_result$predictor,
          Variable = territory_result$variable,
          N_Species = territory_result$n_species,
          Base_Subset_AIC = territory_result$base_subset_aic,
          New_AIC = territory_result$new_aic,
          AIC_Improvement = territory_result$improvement,
          stringsAsFactors = FALSE
        ))
        
        if (is.null(result$best_single_addition) || 
            territory_result$improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- territory_result
        }
      }
    }
    
    # 2. Migration as numeric
    migration_result <- test_predictor(
      "Migration (1-3)",
      as.numeric(full_data$Migration_AVONET),
      "Migration_num"
    )
    
    if (!is.null(migration_result)) {
      result$improvements <- rbind(result$improvements, data.frame(
        Predictor = migration_result$predictor,
        Variable = migration_result$variable,
        N_Species = migration_result$n_species,
        Base_Subset_AIC = migration_result$base_subset_aic,
        New_AIC = migration_result$new_aic,
        AIC_Improvement = migration_result$improvement,
        stringsAsFactors = FALSE
      ))
      
      if (is.null(result$best_single_addition) || 
          migration_result$improvement > result$best_single_addition$improvement) {
        result$best_single_addition <- migration_result
      }
    }
    
    # 3. Absolute Latitude
    if (!grepl("absLat|Latitude", analysis_name)) {
      latitude_result <- test_predictor(
        "Absolute Latitude",
        abs(full_data$Centroid.Latitude_AVONET),
        "abs_Latitude"
      )
      
      if (!is.null(latitude_result)) {
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = latitude_result$predictor,
          Variable = latitude_result$variable,
          N_Species = latitude_result$n_species,
          Base_Subset_AIC = latitude_result$base_subset_aic,
          New_AIC = latitude_result$new_aic,
          AIC_Improvement = latitude_result$improvement,
          stringsAsFactors = FALSE
        ))
        
        if (is.null(result$best_single_addition) || 
            latitude_result$improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- latitude_result
        }
      }
    }
    
    # 4. Wing Dimorphism
    wing_result <- test_predictor(
      "Wing Dimorphism",
      full_data$PercentAbsLogWingDimorphism,
      "PercentAbsLogWingDimorphism"
    )
    
    if (!is.null(wing_result)) {
      result$improvements <- rbind(result$improvements, data.frame(
        Predictor = wing_result$predictor,
        Variable = wing_result$variable,
        N_Species = wing_result$n_species,
        Base_Subset_AIC = wing_result$base_subset_aic,
        New_AIC = wing_result$new_aic,
        AIC_Improvement = wing_result$improvement,
        stringsAsFactors = FALSE
      ))
      
      if (is.null(result$best_single_addition) || 
          wing_result$improvement > result$best_single_addition$improvement) {
        result$best_single_addition <- wing_result
      }
    }
    
    # 5. Familial Living
    familial_result <- test_predictor(
      "Familial Living",
      full_data$Griesser2017FamilialLiving,
      "Griesser2017FamilialLiving"
    )
    
    if (!is.null(familial_result)) {
      result$improvements <- rbind(result$improvements, data.frame(
        Predictor = familial_result$predictor,
        Variable = familial_result$variable,
        N_Species = familial_result$n_species,
        Base_Subset_AIC = familial_result$base_subset_aic,
        New_AIC = familial_result$new_aic,
        AIC_Improvement = familial_result$improvement,
        stringsAsFactors = FALSE
      ))
      
      if (is.null(result$best_single_addition) || 
          familial_result$improvement > result$best_single_addition$improvement) {
        result$best_single_addition <- familial_result
      }
    }
    
    # Show best single predictor
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      cat("\n", paste(rep("-", 40), collapse=""), "\n")
      cat("Best single predictor:", result$best_single_addition$predictor, "\n")
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
      result$final_data = result$best_single_addition$data
      result$final_tree = result$best_single_addition$tree
      result$final_formula = result$best_single_addition$formula
      result$final_aic = result$best_single_addition$new_aic
      result$final_effects = final_effects
      
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
    
    expansion_results[[analysis_name]] <- result
  }
  
  # Save results
  saveRDS(expansion_results, file.path(output_dir, "expansion_results_corrected.rds"))
  
  # Create visualizations
  cat("\n\nCreating visualizations...\n")
  
  # 1. Improvement comparison plot
  all_improvements <- data.frame()
  for (name in names(expansion_results)) {
    imp <- expansion_results[[name]]$improvements
    if (nrow(imp) > 0) {
      imp$Analysis <- name
      all_improvements <- rbind(all_improvements, imp)
    }
  }
  
  if (nrow(all_improvements) > 0) {
    # Add percentage of species retained
    all_improvements$Percent_Species <- 100 * all_improvements$N_Species / 
      ifelse(all_improvements$Analysis == "FS_vs_CB_TerrWS_Mass", 875, 875)
    
    p_improve <- ggplot(all_improvements, 
                       aes(x = reorder(Predictor, AIC_Improvement), 
                           y = AIC_Improvement, 
                           fill = Analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, color = "black") +
    # Add text labels with N species
    geom_text(aes(label = paste0("n=", N_Species)), 
              position = position_dodge(width = 0.9),
              hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_brewer(palette = "Set1", 
                     labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding",
                               "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song")) +
    labs(title = "Corrected Stepwise Model Improvements",
         subtitle = paste("AIC improvements calculated on same species subsets\nThreshold =", aic_threshold),
         x = "", y = "AIC Improvement",
         fill = "Direction") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
    ggsave(file.path(output_dir, "improvements_comparison_corrected.png"),
           p_improve, width = 10, height = 6, dpi = 300)
    
    # Create table showing the correction
    comparison_table <- all_improvements
    comparison_table$Percent_Retained <- round(comparison_table$Percent_Species, 1)
    comparison_table <- comparison_table[order(-comparison_table$AIC_Improvement), 
                                       c("Analysis", "Predictor", "N_Species", 
                                         "Percent_Retained", "AIC_Improvement")]
    
    write.csv(comparison_table, 
              file.path(output_dir, "improvement_comparison_table.csv"),
              row.names = FALSE)
  }
  
  # 3. Summary table
  summary_text <- "CORRECTED STEPWISE EXPANSION SUMMARY\n"
  summary_text <- paste0(summary_text, paste(rep("=", 60), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Generated: ", Sys.Date(), "\n")
  summary_text <- paste0(summary_text, "AIC improvement threshold: ", aic_threshold, "\n\n")
  
  summary_text <- paste0(summary_text, "IMPORTANT: AIC improvements now calculated correctly\n")
  summary_text <- paste0(summary_text, "by refitting base model to same species subset as expanded model\n\n")
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    summary_text <- paste0(summary_text, name, "\n")
    summary_text <- paste0(summary_text, paste(rep("-", nchar(name)), collapse=""), "\n")
    summary_text <- paste0(summary_text, "Original base AIC (", result$original_n_species, 
                          " species): ", round(result$original_base_aic, 2), "\n\n")
    
    if (nrow(result$improvements) > 0) {
      summary_text <- paste0(summary_text, "Predictors tested:\n")
      for (i in 1:nrow(result$improvements)) {
        imp <- result$improvements[i,]
        summary_text <- paste0(summary_text, "  ", imp$Predictor, 
                              ": AIC improvement = ", round(imp$AIC_Improvement, 2),
                              " (", imp$N_Species, " species)\n")
      }
    }
    
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      summary_text <- paste0(summary_text, "\nBest single addition: ", 
                            result$best_single_addition$predictor, "\n")
      summary_text <- paste0(summary_text, "AIC improvement: ", 
                            round(result$best_single_addition$improvement, 2), "\n")
    } else {
      summary_text <- paste0(summary_text, "\nNo predictors improved model by >", 
                            aic_threshold, " AIC units\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "expansion_summary_corrected.txt"))
  cat("\n", summary_text)
  
  cat("\nCorrected results saved to:", output_dir, "\n")
  
  return(expansion_results)
}

# Run the corrected analysis
results <- run_stepwise_expansion_proper()