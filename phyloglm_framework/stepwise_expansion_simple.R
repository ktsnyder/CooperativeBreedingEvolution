# Simplified stepwise expansion that correctly calculates improvements
# This version tests each predictor against the base model independently

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

run_stepwise_expansion_corrected <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_simple_", format(Sys.Date(), "%Y%m%d"))
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
    base_aic <- base_result$comparison$comparison$AIC[1]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    
    cat("Base model:", best_model_name, "\n")
    cat("Base AIC:", round(base_aic, 2), "\n")
    cat("Base formula:", deparse(formula(best_model)), "\n")
    cat("N species:", nrow(base_data), "\n\n")
    
    # Store results
    result <- list(
      base_model = best_model,
      base_aic = base_aic,
      base_data = base_data,
      base_tree = base_tree,
      improvements = data.frame(),
      best_single_addition = NULL
    )
    
    # Test each predictor individually
    cat("Testing predictors individually:\n")
    cat(paste(rep("-", 40), collapse=""), "\n")
    
    # 1. Territory as numeric
    cat("\n1. Territory (1-3):\n")
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
        
        cat("   N species:", nrow(test_data), "\n")
        cat("   AIC improvement:", round(improvement, 2), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "Territory",
          Variable = "Territory_num",
          N_Species = nrow(test_data),
          Base_AIC = base_aic,
          New_AIC = new_aic,
          AIC_Improvement = improvement,
          stringsAsFactors = FALSE
        ))
        
        # Keep track of best improvement
        if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- list(
            predictor = "Territory",
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
    
    # 2. Migration as numeric
    cat("\n2. Migration (1-3):\n")
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
        
        cat("   N species:", nrow(test_data), "\n")
        cat("   AIC improvement:", round(improvement, 2), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "Migration",
          Variable = "Migration_num",
          N_Species = nrow(test_data),
          Base_AIC = base_aic,
          New_AIC = new_aic,
          AIC_Improvement = improvement,
          stringsAsFactors = FALSE
        ))
        
        if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
          result$best_single_addition <- list(
            predictor = "Migration",
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
    
    # 3. Absolute Latitude
    if (!grepl("absLat", analysis_name)) {
      cat("\n3. Absolute Latitude:\n")
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
          
          cat("   N species:", nrow(test_data), "\n")
          cat("   AIC improvement:", round(improvement, 2), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Latitude",
            Variable = "abs_Latitude",
            N_Species = nrow(test_data),
            Base_AIC = base_aic,
            New_AIC = new_aic,
            AIC_Improvement = improvement,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- list(
              predictor = "Latitude",
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
      result$final_data <- result$best_single_addition$data
      result$final_tree <- result$best_single_addition$tree
      result$final_formula <- result$best_single_addition$formula
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
    
    expansion_results[[analysis_name]] <- result
  }
  
  # Save results
  saveRDS(expansion_results, file.path(output_dir, "expansion_results.rds"))
  
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
    p_improve <- ggplot(all_improvements, 
                       aes(x = reorder(Predictor, AIC_Improvement), 
                           y = AIC_Improvement, 
                           fill = Analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_brewer(palette = "Set1", 
                     labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding",
                               "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song")) +
    labs(title = "Stepwise Model Improvements",
         subtitle = paste("AIC improvement threshold =", aic_threshold),
         x = "", y = "AIC Improvement",
         fill = "Direction") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
    ggsave(file.path(output_dir, "improvements_comparison.png"),
           p_improve, width = 10, height = 6, dpi = 300)
  }
  
  # 2. Forest plots for expanded models
  forest_plots <- list()
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    
    if ("final_effects" %in% names(result)) {
      effects <- result$final_effects
      
      # Clean names
      effects$Variable <- effects$Parameter
      effects$Variable <- gsub("HighConfidence_Coop", "Cooperative Breeding", effects$Variable)
      effects$Variable <- gsub("FemaleSong_Agg01", "Female Song", effects$Variable)
      effects$Variable <- gsub("TerritorialityWeakVsStrong", "Territoriality", effects$Variable)
      effects$Variable <- gsub("Strong", " (Strong)", effects$Variable)
      effects$Variable <- gsub("Territory_num", "Territory (1-3)", effects$Variable)
      effects$Variable <- gsub("Migration_num", "Migration (1-3)", effects$Variable)
      effects$Variable <- gsub("abs_Latitude", "Latitude (absolute)", effects$Variable)
      effects$Variable <- gsub("logMass_AVONET", "Body Mass (log)", effects$Variable)
      effects$Variable <- gsub(":", " × ", effects$Variable)
      
      # Create plot
      p <- ggplot(effects, aes(x = OddsRatio, y = reorder(Variable, OddsRatio))) +
        geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
        geom_errorbarh(aes(xmin = OR_CI_lower, xmax = OR_CI_upper), 
                      height = 0.2, color = "gray30") +
        geom_point(aes(color = P.Value < 0.05), size = 4) +
        scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8)) +
        scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray60"),
                          labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
        labs(
          title = ifelse(grepl("FS_vs_CB", name),
                        "Female Song → Cooperative Breeding",
                        "Cooperative Breeding → Female Song"),
          subtitle = paste("Best model: Base +", result$best_single_addition$predictor),
          x = "Odds Ratio (95% CI)",
          y = "",
          color = "Significant"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom"
        )
      
      forest_plots[[name]] <- p
    }
  }
  
  if (length(forest_plots) > 0) {
    combined_forest <- wrap_plots(forest_plots, ncol = 1)
    ggsave(file.path(output_dir, "expanded_models_forest.png"),
           combined_forest, width = 10, height = 6 * length(forest_plots), dpi = 300)
  }
  
  # 3. Summary table
  summary_text <- "STEPWISE EXPANSION SUMMARY\n"
  summary_text <- paste0(summary_text, paste(rep("=", 60), collapse=""), "\n\n")
  
  for (name in names(expansion_results)) {
    result <- expansion_results[[name]]
    summary_text <- paste0(summary_text, name, "\n")
    summary_text <- paste0(summary_text, paste(rep("-", nchar(name)), collapse=""), "\n")
    summary_text <- paste0(summary_text, "Base AIC: ", round(result$base_aic, 2), "\n")
    
    # Show all improvements
    if (nrow(result$improvements) > 0) {
      summary_text <- paste0(summary_text, "\nPredictors tested:\n")
      for (i in 1:nrow(result$improvements)) {
        imp <- result$improvements[i,]
        summary_text <- paste0(summary_text, "  ", imp$Predictor, 
                              ": AIC improvement = ", round(imp$AIC_Improvement, 2),
                              " (", imp$N_Species, " species)\n")
      }
    }
    
    if ("final_aic" %in% names(result)) {
      summary_text <- paste0(summary_text, "\nBest single addition: ", 
                            result$best_single_addition$predictor, "\n")
      summary_text <- paste0(summary_text, "Final AIC: ", round(result$final_aic, 2), "\n")
      summary_text <- paste0(summary_text, "Total improvement: ", 
                            round(result$best_single_addition$improvement, 2), "\n")
    } else {
      summary_text <- paste0(summary_text, "No predictors improved model by >", 
                            aic_threshold, " AIC units\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "expansion_summary.txt"))
  cat("\n", summary_text)
  
  cat("\nResults saved to:", output_dir, "\n")
  
  return(expansion_results)
}

# Run the analysis
results <- run_stepwise_expansion_corrected()