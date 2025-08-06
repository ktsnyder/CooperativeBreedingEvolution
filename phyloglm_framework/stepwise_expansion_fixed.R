# Fixed stepwise expansion with bootstrap support
# This version correctly handles the model building and visualization

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load batch runner helpers for effect calculations
source("phyloglm_framework/batch_runner_helpers.R")

#' Run bootstrap for a phyloglm model
#' @param model phyloglm model object
#' @param data data frame used to fit the model
#' @param tree phylogenetic tree
#' @param n_boot number of bootstrap iterations
#' @return matrix of bootstrap coefficients
bootstrap_model <- function(model, data, tree, n_boot = 100) {
  
  formula <- formula(model)
  n_species <- nrow(data)
  
  # Initialize storage
  boot_coefs <- matrix(NA, nrow = n_boot, ncol = length(coef(model)))
  colnames(boot_coefs) <- names(coef(model))
  
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  
  for (i in 1:n_boot) {
    # Resample species with replacement
    boot_idx <- sample(1:n_species, n_species, replace = TRUE)
    boot_data <- data[boot_idx, ]
    boot_tree <- keep.tip(tree, unique(rownames(boot_data)))
    
    # Match data to tree
    boot_data <- boot_data[rownames(boot_data) %in% boot_tree$tip.label, ]
    boot_data <- boot_data[match(boot_tree$tip.label, rownames(boot_data)), ]
    
    # Fit model
    tryCatch({
      boot_model <- phyloglm(
        formula = formula,
        data = boot_data,
        phy = boot_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      boot_coefs[i, ] <- coef(boot_model)
    }, error = function(e) {
      # Leave as NA if failed
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(boot_coefs)
}

#' Calculate bootstrap confidence intervals
#' @param boot_coefs matrix of bootstrap coefficients
#' @param conf_level confidence level (default 0.95)
#' @return data frame with CIs
calculate_bootstrap_ci <- function(boot_coefs, conf_level = 0.95) {
  # Remove failed bootstraps
  boot_coefs <- boot_coefs[complete.cases(boot_coefs), ]
  
  if (nrow(boot_coefs) < 10) {
    warning("Less than 10 successful bootstrap iterations")
  }
  
  # Calculate percentile CIs
  alpha <- 1 - conf_level
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2
  
  ci_df <- data.frame(
    Parameter = colnames(boot_coefs),
    Boot_Mean = colMeans(boot_coefs),
    Boot_SD = apply(boot_coefs, 2, sd),
    CI_lower = apply(boot_coefs, 2, quantile, probs = lower_q),
    CI_upper = apply(boot_coefs, 2, quantile, probs = upper_q),
    n_boot_success = nrow(boot_coefs),
    stringsAsFactors = FALSE
  )
  
  return(ci_df)
}

#' Main stepwise expansion function
run_stepwise_expansion_complete <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  n_bootstrap = 100,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_complete_", format(Sys.Date(), "%Y%m%d"))
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
      expanded_models = list()
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
        cat("   Status:", ifelse(improvement > aic_threshold, "KEPT", "Rejected"), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "Territory",
          Variable = "Territory_num",
          N_Species = nrow(test_data),
          AIC_Improvement = improvement,
          Kept = improvement > aic_threshold,
          stringsAsFactors = FALSE
        ))
        
        if (improvement > aic_threshold) {
          result$expanded_models$Territory <- list(
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
        cat("   Status:", ifelse(improvement > aic_threshold, "KEPT", "Rejected"), "\n")
        
        result$improvements <- rbind(result$improvements, data.frame(
          Predictor = "Migration",
          Variable = "Migration_num",
          N_Species = nrow(test_data),
          AIC_Improvement = improvement,
          Kept = improvement > aic_threshold,
          stringsAsFactors = FALSE
        ))
        
        if (improvement > aic_threshold) {
          result$expanded_models$Migration <- list(
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
          cat("   Status:", ifelse(improvement > aic_threshold, "KEPT", "Rejected"), "\n")
          
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = "Latitude",
            Variable = "abs_Latitude",
            N_Species = nrow(test_data),
            AIC_Improvement = improvement,
            Kept = improvement > aic_threshold,
            stringsAsFactors = FALSE
          ))
          
          if (improvement > aic_threshold) {
            result$expanded_models$Latitude <- list(
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
    
    # Build final model with all kept predictors
    kept_predictors <- result$improvements$Predictor[result$improvements$Kept]
    
    if (length(kept_predictors) > 0) {
      cat("\n", paste(rep("-", 40), collapse=""), "\n")
      cat("Building final model with:", paste(kept_predictors, collapse=", "), "\n")
      
      # Use the data from the best individual model
      if (length(kept_predictors) == 1) {
        # If only one predictor, use its data directly
        best_individual <- result$expanded_models[[kept_predictors[1]]]
        final_data <- best_individual$data
        final_tree <- best_individual$tree
        final_formula <- best_individual$formula
      } else {
        # If multiple predictors, start with the first and add others
        best_individual <- result$expanded_models[[kept_predictors[1]]]
        final_data <- best_individual$data
        final_tree <- best_individual$tree
        
        # Add other predictors
        for (i in 2:length(kept_predictors)) {
          pred <- kept_predictors[i]
          if (pred == "Territory" && !("Territory_num" %in% names(final_data))) {
            final_data$Territory_num <- as.numeric(full_data$Territory[match(final_data$species, full_data$species)])
          } else if (pred == "Migration" && !("Migration_num" %in% names(final_data))) {
            final_data$Migration_num <- as.numeric(full_data$Migration_AVONET[match(final_data$species, full_data$species)])
          } else if (pred == "Latitude" && !("abs_Latitude" %in% names(final_data))) {
            final_data$abs_Latitude <- abs(full_data$Centroid.Latitude_AVONET[match(final_data$species, full_data$species)])
          }
        }
      }
      
      if (length(kept_predictors) > 1) {
        # For multiple predictors, we need to filter data
        cat("Variables in final_data:", paste(names(final_data), collapse=", "), "\n")
        
        # Get the new variables we added
        new_vars <- c()
        if ("Territory" %in% kept_predictors) new_vars <- c(new_vars, "Territory_num")
        if ("Migration" %in% kept_predictors) new_vars <- c(new_vars, "Migration_num")
        if ("Latitude" %in% kept_predictors) new_vars <- c(new_vars, "abs_Latitude")
        
        # Check for NAs in new variables
        for (var in new_vars) {
          n_na <- sum(is.na(final_data[[var]]))
          if (n_na > 0) {
            cat("  Removing", n_na, "species with NA in", var, "\n")
            final_data <- final_data[!is.na(final_data[[var]]), ]
          }
        }
        
        cat("N species after removing NAs:", nrow(final_data), "\n")
        
        if (nrow(final_data) == 0) {
          cat("ERROR: No complete cases in final data!\n")
          next
        }
        
        # Update tree to match
        final_tree <- keep.tip(final_tree, final_data$species)
        final_data <- final_data[match(final_tree$tip.label, final_data$species), ]
        
        # Build formula by updating the base formula
        additions <- result$improvements$Variable[result$improvements$Kept]
        final_formula <- formula(best_model)
        for (var in additions) {
          final_formula <- update(final_formula, as.formula(paste("~ . +", var)))
        }
      }
      
      cat("Final formula:", deparse(final_formula), "\n")
      cat("N species:", nrow(final_data), "\n")
      
      # Fit final model
      final_model <- phyloglm(
        formula = final_formula,
        data = final_data,
        phy = final_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
      
      final_aic <- -2 * final_model$logLik + 2 * final_model$d
      total_improvement <- base_aic - final_aic
      
      cat("Final AIC:", round(final_aic, 2), "\n")
      cat("Total improvement:", round(total_improvement, 2), "\n")
      
      # Extract coefficients
      coef_summary <- summary(final_model)$coefficients
      final_coefs <- data.frame(
        Parameter = rownames(coef_summary),
        Estimate = coef_summary[, "Estimate"],
        StdErr = coef_summary[, "StdErr"],
        z_value = coef_summary[, "z.value"],
        P.Value = coef_summary[, "p.value"],
        stringsAsFactors = FALSE
      )
      
      # Bootstrap if requested
      if (n_bootstrap > 0) {
        cat("\nRunning bootstrap (", n_bootstrap, "iterations)...\n")
        boot_coefs <- bootstrap_model(final_model, final_data, final_tree, n_bootstrap)
        boot_ci <- calculate_bootstrap_ci(boot_coefs)
        
        # Merge with coefficients
        final_coefs <- merge(final_coefs, 
                           boot_ci[, c("Parameter", "CI_lower", "CI_upper")],
                           by = "Parameter", all.x = TRUE)
        
        cat("Bootstrap complete:", boot_ci$n_boot_success[1], "successful iterations\n")
      }
      
      # Calculate effects
      final_effects <- calculate_effect_sizes(final_coefs)
      
      # Store results
      result$final_model <- final_model
      result$final_data <- final_data
      result$final_tree <- final_tree
      result$final_formula <- final_formula
      result$final_aic <- final_aic
      result$total_improvement <- total_improvement
      result$final_coefficients <- final_coefs
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
    # Ensure we have all required columns
    if (!("Predictor" %in% names(all_improvements))) {
      cat("WARNING: No 'Predictor' column found in improvements data\n")
      cat("Columns found:", paste(names(all_improvements), collapse=", "), "\n")
    } else {
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
          subtitle = paste("Expanded model with", 
                          paste(result$improvements$Predictor[result$improvements$Kept], 
                                collapse = ", ")),
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
    
    if ("final_aic" %in% names(result)) {
      summary_text <- paste0(summary_text, "Final AIC: ", round(result$final_aic, 2), "\n")
      summary_text <- paste0(summary_text, "Total improvement: ", 
                            round(result$total_improvement, 2), "\n")
      summary_text <- paste0(summary_text, "Predictors added: ",
                            paste(result$improvements$Predictor[result$improvements$Kept],
                                  collapse = ", "), "\n")
    } else {
      summary_text <- paste0(summary_text, "No predictors improved model\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "expansion_summary.txt"))
  cat("\n", summary_text)
  
  cat("\nResults saved to:", output_dir, "\n")
  
  return(expansion_results)
}

# Run the analysis
results <- run_stepwise_expansion_complete()