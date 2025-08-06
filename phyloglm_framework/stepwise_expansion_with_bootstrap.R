# Stepwise expansion with bootstrap confidence intervals
# This version runs bootstrap iterations for each model to get proper CIs

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)
library(parallel)

source("phyloglm_framework/batch_runner_helpers.R")

# Function to run bootstrap for a single model
run_bootstrap_model <- function(formula, data, tree, n_boot = 1000, 
                               method = "logistic_MPLE", ncores = 1) {
  
  cat("      Running", n_boot, "bootstrap iterations...\n")
  
  # Fit original model
  original_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4
  )
  
  # Extract original coefficients
  original_coef <- coef(original_fit)
  n_params <- length(original_coef)
  
  # Bootstrap function
  boot_fun <- function(i) {
    tryCatch({
      boot_fit <- phyloglm(
        formula = formula,
        data = data,
        phy = tree,
        method = method,
        btol = 50,
        log.alpha.bound = 4,
        boot = 1
      )
      return(coef(boot_fit))
    }, error = function(e) {
      return(rep(NA, n_params))
    })
  }
  
  # Run bootstrap (simplified - no parallelization for now)
  boot_results <- lapply(1:n_boot, boot_fun)
  
  # Convert to matrix
  boot_matrix <- do.call(rbind, boot_results)
  
  # Calculate statistics
  boot_means <- colMeans(boot_matrix, na.rm = TRUE)
  boot_sds <- apply(boot_matrix, 2, sd, na.rm = TRUE)
  boot_lower <- apply(boot_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
  boot_upper <- apply(boot_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  # Create summary
  coef_summary <- data.frame(
    Parameter = names(original_coef),
    Estimate = original_coef,
    Boot_Mean = boot_means,
    Boot_SD = boot_sds,
    CI_Lower = boot_lower,
    CI_Upper = boot_upper,
    stringsAsFactors = FALSE
  )
  
  return(list(
    fit = original_fit,
    coefficients = coef_summary,
    bootstrap_matrix = boot_matrix,
    n_successful_boots = sum(!is.na(boot_matrix[,1]))
  ))
}

run_stepwise_expansion_bootstrap <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2,
  n_bootstrap = 1000,
  ncores = parallel::detectCores() - 1,
  test_mode = FALSE  # If TRUE, only run 10 bootstraps for testing
) {
  
  if (test_mode) {
    n_bootstrap <- 10
    cat("*** TEST MODE: Only running", n_bootstrap, "bootstraps ***\n\n")
  }
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_bootstrap_", format(Sys.Date(), "%Y%m%d"))
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
      predictor_results = list(),
      best_single_addition = NULL
    )
    
    # Test each predictor individually
    cat("Testing predictors with bootstrap CIs:\n")
    cat(paste(rep("-", 40), collapse=""), "\n")
    
    # Helper function to test a predictor with bootstrap
    test_predictor_bootstrap <- function(predictor_name, predictor_values, var_name) {
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
      
      # Refit base model to this subset WITH BOOTSTRAP
      cat("   Fitting base model on subset with bootstrap...\n")
      base_subset_boot <- run_bootstrap_model(
        formula = base_formula,
        data = test_data,
        tree = test_tree,
        n_boot = n_bootstrap,
        ncores = ncores
      )
      base_subset_fit <- base_subset_boot$fit
      base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
      cat("   Base model AIC on subset:", round(base_subset_aic, 2), "\n")
      cat("   Bootstrap successful for", base_subset_boot$n_successful_boots, "/", n_bootstrap, "iterations\n")
      
      # Now fit model with predictor WITH BOOTSTRAP
      new_formula <- update(base_formula, paste("~ . +", var_name))
      cat("   Fitting expanded model with bootstrap...\n")
      new_boot <- run_bootstrap_model(
        formula = new_formula,
        data = test_data,
        tree = test_tree,
        n_boot = n_bootstrap,
        ncores = ncores
      )
      new_fit <- new_boot$fit
      new_aic <- -2 * new_fit$logLik + 2 * new_fit$d
      
      # Calculate improvement
      improvement <- base_subset_aic - new_aic
      
      cat("   Model with", predictor_name, "AIC:", round(new_aic, 2), "\n")
      cat("   AIC improvement:", round(improvement, 2), "\n")
      cat("   Bootstrap successful for", new_boot$n_successful_boots, "/", n_bootstrap, "iterations\n")
      
      # Extract main effect changes
      response_var <- all.vars(base_formula)[1]
      if (response_var == "FemaleSong_Agg01") {
        main_pred <- "HighConfidence_Coop"
      } else {
        main_pred <- "FemaleSong_Agg01"
      }
      
      # Get main effect from both models
      base_main <- base_subset_boot$coefficients[base_subset_boot$coefficients$Parameter == main_pred, ]
      new_main <- new_boot$coefficients[new_boot$coefficients$Parameter == main_pred, ]
      
      main_effect_change <- NA
      main_effect_change_ci <- c(NA, NA)
      
      if (nrow(base_main) > 0 && nrow(new_main) > 0) {
        main_effect_change <- new_main$Estimate - base_main$Estimate
        
        # Calculate CI for the change using bootstrap samples
        if (main_pred %in% colnames(base_subset_boot$bootstrap_matrix) && 
            main_pred %in% colnames(new_boot$bootstrap_matrix)) {
          base_main_idx <- which(colnames(base_subset_boot$bootstrap_matrix) == main_pred)
          new_main_idx <- which(colnames(new_boot$bootstrap_matrix) == main_pred)
          
          change_samples <- new_boot$bootstrap_matrix[, new_main_idx] - 
                           base_subset_boot$bootstrap_matrix[, base_main_idx]
          main_effect_change_ci <- quantile(change_samples, c(0.025, 0.975), na.rm = TRUE)
        }
      }
      
      return(list(
        predictor = predictor_name,
        variable = var_name,
        n_species = nrow(test_data),
        base_subset_aic = base_subset_aic,
        new_aic = new_aic,
        improvement = improvement,
        base_bootstrap = base_subset_boot,
        new_bootstrap = new_boot,
        main_effect_change = main_effect_change,
        main_effect_change_ci = main_effect_change_ci,
        data = test_data,
        tree = test_tree,
        formula = new_formula
      ))
    }
    
    # Test predictors
    predictors_to_test <- list(
      list(name = "Territory (1-3)", var = "Territory_num", 
           values = as.numeric(full_data$Territory), 
           skip_if = grepl("Terr3", analysis_name)),
      list(name = "Migration (1-3)", var = "Migration_num",
           values = as.numeric(full_data$Migration_AVONET), 
           skip_if = FALSE),
      list(name = "Absolute Latitude", var = "abs_Latitude",
           values = abs(full_data$Centroid.Latitude_AVONET), 
           skip_if = grepl("absLat|Latitude", analysis_name)),
      list(name = "Wing Dimorphism", var = "PercentAbsLogWingDimorphism",
           values = full_data$PercentAbsLogWingDimorphism, 
           skip_if = FALSE),
      list(name = "Familial Living", var = "Griesser2017FamilialLiving",
           values = full_data$Griesser2017FamilialLiving, 
           skip_if = FALSE),
      list(name = "Plumage Dimorphism", var = "logMaleFemalePlumageDiffAbs",
           values = log(full_data$MaleFemalePlumageDiffAbs + 1), 
           skip_if = FALSE),
      list(name = "Geographic Region", var = "GeographicRegion_Jetz",
           values = full_data$GeographicRegion_Jetz, 
           skip_if = FALSE)
    )
    
    for (pred in predictors_to_test) {
      if (!pred$skip_if) {
        pred_result <- test_predictor_bootstrap(pred$name, pred$values, pred$var)
        
        if (!is.null(pred_result)) {
          # Store full results
          result$predictor_results[[pred$name]] <- pred_result
          
          # Add to summary table
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor = pred_result$predictor,
            Variable = pred_result$variable,
            N_Species = pred_result$n_species,
            Base_Subset_AIC = pred_result$base_subset_aic,
            New_AIC = pred_result$new_aic,
            AIC_Improvement = pred_result$improvement,
            Main_Effect_Change = pred_result$main_effect_change,
            Main_Effect_Change_CI_Lower = pred_result$main_effect_change_ci[1],
            Main_Effect_Change_CI_Upper = pred_result$main_effect_change_ci[2],
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || 
              pred_result$improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- pred_result
          }
        }
      }
    }
    
    # Show best single predictor
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      cat("\n", paste(rep("-", 40), collapse=""), "\n")
      cat("Best single predictor:", result$best_single_addition$predictor, "\n")
      cat("AIC improvement:", round(result$best_single_addition$improvement, 2), "\n")
      
      # Show bootstrap coefficients
      cat("\nExpanded model coefficients (with 95% bootstrap CIs):\n")
      coefs <- result$best_single_addition$new_bootstrap$coefficients
      for (i in 1:nrow(coefs)) {
        if (abs(coefs$Estimate[i]) > 0.001) {  # Skip very small coefficients
          cat(sprintf("  %-30s: %6.3f [%6.3f, %6.3f]\n",
                      coefs$Parameter[i],
                      coefs$Estimate[i],
                      coefs$CI_Lower[i],
                      coefs$CI_Upper[i]))
        }
      }
    } else {
      cat("\nNo predictors improved the model by >", aic_threshold, "AIC units\n")
    }
    
    expansion_results[[analysis_name]] <- result
  }
  
  # Save results
  saveRDS(expansion_results, file.path(output_dir, "expansion_results_bootstrap.rds"))
  
  # Create comprehensive table with bootstrap CIs
  cat("\n\nCreating comprehensive results table...\n")
  
  comprehensive_table <- data.frame()
  
  for (analysis_name in names(expansion_results)) {
    result <- expansion_results[[analysis_name]]
    
    for (pred_name in names(result$predictor_results)) {
      pred_result <- result$predictor_results[[pred_name]]
      
      # Get all coefficients from expanded model
      expanded_coefs <- pred_result$new_bootstrap$coefficients
      
      # Create row for each coefficient
      for (i in 1:nrow(expanded_coefs)) {
        comprehensive_table <- rbind(comprehensive_table, data.frame(
          Analysis = analysis_name,
          Predictor_Added = pred_name,
          Variable_Name = pred_result$variable,
          N_Species = pred_result$n_species,
          Base_AIC_Subset = pred_result$base_subset_aic,
          Expanded_AIC = pred_result$new_aic,
          AIC_Improvement = pred_result$improvement,
          Parameter = expanded_coefs$Parameter[i],
          Estimate = expanded_coefs$Estimate[i],
          Bootstrap_Mean = expanded_coefs$Boot_Mean[i],
          Bootstrap_SD = expanded_coefs$Boot_SD[i],
          CI_Lower_95 = expanded_coefs$CI_Lower[i],
          CI_Upper_95 = expanded_coefs$CI_Upper[i],
          Significant = (expanded_coefs$CI_Lower[i] * expanded_coefs$CI_Upper[i]) > 0,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  write.csv(comprehensive_table, 
            file.path(output_dir, "stepwise_bootstrap_comprehensive_table.csv"),
            row.names = FALSE)
  
  # Create summary of main effects
  main_effects_summary <- expansion_results[[1]]$improvements
  for (i in 2:length(expansion_results)) {
    if (nrow(expansion_results[[i]]$improvements) > 0) {
      df <- expansion_results[[i]]$improvements
      df$Analysis <- names(expansion_results)[i]
      main_effects_summary$Analysis <- names(expansion_results)[1]
      main_effects_summary <- rbind(main_effects_summary, df)
    }
  }
  
  write.csv(main_effects_summary,
            file.path(output_dir, "main_effects_changes_with_ci.csv"),
            row.names = FALSE)
  
  # Create visualization with CIs
  create_bootstrap_plots(expansion_results, output_dir, aic_threshold)
  
  # Summary text
  summary_text <- "STEPWISE EXPANSION WITH BOOTSTRAP CONFIDENCE INTERVALS\n"
  summary_text <- paste0(summary_text, paste(rep("=", 60), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Generated: ", Sys.Date(), "\n")
  summary_text <- paste0(summary_text, "Bootstrap iterations: ", n_bootstrap, "\n")
  summary_text <- paste0(summary_text, "AIC improvement threshold: ", aic_threshold, "\n\n")
  
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
        if (!is.na(imp$Main_Effect_Change)) {
          summary_text <- paste0(summary_text, "     Main effect change: ", 
                                sprintf("%.3f [%.3f, %.3f]\n", 
                                        imp$Main_Effect_Change,
                                        imp$Main_Effect_Change_CI_Lower,
                                        imp$Main_Effect_Change_CI_Upper))
        }
      }
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "bootstrap_summary.txt"))
  cat("\n", summary_text)
  
  cat("\nResults saved to:", output_dir, "\n")
  
  return(expansion_results)
}

# Function to create plots with bootstrap CIs
create_bootstrap_plots <- function(expansion_results, output_dir, aic_threshold) {
  
  # 1. Main effects change plot with CIs
  plot_data <- data.frame()
  
  for (analysis_name in names(expansion_results)) {
    result <- expansion_results[[analysis_name]]
    if (nrow(result$improvements) > 0) {
      df <- result$improvements
      df$Analysis <- analysis_name
      df$Significant_Improvement <- df$AIC_Improvement > aic_threshold
      plot_data <- rbind(plot_data, df)
    }
  }
  
  if (nrow(plot_data) > 0 && any(!is.na(plot_data$Main_Effect_Change))) {
    # Filter to only predictors with main effect changes
    plot_data_effects <- plot_data[!is.na(plot_data$Main_Effect_Change), ]
    
    p1 <- ggplot(plot_data_effects, 
                 aes(x = reorder(Predictor, AIC_Improvement), 
                     y = Main_Effect_Change,
                     color = Analysis)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = Main_Effect_Change_CI_Lower,
                        ymax = Main_Effect_Change_CI_Upper),
                    position = position_dodge(width = 0.5),
                    width = 0.2) +
      coord_flip() +
      scale_color_brewer(palette = "Set1",
                        labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding",
                                  "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song")) +
      labs(title = "Change in Main Association When Adding Predictors",
           subtitle = "With 95% bootstrap confidence intervals",
           x = "Added Predictor",
           y = "Change in Main Effect Coefficient",
           color = "Direction") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "main_effect_changes_with_ci.png"),
           p1, width = 10, height = 8, dpi = 300)
  }
  
  # 2. AIC improvement plot
  if (nrow(plot_data) > 0) {
    p2 <- ggplot(plot_data, 
                 aes(x = reorder(Predictor, AIC_Improvement), 
                     y = AIC_Improvement,
                     fill = Analysis)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_text(aes(label = paste0("n=", N_Species)), 
                position = position_dodge(width = 0.9),
                hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set1",
                       labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding",
                                 "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song")) +
      labs(title = "Stepwise Model Improvements",
           subtitle = paste("Based on", unique(plot_data$N_Species), "bootstrap iterations per model"),
           x = "",
           y = "AIC Improvement",
           fill = "Direction") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "aic_improvements_bootstrap.png"),
           p2, width = 10, height = 6, dpi = 300)
  }
}

# Run the analysis
# For testing, use test_mode = TRUE
# results <- run_stepwise_expansion_bootstrap(test_mode = TRUE)

# For full analysis with 1000 bootstraps
# results <- run_stepwise_expansion_bootstrap(n_bootstrap = 1000)