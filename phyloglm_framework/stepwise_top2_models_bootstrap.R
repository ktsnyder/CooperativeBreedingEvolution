# Stepwise expansion for top 2 models with bootstrap confidence intervals
# This version includes proper uncertainty quantification

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

# Bootstrap function for a single model
run_model_with_bootstrap <- function(formula, data, tree, n_boot = 1000, 
                                   method = "logistic_MPLE") {
  
  # Fit original model
  original_fit <- tryCatch({
    phyloglm(
      formula = formula,
      data = data,
      phy = tree,
      method = method,
      btol = 50,
      log.alpha.bound = 4
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(original_fit)) {
    return(NULL)
  }
  
  # Extract original metrics
  original_aic <- -2 * original_fit$logLik + 2 * original_fit$d
  original_coef <- coef(original_fit)
  n_params <- length(original_coef)
  
  # Bootstrap
  boot_coefs <- matrix(NA, nrow = n_boot, ncol = n_params)
  colnames(boot_coefs) <- names(original_coef)
  
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  
  for (i in 1:n_boot) {
    boot_fit <- tryCatch({
      phyloglm(
        formula = formula,
        data = data,
        phy = tree,
        method = method,
        btol = 50,
        log.alpha.bound = 4,
        boot = 1
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(boot_fit)) {
      boot_coefs[i, ] <- coef(boot_fit)
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate bootstrap statistics
  boot_means <- colMeans(boot_coefs, na.rm = TRUE)
  boot_sds <- apply(boot_coefs, 2, sd, na.rm = TRUE)
  boot_lower <- apply(boot_coefs, 2, quantile, probs = 0.025, na.rm = TRUE)
  boot_upper <- apply(boot_coefs, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  # Create coefficient summary
  coef_summary <- data.frame(
    Parameter = names(original_coef),
    Estimate = original_coef,
    Boot_Mean = boot_means,
    Boot_SD = boot_sds,
    CI_Lower = boot_lower,
    CI_Upper = boot_upper,
    Significant = (boot_lower * boot_upper) > 0,
    stringsAsFactors = FALSE
  )
  
  return(list(
    fit = original_fit,
    aic = original_aic,
    logLik = original_fit$logLik,
    alpha = original_fit$alpha,
    coefficients = coef_summary,
    n_successful_boots = sum(!is.na(boot_coefs[,1])),
    bootstrap_matrix = boot_coefs
  ))
}

# Extract main effects with bootstrap CIs
extract_main_effects_bootstrap <- function(boot_result, analysis_direction) {
  if (is.null(boot_result)) return(NULL)
  
  coef_summary <- boot_result$coefficients
  
  # Determine main predictor based on analysis direction
  if (grepl("^FS_vs_CB", analysis_direction)) {
    main_pred <- "HighConfidence_Coop"
    outcome <- "Female Song"
  } else {
    main_pred <- "FemaleSong_Agg01"
    outcome <- "Cooperative Breeding"
  }
  
  # Find main effect row
  main_rows <- grep(paste0("^", main_pred, "$"), coef_summary$Parameter)
  
  if (length(main_rows) > 0) {
    main_row <- coef_summary[main_rows[1], ]
    
    # Calculate odds ratio with CIs
    or <- exp(main_row$Estimate)
    or_lower <- exp(main_row$CI_Lower)
    or_upper <- exp(main_row$CI_Upper)
    
    return(data.frame(
      Main_Predictor = main_pred,
      Outcome = outcome,
      Coefficient = main_row$Estimate,
      Boot_Mean = main_row$Boot_Mean,
      Boot_SD = main_row$Boot_SD,
      CI_Lower = main_row$CI_Lower,
      CI_Upper = main_row$CI_Upper,
      P_Value = 2 * pnorm(-abs(main_row$Estimate / main_row$Boot_SD)),
      Odds_Ratio = or,
      OR_CI_Lower = or_lower,
      OR_CI_Upper = or_upper,
      Significant = main_row$Significant,
      stringsAsFactors = FALSE
    ))
  }
  
  return(NULL)
}

# Main function
run_stepwise_top2_models_bootstrap <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2,
  n_bootstrap = 1000,
  test_mode = FALSE
) {
  
  if (test_mode) {
    n_bootstrap <- 10
    cat("*** TEST MODE: Only running", n_bootstrap, "bootstraps ***\n\n")
  }
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_top2_bootstrap_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Storage for all results
  all_expansion_results <- list()
  comprehensive_table <- data.frame()
  main_coefficients_tracking <- data.frame()
  
  # Process each analysis
  analyses <- c("FS_vs_CB_Terr3_Mass", "FS_vs_CB_TerrWS_Mass", 
                "CB_vs_FS_Terr3_Mass", "CB_vs_FS_TerrWS_Mass")
  
  for (analysis_name in analyses) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("Analyzing:", analysis_name, "\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
    
    # Get top 2 models
    base_result <- all_results[[analysis_name]]
    model_comparison <- base_result$comparison$comparison
    top_models <- head(model_comparison, 2)
    
    # Process each of the top 2 models
    for (model_rank in 1:nrow(top_models)) {
      model_name <- top_models$Model[model_rank]
      model_aic <- top_models$AIC[model_rank]
      
      cat("\nModel", model_rank, ":", model_name, "(AIC =", round(model_aic, 2), ")\n")
      cat(paste(rep("-", 50), collapse=""), "\n")
      
      # Get base model
      base_model <- base_result$models$models[[model_name]]
      base_data <- base_result$prepared_data$data
      base_tree <- base_result$prepared_data$tree
      base_formula <- formula(base_model)
      
      # Track main effect for base model
      cat("Fitting base model with bootstrap...\n")
      base_boot <- run_model_with_bootstrap(
        formula = base_formula,
        data = base_data,
        tree = base_tree,
        n_boot = n_bootstrap
      )
      
      base_main_effect <- extract_main_effects_bootstrap(base_boot, analysis_name)
      if (!is.null(base_main_effect)) {
        base_main_effect$Analysis <- analysis_name
        base_main_effect$Model_Rank <- model_rank
        base_main_effect$Model_Name <- model_name
        base_main_effect$Predictor_Added <- "Base"
        base_main_effect$N_Species <- nrow(base_data)
        base_main_effect$AIC <- model_aic
        base_main_effect$AIC_Improvement <- 0
        main_coefficients_tracking <- rbind(main_coefficients_tracking, base_main_effect)
      }
      
      # Test each predictor
      cat("\nTesting additional predictors:\n")
      
      predictors_to_test <- list(
        list(name = "Territory (1-3)", var = "Territory_num", 
             values = as.numeric(full_data$Territory)),
        list(name = "Migration (1-3)", var = "Migration_num",
             values = as.numeric(full_data$Migration_AVONET)),
        list(name = "Absolute Latitude", var = "abs_Latitude",
             values = abs(full_data$Centroid.Latitude_AVONET)),
        list(name = "Wing Dimorphism", var = "PercentAbsLogWingDimorphism",
             values = full_data$PercentAbsLogWingDimorphism),
        list(name = "Familial Living", var = "Griesser2017FamilialLiving",
             values = full_data$Griesser2017FamilialLiving),
        list(name = "Plumage Dimorphism", var = "logMaleFemalePlumageDiffAbs",
             values = log(full_data$MaleFemalePlumageDiffAbs + 1)),
        list(name = "Geographic Region", var = "GeographicRegion_Jetz",
             values = full_data$GeographicRegion_Jetz)
      )
      
      # Remove predictors already in the model
      base_vars <- all.vars(base_formula)
      predictors_to_test <- predictors_to_test[!sapply(predictors_to_test, 
                                                       function(p) p$var %in% base_vars)]
      
      model_improvements <- data.frame()
      
      for (pred in predictors_to_test) {
        cat("\n  Testing", pred$name, "... ")
        
        # Add predictor to data
        test_data <- base_data
        test_data[[pred$var]] <- pred$values[match(test_data$species, full_data$species)]
        
        # Remove NAs
        test_data <- test_data[!is.na(test_data[[pred$var]]), ]
        
        if (nrow(test_data) < 100) {
          cat("insufficient data (n =", nrow(test_data), ")\n")
          next
        }
        
        # Match tree
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        cat("n =", nrow(test_data), "\n")
        
        # Refit base model on subset
        cat("    Refitting base model on subset with bootstrap...\n")
        base_subset_boot <- run_model_with_bootstrap(
          formula = base_formula,
          data = test_data,
          tree = test_tree,
          n_boot = n_bootstrap
        )
        
        if (is.null(base_subset_boot)) {
          cat("    Base model failed on subset\n")
          next
        }
        
        # Fit expanded model
        new_formula <- update(base_formula, paste("~ . +", pred$var))
        cat("    Fitting expanded model with bootstrap...\n")
        expanded_boot <- run_model_with_bootstrap(
          formula = new_formula,
          data = test_data,
          tree = test_tree,
          n_boot = n_bootstrap
        )
        
        if (is.null(expanded_boot)) {
          cat("    Expanded model failed\n")
          next
        }
        
        # Calculate improvement
        improvement <- base_subset_boot$aic - expanded_boot$aic
        cat("    AIC improvement:", round(improvement, 2), "\n")
        
        # Extract main effects
        expanded_main <- extract_main_effects_bootstrap(expanded_boot, analysis_name)
        
        # Track coefficient changes
        if (!is.null(expanded_main)) {
          expanded_main$Analysis <- analysis_name
          expanded_main$Model_Rank <- model_rank
          expanded_main$Model_Name <- model_name
          expanded_main$Predictor_Added <- pred$name
          expanded_main$N_Species <- nrow(test_data)
          expanded_main$AIC <- expanded_boot$aic
          expanded_main$AIC_Improvement <- improvement
          main_coefficients_tracking <- rbind(main_coefficients_tracking, expanded_main)
        }
        
        # Add to comprehensive table
        comprehensive_row <- data.frame(
          Analysis = analysis_name,
          Model_Rank = model_rank,
          Base_Model_Name = model_name,
          Base_Model_Formula = paste(deparse(base_formula), collapse = " "),
          Original_Base_AIC = model_aic,
          Predictor_Added = pred$name,
          Predictor_Variable_Name = pred$var,
          N_Species_Original = nrow(base_data),
          N_Species_Subset = nrow(test_data),
          Base_AIC_Subset = base_subset_boot$aic,
          Expanded_Model_AIC = expanded_boot$aic,
          AIC_Improvement = improvement,
          Expanded_Model_Formula = paste(deparse(new_formula), collapse = " "),
          Phylogenetic_Signal_Base = base_subset_boot$fit$alpha,
          Phylogenetic_Signal_Expanded = expanded_boot$fit$alpha,
          LogLik_Base = base_subset_boot$logLik,
          LogLik_Expanded = expanded_boot$logLik,
          N_Bootstrap = n_bootstrap,
          N_Successful_Base = base_subset_boot$n_successful_boots,
          N_Successful_Expanded = expanded_boot$n_successful_boots,
          stringsAsFactors = FALSE
        )
        
        comprehensive_table <- rbind(comprehensive_table, comprehensive_row)
        
        # Store improvement info
        model_improvements <- rbind(model_improvements, data.frame(
          Predictor = pred$name,
          AIC_Improvement = improvement,
          Main_Effect_Change = ifelse(!is.null(expanded_main) && !is.null(base_main_effect),
                                     expanded_main$Coefficient - base_main_effect$Coefficient,
                                     NA),
          stringsAsFactors = FALSE
        ))
      }
      
      # Store results for this model
      model_key <- paste0(analysis_name, "_Model", model_rank)
      all_expansion_results[[model_key]] <- list(
        improvements = model_improvements,
        base_formula = base_formula,
        n_species = nrow(base_data)
      )
    }
  }
  
  # Save all results
  write.csv(comprehensive_table, 
            file.path(output_dir, "top2_models_comprehensive_table_bootstrap.csv"),
            row.names = FALSE)
  
  write.csv(main_coefficients_tracking,
            file.path(output_dir, "main_coefficients_comparison_bootstrap.csv"),
            row.names = FALSE)
  
  saveRDS(all_expansion_results,
          file.path(output_dir, "top2_models_expansion_results_bootstrap.rds"))
  
  # Create visualizations
  create_top2_bootstrap_plots(main_coefficients_tracking, comprehensive_table, 
                             output_dir, aic_threshold)
  
  # Summary
  cat("\n\nResults saved to:", output_dir, "\n")
  
  return(list(
    comprehensive_table = comprehensive_table,
    main_coefficients = main_coefficients_tracking,
    expansion_results = all_expansion_results
  ))
}

# Create plots with bootstrap CIs
create_top2_bootstrap_plots <- function(main_coefs, comp_table, output_dir, aic_threshold) {
  
  # 1. Main coefficients comparison with CIs
  plot_data <- main_coefs %>%
    filter(Predictor_Added != "Base") %>%
    mutate(
      Analysis_Clean = case_when(
        grepl("CB_vs_FS_TerrWS", Analysis) ~ "CB vs FS Terr W/S",
        grepl("CB_vs_FS_Terr3", Analysis) ~ "CB vs FS Terr12vs3", 
        grepl("FS_vs_CB_TerrWS", Analysis) ~ "FS vs CB Terr W/S",
        grepl("FS_vs_CB_Terr3", Analysis) ~ "FS vs CB Terr12vs3",
        TRUE ~ Analysis
      ),
      Panel_Title = paste0(Model_Name, " (Model ", Model_Rank, ")"),
      Significant_Improvement = AIC_Improvement > aic_threshold
    ) %>%
    arrange(Analysis, Model_Rank, desc(AIC_Improvement))
  
  # Create faceted plot
  p1 <- ggplot(plot_data, aes(x = reorder(Predictor_Added, AIC_Improvement), 
                               y = Coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(aes(color = Significant_Improvement), size = 3) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper,
                      color = Significant_Improvement),
                  width = 0.2) +
    facet_grid(Analysis_Clean ~ paste0("Model Rank ", Model_Rank),
               scales = "free_x") +
    coord_flip() +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray50"),
                       labels = c("TRUE" = "Significant (>2 AIC)", 
                                 "FALSE" = "Not significant"),
                       name = "AIC Improvement") +
    labs(title = "Main Association Coefficients with 95% Bootstrap CIs",
         subtitle = paste0("Based on ", unique(plot_data$N_Bootstrap)[1], " bootstrap iterations"),
         x = "Added Predictor",
         y = "Main Association Coefficient") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "main_coefficients_comparison_bootstrap.png"),
         p1, width = 14, height = 10, dpi = 300)
  
  # 2. AIC improvements plot
  imp_data <- comp_table %>%
    mutate(
      Analysis_Clean = case_when(
        grepl("CB_vs_FS_TerrWS", Analysis) ~ "CB vs FS Terr W/S",
        grepl("CB_vs_FS_Terr3", Analysis) ~ "CB vs FS Terr12vs3",
        grepl("FS_vs_CB_TerrWS", Analysis) ~ "FS vs CB Terr W/S", 
        grepl("FS_vs_CB_Terr3", Analysis) ~ "FS vs CB Terr12vs3",
        TRUE ~ Analysis
      ),
      Significant = AIC_Improvement > aic_threshold
    )
  
  p2 <- ggplot(imp_data, aes(x = reorder(Predictor_Added, AIC_Improvement),
                              y = AIC_Improvement,
                              fill = factor(Model_Rank))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
    geom_text(aes(label = paste0("n=", N_Species_Subset)),
              position = position_dodge(width = 0.9),
              hjust = -0.1, size = 3) +
    facet_wrap(~ Analysis_Clean, scales = "free", ncol = 2) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2", name = "Model Rank") +
    labs(title = "AIC Improvements for Top 2 Models per Analysis",
         subtitle = "Red line indicates threshold of 2 AIC units",
         x = "Added Predictor",
         y = "AIC Improvement") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "aic_improvements_by_model_bootstrap.png"),
         p2, width = 12, height = 10, dpi = 300)
}

# Run the analysis
# For testing:
# results <- run_stepwise_top2_models_bootstrap(test_mode = TRUE)

# For full analysis:
# results <- run_stepwise_top2_models_bootstrap(n_bootstrap = 1000)