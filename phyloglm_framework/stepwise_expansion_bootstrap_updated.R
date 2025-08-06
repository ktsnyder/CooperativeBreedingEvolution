# Stepwise expansion with bootstrap confidence intervals - UPDATED VERSION
# This version implements all proposed improvements:
# - Removes invalid change CIs
# - Includes all 4 analyses
# - Saves bootstrap matrices as RDS files
# - Improves column names
# - Tracks convergence
# - Handles categorical variables properly
# - Uses message() instead of cat() for reliable console output

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

# Utility function to load bootstrap matrices when needed
load_bootstrap_matrix <- function(matrix_file) {
  if (!is.null(matrix_file) && file.exists(matrix_file)) {
    return(readRDS(matrix_file))
  }
  return(NULL)
}

# Function to run bootstrap for a single model
run_bootstrap_model <- function(formula, data, tree, n_boot = 1000, 
                               method = "logistic_MPLE", save_prefix = NULL,
                               save_matrices = TRUE, matrix_dir = NULL) {
  
  message("      Running ", n_boot, " bootstrap iterations...")
  
  # Fit original model
  original_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4
  )
  
  # Get bootstrap results directly from phyloglm
  # boot parameter in phyloglm specifies number of bootstrap replicates
  message("      Fitting model with bootstrap=", n_boot, "...")
  boot_fit <- phyloglm(
    formula = formula,
    data = data,
    phy = tree,
    method = method,
    btol = 50,
    log.alpha.bound = 4,
    boot = n_boot  # This will run n_boot bootstrap iterations internally
  )
  
  # Extract bootstrap results
  if (!is.null(boot_fit$bootstrap)) {
    boot_matrix <- boot_fit$bootstrap
    
    # Save matrix to RDS if requested
    matrix_file <- NULL
    if (save_matrices && !is.null(save_prefix) && !is.null(matrix_dir)) {
      matrix_file <- file.path(matrix_dir, paste0(save_prefix, "_boot_matrix.rds"))
      saveRDS(boot_matrix, matrix_file, compress = TRUE)
      message("      Saved bootstrap matrix to: ", basename(matrix_file))
    }
    
    # Calculate statistics
    original_coef <- coef(original_fit)
    
    # Make sure we're only working with regression coefficients
    # (bootstrap matrix might include alpha or other parameters)
    coef_names <- names(original_coef)
    n_coef <- length(coef_names)
    
    # Extract only the columns that correspond to regression coefficients
    if (ncol(boot_matrix) > n_coef) {
      boot_matrix_coef <- boot_matrix[, 1:n_coef, drop = FALSE]
    } else {
      boot_matrix_coef <- boot_matrix
    }
    
    boot_means <- colMeans(boot_matrix_coef, na.rm = TRUE)
    boot_sds <- apply(boot_matrix_coef, 2, sd, na.rm = TRUE)
    boot_lower <- apply(boot_matrix_coef, 2, quantile, probs = 0.025, na.rm = TRUE)
    boot_upper <- apply(boot_matrix_coef, 2, quantile, probs = 0.975, na.rm = TRUE)
    
    # Create summary
    coef_summary <- data.frame(
      Parameter = coef_names,
      Estimate = original_coef,
      Boot_Mean = boot_means,
      Boot_SD = boot_sds,
      CI_Lower = boot_lower,
      CI_Upper = boot_upper,
      stringsAsFactors = FALSE
    )
    
    n_successful <- sum(complete.cases(boot_matrix_coef))
    n_converged <- sum(!is.na(boot_matrix_coef[,1]))  # Check first coefficient column
    convergence_rate <- n_converged / n_boot
    
  } else {
    # Fallback if bootstrap didn't work
    message("      Warning: Bootstrap results not available, using standard errors")
    coef_summary <- summary(original_fit)$coefficients
    coef_summary <- data.frame(
      Parameter = rownames(coef_summary),
      Estimate = coef_summary[, "Estimate"],
      Boot_Mean = coef_summary[, "Estimate"],
      Boot_SD = coef_summary[, "StdErr"],
      CI_Lower = coef_summary[, "Estimate"] - 1.96 * coef_summary[, "StdErr"],
      CI_Upper = coef_summary[, "Estimate"] + 1.96 * coef_summary[, "StdErr"],
      stringsAsFactors = FALSE
    )
    matrix_file <- NULL
    n_successful <- 1
    n_converged <- 1
    convergence_rate <- 1
  }
  
  # Save interim results if prefix provided
  if (!is.null(save_prefix) && !save_matrices) {
    interim_file <- file.path(matrix_dir, paste0(save_prefix, "_bootstrap_result.rds"))
    saveRDS(list(
      fit = original_fit,
      coefficients = coef_summary,
      bootstrap_matrix_file = matrix_file
    ), interim_file)
  }
  
  return(list(
    fit = original_fit,
    coefficients = coef_summary,
    bootstrap_matrix_file = matrix_file,  # Store path instead of matrix
    n_successful_boots = n_successful,
    n_converged = n_converged,
    convergence_rate = convergence_rate
  ))
}

run_stepwise_expansion_bootstrap <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2,
  n_bootstrap = 1000,
  ncores = 1,
  test_mode = FALSE,
  save_interim = TRUE,
  save_boot_matrices = TRUE,      # Whether to save bootstrap matrices
  boot_matrix_dir = NULL          # Subdirectory for matrix files
) {
  
  message("Starting stepwise expansion with bootstrap...")
  
  if (test_mode) {
    n_bootstrap <- 100  # Use 100 for test mode to get some variation
    message("*** TEST MODE: Only running ", n_bootstrap, " bootstraps ***\n")
  }
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_boot", n_bootstrap, "_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create subdirectory for interim saves
  if (save_interim) {
    interim_dir <- file.path(output_dir, "interim_results")
    dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create subdirectory for bootstrap matrices
  if (save_boot_matrices) {
    if (is.null(boot_matrix_dir)) {
      boot_matrix_dir <- file.path(output_dir, "bootstrap_matrices")
    }
    dir.create(boot_matrix_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Load data
  message("Loading data and results...")
  
  message("  Loading results from: ", results_path)
  all_results <- readRDS(results_path)
  
  message("  Loading data from: ", data_path)
  full_data <- read.csv(data_path)
  
  message("  Loading tree from: ", tree_path)
  tree <- read.nexus(tree_path)
  
  message("  Data loaded successfully!")
  
  # Results storage
  expansion_results <- list()
  
  # Analyze each direction - NOW INCLUDING ALL 4 ANALYSES
  for (analysis_name in c("CB_vs_FS_TerrWS_Mass", "FS_vs_CB_TerrWS_Mass", 
                          "CB_vs_FS_Terr3_Mass", "FS_vs_CB_Terr3_Mass")) {
    
    message("\n", paste(rep("=", 60), collapse=""))
    message("Analyzing: ", analysis_name)
    message(paste(rep("=", 60), collapse=""), "\n")
    
    # Get base model
    base_result <- all_results[[analysis_name]]
    if (is.null(base_result)) {
      message("Analysis ", analysis_name, " not found. Skipping.")
      next
    }
    
    best_model_name <- base_result$comparison$comparison$Model[1]
    best_model <- base_result$models$models[[best_model_name]]
    original_base_aic <- base_result$comparison$comparison$AIC[1]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    base_formula <- formula(best_model)
    
    message("Base model: ", best_model_name)
    message("Original AIC (", nrow(base_data), " species): ", round(original_base_aic, 2))
    message("Base formula: ", deparse(base_formula), "\n")
    
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
    message("Testing predictors with bootstrap CIs:")
    message(paste(rep("-", 40), collapse=""))
    
    # Helper function to test a predictor with bootstrap
    test_predictor_bootstrap <- function(predictor_name, predictor_values, var_name) {
      message("\n", predictor_name, ":")
      
      # Add predictor to base data
      test_data <- base_data
      test_data[[var_name]] <- predictor_values[match(test_data$species, full_data$species)]
      
      # Convert categorical variables to factors
      if (var_name == "GeographicRegion_Jetz") {
        test_data[[var_name]] <- as.factor(test_data[[var_name]])
      }
      
      # Remove NAs
      test_data <- test_data[!is.na(test_data[[var_name]]), ]
      
      if (nrow(test_data) < 100) {
        message("   Not enough species with data (", nrow(test_data), ")")
        return(NULL)
      }
      
      # Match tree
      test_tree <- keep.tip(base_tree, test_data$species)
      test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
      
      message("   N species with ", predictor_name, " data: ", nrow(test_data))
      
      # Create save prefix for interim results
      save_prefix <- NULL
      if (save_interim) {
        save_prefix <- paste0(analysis_name, "_", 
                             gsub("[^A-Za-z0-9]", "", predictor_name))
      }
      
      # Refit base model to this subset WITH BOOTSTRAP
      message("   Fitting base model on subset with bootstrap...")
      base_subset_boot <- run_bootstrap_model(
        formula = base_formula,
        data = test_data,
        tree = test_tree,
        n_boot = n_bootstrap,
        save_prefix = if(save_interim) paste0(save_prefix, "_base") else NULL,
        save_matrices = save_boot_matrices,
        matrix_dir = boot_matrix_dir
      )
      base_subset_fit <- base_subset_boot$fit
      base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
      message("   Base model AIC on subset: ", round(base_subset_aic, 2))
      message("   Bootstrap convergence rate: ", sprintf("%.1f%%", base_subset_boot$convergence_rate * 100))
      
      # Now fit model with predictor WITH BOOTSTRAP
      new_formula <- update(base_formula, paste("~ . +", var_name))
      message("   Fitting expanded model with bootstrap...")
      new_boot <- run_bootstrap_model(
        formula = new_formula,
        data = test_data,
        tree = test_tree,
        n_boot = n_bootstrap,
        save_prefix = if(save_interim) paste0(save_prefix, "_expanded") else NULL,
        save_matrices = save_boot_matrices,
        matrix_dir = boot_matrix_dir
      )
      new_fit <- new_boot$fit
      new_aic <- -2 * new_fit$logLik + 2 * new_fit$d
      
      # Calculate improvement
      improvement <- base_subset_aic - new_aic
      
      message("   Model with ", predictor_name, " AIC: ", round(new_aic, 2))
      message("   AIC improvement: ", round(improvement, 2))
      message("   Bootstrap convergence rate: ", sprintf("%.1f%%", new_boot$convergence_rate * 100))
      
      # Extract main effect changes - UPDATED SECTION
      response_var <- all.vars(base_formula)[1]
      if (response_var == "FemaleSong_Agg01") {
        main_pred <- "HighConfidence_Coop"
      } else {
        main_pred <- "FemaleSong_Agg01"
      }
      
      # Get main effect from both models
      base_main <- base_subset_boot$coefficients[base_subset_boot$coefficients$Parameter == main_pred, ]
      new_main <- new_boot$coefficients[new_boot$coefficients$Parameter == main_pred, ]
      
      # Initialize with NA values
      main_effect_change <- NA
      base_main_coef <- NA
      base_main_ci <- c(NA, NA)
      expanded_main_coef <- NA
      expanded_main_ci <- c(NA, NA)
      
      if (nrow(base_main) > 0 && nrow(new_main) > 0) {
        # Store coefficients separately
        base_main_coef <- base_main$Estimate
        base_main_ci <- c(base_main$CI_Lower, base_main$CI_Upper)
        expanded_main_coef <- new_main$Estimate
        expanded_main_ci <- c(new_main$CI_Lower, new_main$CI_Upper)
        
        # Calculate point estimate of change (for plotting/sorting purposes only)
        main_effect_change <- expanded_main_coef - base_main_coef
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
        base_main_coef = base_main_coef,
        base_main_ci = base_main_ci,
        expanded_main_coef = expanded_main_coef,
        expanded_main_ci = expanded_main_ci,
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
          
          # Add to summary table - UPDATED COLUMN NAMES
          result$improvements <- rbind(result$improvements, data.frame(
            Predictor_Name = pred_result$predictor,
            Predictor_Variable = pred_result$variable,
            N_Species = pred_result$n_species,
            Base_Model_AIC_Subset = pred_result$base_subset_aic,
            Expanded_Model_AIC = pred_result$new_aic,
            AIC_Improvement = pred_result$improvement,
            Improved_Model = pred_result$improvement > aic_threshold,
            Main_Effect_Base_Coef = pred_result$base_main_coef,
            Main_Effect_Base_CI_Lower = pred_result$base_main_ci[1],
            Main_Effect_Base_CI_Upper = pred_result$base_main_ci[2],
            Main_Effect_Expanded_Coef = pred_result$expanded_main_coef,
            Main_Effect_Expanded_CI_Lower = pred_result$expanded_main_ci[1],
            Main_Effect_Expanded_CI_Upper = pred_result$expanded_main_ci[2],
            Main_Effect_Change_Estimate = pred_result$main_effect_change,
            stringsAsFactors = FALSE
          ))
          
          if (is.null(result$best_single_addition) || 
              pred_result$improvement > result$best_single_addition$improvement) {
            result$best_single_addition <- pred_result
          }
          
          # Save interim summary after each predictor
          if (save_interim) {
            saveRDS(result, file.path(interim_dir, paste0(analysis_name, "_current_results.rds")))
          }
        }
      }
    }
    
    # Show best single predictor
    if (!is.null(result$best_single_addition) && result$best_single_addition$improvement > aic_threshold) {
      message("\n", paste(rep("-", 40), collapse=""))
      message("Best single predictor: ", result$best_single_addition$predictor)
      message("AIC improvement: ", round(result$best_single_addition$improvement, 2))
      
      # Show bootstrap coefficients
      message("\nExpanded model coefficients (with 95% bootstrap CIs):")
      coefs <- result$best_single_addition$new_bootstrap$coefficients
      for (i in 1:nrow(coefs)) {
        if (abs(coefs$Estimate[i]) > 0.001) {  # Skip very small coefficients
          message(sprintf("  %-30s: %6.3f [%6.3f, %6.3f] (SD=%6.3f)",
                      coefs$Parameter[i],
                      coefs$Estimate[i],
                      coefs$CI_Lower[i],
                      coefs$CI_Upper[i],
                      coefs$Boot_SD[i]))
        }
      }
    } else {
      message("\nNo predictors improved the model by >", aic_threshold, " AIC units")
    }
    
    expansion_results[[analysis_name]] <- result
    
    # Save after each analysis completes
    if (save_interim) {
      saveRDS(expansion_results, file.path(output_dir, paste0("expansion_results_boot", n_bootstrap, "_partial.rds")))
    }
  }
  
  # Save final results
  saveRDS(expansion_results, file.path(output_dir, paste0("expansion_results_boot", n_bootstrap, ".rds")))
  
  # Create comprehensive table with bootstrap CIs - UPDATED COLUMN NAMES
  message("\n\nCreating comprehensive results table...")
  
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
          Analysis_Name = analysis_name,
          Analysis_Direction = ifelse(grepl("^FS_vs_CB", analysis_name), "FS→CB", "CB→FS"),
          Territory_Model = ifelse(grepl("Terr3", analysis_name), "Terr_1/3", "Terr_W/S"),
          Predictor_Added = pred_name,
          Predictor_Variable = pred_result$variable,
          N_Species_In_Model = pred_result$n_species,
          Base_Model_AIC = pred_result$base_subset_aic,
          Expanded_Model_AIC = pred_result$new_aic,
          AIC_Improvement = pred_result$improvement,
          Model_Improved = pred_result$improvement > aic_threshold,
          Coefficient_Name = expanded_coefs$Parameter[i],
          Coefficient_Estimate = expanded_coefs$Estimate[i],
          Bootstrap_Mean_Estimate = expanded_coefs$Boot_Mean[i],
          Bootstrap_SD = expanded_coefs$Boot_SD[i],
          CI_Lower_2.5 = expanded_coefs$CI_Lower[i],
          CI_Upper_97.5 = expanded_coefs$CI_Upper[i],
          Is_Significant = (expanded_coefs$CI_Lower[i] * expanded_coefs$CI_Upper[i]) > 0,
          Base_Bootstrap_Matrix_File = pred_result$base_bootstrap$bootstrap_matrix_file,
          Expanded_Bootstrap_Matrix_File = pred_result$new_bootstrap$bootstrap_matrix_file,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  write.csv(comprehensive_table, 
            file.path(output_dir, paste0("stepwise_boot", n_bootstrap, "_comprehensive_table.csv")),
            row.names = FALSE)
  
  # Create summary of main effects
  main_effects_summary <- expansion_results[[1]]$improvements
  main_effects_summary$Analysis = names(expansion_results)[1]
  for (i in 2:length(expansion_results)) {
    if (nrow(expansion_results[[i]]$improvements) > 0) {
      df <- expansion_results[[i]]$improvements
      df$Analysis <- names(expansion_results)[i]
      main_effects_summary <- rbind(main_effects_summary, df)
    }
  }
  
  write.csv(main_effects_summary,
            file.path(output_dir, paste0("main_effects_changes_boot", n_bootstrap, "_with_ci.csv")),
            row.names = FALSE)
  
  # Create visualization with CIs
  create_bootstrap_plots(expansion_results, output_dir, aic_threshold, n_bootstrap)
  
  # Summary text - UPDATED TO REFLECT SEPARATE CIs
  summary_text <- paste0("STEPWISE EXPANSION WITH BOOTSTRAP CONFIDENCE INTERVALS\n")
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
        summary_text <- paste0(summary_text, "  ", imp$Predictor_Name, 
                              ": AIC improvement = ", round(imp$AIC_Improvement, 2),
                              " (", imp$N_Species, " species)\n")
        if (!is.na(imp$Main_Effect_Base_Coef)) {
          summary_text <- paste0(summary_text, 
            sprintf("     Main effect in base model: %.3f [%.3f, %.3f]\n", 
                    imp$Main_Effect_Base_Coef,
                    imp$Main_Effect_Base_CI_Lower,
                    imp$Main_Effect_Base_CI_Upper))
          summary_text <- paste0(summary_text,
            sprintf("     Main effect in expanded model: %.3f [%.3f, %.3f]\n", 
                    imp$Main_Effect_Expanded_Coef,
                    imp$Main_Effect_Expanded_CI_Lower,
                    imp$Main_Effect_Expanded_CI_Upper))
          summary_text <- paste0(summary_text,
            sprintf("     Change in main effect: %.3f (no CI available)\n", 
                    imp$Main_Effect_Change_Estimate))
        }
      }
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  writeLines(summary_text, file.path(output_dir, paste0("bootstrap_summary_boot", n_bootstrap, ".txt")))
  message("\n", summary_text)
  
  message("\nResults saved to: ", output_dir)
  
  return(expansion_results)
}

# Function to create plots with bootstrap CIs - UPDATED
create_bootstrap_plots <- function(expansion_results, output_dir, aic_threshold, n_bootstrap) {
  
  # 1. Main effects plot showing base and expanded coefficients
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
  
  if (nrow(plot_data) > 0 && any(!is.na(plot_data$Main_Effect_Change_Estimate))) {
    # Filter to only predictors with main effect changes
    plot_data_effects <- plot_data[!is.na(plot_data$Main_Effect_Change_Estimate), ]
    
    # Reshape data for plotting both base and expanded coefficients
    plot_data_long <- rbind(
      data.frame(
        Analysis = plot_data_effects$Analysis,
        Predictor = plot_data_effects$Predictor_Name,
        Model = "Base",
        Coefficient = plot_data_effects$Main_Effect_Base_Coef,
        CI_Lower = plot_data_effects$Main_Effect_Base_CI_Lower,
        CI_Upper = plot_data_effects$Main_Effect_Base_CI_Upper,
        AIC_Improvement = plot_data_effects$AIC_Improvement
      ),
      data.frame(
        Analysis = plot_data_effects$Analysis,
        Predictor = plot_data_effects$Predictor_Name,
        Model = "Expanded",
        Coefficient = plot_data_effects$Main_Effect_Expanded_Coef,
        CI_Lower = plot_data_effects$Main_Effect_Expanded_CI_Lower,
        CI_Upper = plot_data_effects$Main_Effect_Expanded_CI_Upper,
        AIC_Improvement = plot_data_effects$AIC_Improvement
      )
    )
    
    p1 <- ggplot(plot_data_long, 
                 aes(x = reorder(Predictor, AIC_Improvement), 
                     y = Coefficient,
                     color = Analysis,
                     shape = Model)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                    position = position_dodge(width = 0.5),
                    width = 0.2) +
      coord_flip() +
      scale_color_brewer(palette = "Set1",
                        labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding (Terr W/S)",
                                  "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song (Terr W/S)",
                                  "FS_vs_CB_Terr3_Mass" = "Female Song → Coop. Breeding (Terr 1/3)",
                                  "CB_vs_FS_Terr3_Mass" = "Coop. Breeding → Female Song (Terr 1/3)")) +
      scale_shape_manual(values = c("Base" = 1, "Expanded" = 16)) +
      labs(title = "Main Association Coefficients: Base vs Expanded Models",
           subtitle = paste0("Showing base (○) and expanded (●) model coefficients with 95% bootstrap CIs (n=", 
                            n_bootstrap, ")"),
           x = "Added Predictor",
           y = "Main Effect Coefficient",
           color = "Direction",
           shape = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, paste0("main_effect_changes_boot", n_bootstrap, "_with_ci.png")),
           p1, width = 12, height = 8, dpi = 300)
  }
  
  # 2. AIC improvement plot
  if (nrow(plot_data) > 0) {
    p2 <- ggplot(plot_data, 
                 aes(x = reorder(Predictor_Name, AIC_Improvement), 
                     y = AIC_Improvement,
                     fill = Analysis)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_text(aes(label = paste0("n=", N_Species)), 
                position = position_dodge(width = 0.9),
                hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set1",
                       labels = c("FS_vs_CB_TerrWS_Mass" = "Female Song → Coop. Breeding (Terr W/S)",
                                 "CB_vs_FS_TerrWS_Mass" = "Coop. Breeding → Female Song (Terr W/S)",
                                 "FS_vs_CB_Terr3_Mass" = "Female Song → Coop. Breeding (Terr 1/3)",
                                 "CB_vs_FS_Terr3_Mass" = "Coop. Breeding → Female Song (Terr 1/3)")) +
      labs(title = "Stepwise Model Improvements",
           subtitle = paste0("Based on ", n_bootstrap, " bootstrap iterations per model"),
           x = "",
           y = "AIC Improvement",
           fill = "Direction") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, paste0("aic_improvements_boot", n_bootstrap, ".png")),
           p2, width = 10, height = 6, dpi = 300)
  }
}

# Run the analysis
# For testing, use test_mode = TRUE
# results <- run_stepwise_expansion_bootstrap(test_mode = TRUE)

# For full analysis with 1000 bootstraps
# results <- run_stepwise_expansion_bootstrap(n_bootstrap = 1000)