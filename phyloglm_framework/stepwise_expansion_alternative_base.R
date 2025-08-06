# Stepwise expansion using alternative base model for FS_vs_CB_Terr3_Mass
# This addresses the case where two models are within 0.3 AIC units

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)

source("phyloglm_framework/batch_runner_helpers.R")

run_stepwise_with_alternative_base <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_alternative_base_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Focus on FS_vs_CB_Terr3_Mass with alternative base
  analysis_name <- "FS_vs_CB_Terr3_Mass"
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("Analyzing:", analysis_name, "with ALTERNATIVE BASE MODEL\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  # Get the second-best model as base
  base_result <- all_results[[analysis_name]]
  model_comparison <- base_result$comparison$comparison
  
  cat("Model comparison (top 5):\n")
  print(model_comparison[1:5,])
  cat("\n")
  
  # Use MainPred_Terr_Mass as base
  best_model_name <- "MainPred_Terr_Mass"
  best_model <- base_result$models$models[[best_model_name]]
  original_base_aic <- model_comparison$AIC[model_comparison$Model == best_model_name]
  base_data <- base_result$prepared_data$data
  base_tree <- base_result$prepared_data$tree
  base_formula <- formula(best_model)
  base_formula_str <- paste(deparse(base_formula), collapse = " ")
  
  cat("Using alternative base model:", best_model_name, "\n")
  cat("Original AIC (", nrow(base_data), " species):", round(original_base_aic, 2), "\n")
  cat("Base formula:", base_formula_str, "\n\n")
  
  # Original sample size
  original_n <- nrow(base_data)
  
  # Results storage
  comprehensive_table <- data.frame()
  
  # Helper function to test a predictor
  test_predictor <- function(predictor_name, predictor_values, var_name, var_type = "continuous") {
    cat("\n", predictor_name, ":\n", sep="")
    
    # Add predictor to base data
    test_data <- base_data
    test_data[[var_name]] <- predictor_values[match(test_data$species, full_data$species)]
    
    # Count missing data
    n_missing_predictor <- sum(is.na(test_data[[var_name]]))
    
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
    
    # Initialize row data
    row_data <- data.frame(
      Analysis = paste0(analysis_name, "_AltBase"),
      Predictor = predictor_name,
      Base_Model_Name = best_model_name,
      Base_Formula = base_formula_str,
      Base_Original_AIC = original_base_aic,
      Base_Original_N = original_n,
      N_Species = nrow(test_data),
      N_Missing_Predictor = n_missing_predictor,
      Percent_Retained = round(100 * nrow(test_data) / original_n, 1),
      stringsAsFactors = FALSE
    )
    
    # CRITICAL: Refit base model to this subset
    cat("   Refitting base model to same species subset...\n")
    base_subset_fit <- tryCatch({
      phyloglm(
        formula = base_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
    }, error = function(e) {
      cat("   Error fitting base model:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(base_subset_fit)) {
      row_data$Note <- "Base model failed on subset"
      return(row_data)
    }
    
    # Extract base model info
    base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
    row_data$Base_Subset_AIC <- base_subset_aic
    row_data$Base_LogLik <- base_subset_fit$logLik
    row_data$Phylogenetic_Signal_Base <- base_subset_fit$alpha
    
    cat("   Base model AIC on subset:", round(base_subset_aic, 2), "\n")
    
    # Now fit model with predictor
    new_formula <- update(base_formula, paste("~ . +", var_name))
    new_formula_str <- paste(deparse(new_formula), collapse = " ")
    row_data$Predictor_Variable_Name <- var_name
    row_data$Predictor_Type <- var_type
    row_data$Expanded_Formula <- new_formula_str
    
    expanded_fit <- tryCatch({
      phyloglm(
        formula = new_formula,
        data = test_data,
        phy = test_tree,
        method = "logistic_MPLE",
        btol = 50,
        log.alpha.bound = 4
      )
    }, error = function(e) {
      cat("   Error fitting expanded model:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(expanded_fit)) {
      row_data$Note <- "Expanded model failed"
      return(row_data)
    }
    
    # Extract expanded model info
    new_aic <- -2 * expanded_fit$logLik + 2 * expanded_fit$d
    improvement <- base_subset_aic - new_aic
    
    row_data$Expanded_Model_AIC <- new_aic
    row_data$AIC_Improvement <- improvement
    row_data$Expanded_LogLik <- expanded_fit$logLik
    row_data$Phylogenetic_Signal_Expanded <- expanded_fit$alpha
    
    # Extract predictor coefficient info
    coef_summary <- summary(expanded_fit)$coefficients
    
    # For categorical variables, look for any row that starts with the variable name
    if (var_type == "categorical") {
      pred_rows <- grep(paste0("^", var_name), rownames(coef_summary))
    } else {
      pred_rows <- which(rownames(coef_summary) == var_name)
    }
    
    if (length(pred_rows) > 0) {
      # For categorical, take the first level or the most significant
      if (length(pred_rows) > 1) {
        # Find most significant p-value
        p_values <- coef_summary[pred_rows, "p.value"]
        pred_row <- pred_rows[which.min(p_values)]
        row_data$Note <- paste0("Cat(", length(pred_rows), " levels)")
      } else {
        pred_row <- pred_rows[1]
      }
      
      row_data$Predictor_Coefficient <- coef_summary[pred_row, "Estimate"]
      row_data$Predictor_SE <- coef_summary[pred_row, "StdErr"]
      row_data$Predictor_p_value <- coef_summary[pred_row, "p.value"]
      
      # Significance
      pval <- row_data$Predictor_p_value
      if (pval < 0.001) {
        sig <- "***"
      } else if (pval < 0.01) {
        sig <- "**"
      } else if (pval < 0.05) {
        sig <- "*"
      } else {
        sig <- "ns"
      }
      if (!("Note" %in% names(row_data)) || is.na(row_data$Note)) {
        row_data$Note <- sig
      }
    }
    
    cat("   Model with", predictor_name, "AIC:", round(new_aic, 2), "\n")
    cat("   TRUE AIC improvement:", round(improvement, 2), "\n")
    
    return(row_data)
  }
  
  # Test each predictor
  predictors_to_test <- list(
    # Note: CB is already in the alternative base model
    list(name = "Migration (1-3)", 
         values = as.numeric(full_data$Migration_AVONET),
         var_name = "Migration_num", var_type = "ordinal"),
    list(name = "Absolute Latitude",
         values = abs(full_data$Centroid.Latitude_AVONET),
         var_name = "abs_Latitude", var_type = "continuous"),
    list(name = "Geographic Region",
         values = full_data$GeographicRegion_Jetz,
         var_name = "GeographicRegion_Jetz", var_type = "categorical"),
    list(name = "Wing Dimorphism",
         values = full_data$PercentAbsLogWingDimorphism,
         var_name = "PercentAbsLogWingDimorphism", var_type = "continuous"),
    list(name = "Plumage Dimorphism",
         values = full_data$logMaleFemalePlumageDiffAbs,
         var_name = "logMaleFemalePlumageDiffAbs", var_type = "continuous"),
    list(name = "Familial Living",
         values = full_data$Griesser2017FamilialLiving,
         var_name = "Griesser2017FamilialLiving", var_type = "binary")
  )
  
  # Test each predictor
  for (pred in predictors_to_test) {
    row_result <- test_predictor(pred$name, pred$values, pred$var_name, pred$var_type)
    if (!is.null(row_result)) {
      comprehensive_table <- rbind(comprehensive_table, row_result)
    }
  }
  
  # Sort table by AIC improvement
  comprehensive_table <- comprehensive_table %>%
    arrange(desc(AIC_Improvement))
  
  # Write results
  write.csv(comprehensive_table, 
            file.path(output_dir, "alternative_base_comparison.csv"),
            row.names = FALSE)
  
  # Create visualization
  if (nrow(comprehensive_table) > 0) {
    p <- ggplot(comprehensive_table, 
                aes(x = reorder(Predictor, AIC_Improvement), 
                    y = AIC_Improvement)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, color = "black") +
      geom_text(aes(label = paste0("n=", N_Species)), 
                hjust = -0.1, size = 3) +
      coord_flip() +
      labs(title = paste0("Stepwise Improvements: ", analysis_name, " (Alternative Base)"),
           subtitle = paste0("Base: ", base_formula_str),
           x = "", y = "AIC Improvement") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "alternative_base_improvements.png"),
           p, width = 10, height = 6, dpi = 300)
  }
  
  # Create summary
  summary_text <- "ALTERNATIVE BASE MODEL STEPWISE EXPANSION\n"
  summary_text <- paste0(summary_text, paste(rep("=", 50), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Analysis: ", analysis_name, "\n")
  summary_text <- paste0(summary_text, "Base model: ", best_model_name, 
                         " (second-best in original comparison)\n")
  summary_text <- paste0(summary_text, "Base AIC: ", round(original_base_aic, 2), "\n")
  summary_text <- paste0(summary_text, "Base formula: ", base_formula_str, "\n\n")
  
  summary_text <- paste0(summary_text, "RESULTS:\n")
  for (i in 1:nrow(comprehensive_table)) {
    summary_text <- paste0(summary_text, 
                          comprehensive_table$Predictor[i], ": ",
                          "AIC improvement = ", round(comprehensive_table$AIC_Improvement[i], 2),
                          " (", comprehensive_table$Note[i], ")\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "alternative_base_summary.txt"))
  
  cat("\n\nResults saved to:", output_dir, "\n")
  
  return(comprehensive_table)
}

# Run the analysis
results <- run_stepwise_with_alternative_base()