# Enhanced stepwise expansion with comprehensive output table
# This version creates a detailed comparison table with all necessary information

library(phylolm)
library(ape)
library(dplyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

# Helper function to safely extract model information
safe_extract <- function(value, default = NA) {
  if (is.null(value) || length(value) == 0) return(default)
  return(value)
}

run_stepwise_expansion_comprehensive_table <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_comprehensive_table_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Results storage
  expansion_results <- list()
  comprehensive_table <- data.frame()
  
  # Analyze each direction
  for (analysis_name in c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass", 
                          "FS_vs_CB_Terr3_Mass", "CB_vs_FS_Terr3_Mass")) {
    
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("Analyzing:", analysis_name, "\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n")
    
    # Get base model
    base_result <- all_results[[analysis_name]]
    if (is.null(base_result)) {
      cat("Analysis", analysis_name, "not found. Skipping.\n")
      next
    }
    
    best_model_name <- base_result$comparison$comparison$Model[1]
    best_model <- base_result$models$models[[best_model_name]]
    original_base_aic <- base_result$comparison$comparison$AIC[1]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    base_formula <- formula(best_model)
    base_formula_str <- paste(deparse(base_formula), collapse = " ")
    
    cat("Base model:", best_model_name, "\n")
    cat("Original AIC (", nrow(base_data), " species):", round(original_base_aic, 2), "\n")
    cat("Base formula:", base_formula_str, "\n\n")
    
    # Original sample size
    original_n <- nrow(base_data)
    
    # Test each predictor individually
    cat("Testing predictors (with proper subsetting):\n")
    cat(paste(rep("-", 40), collapse=""), "\n")
    
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
        
        # Still add to table with NA values
        row_data <- data.frame(
          Analysis = analysis_name,
          Predictor = predictor_name,
          Base_Model_Name = best_model_name,
          Base_Formula = base_formula_str,
          Base_Original_AIC = original_base_aic,
          Base_Original_N = original_n,
          N_Species = nrow(test_data),
          N_Missing_Predictor = n_missing_predictor,
          Percent_Retained = round(100 * nrow(test_data) / original_n, 1),
          Base_Subset_AIC = NA,
          Predictor_Variable_Name = var_name,
          Predictor_Type = var_type,
          Expanded_Formula = NA,
          Expanded_Model_AIC = NA,
          AIC_Improvement = NA,
          Base_LogLik = NA,
          Expanded_LogLik = NA,
          LRT_Statistic = NA,
          LRT_df = NA,
          LRT_p_value = NA,
          Predictor_Coefficient = NA,
          Predictor_SE = NA,
          Predictor_z_value = NA,
          Predictor_p_value = NA,
          Predictor_OR = NA,
          Predictor_OR_CI_lower = NA,
          Predictor_OR_CI_upper = NA,
          Phylogenetic_Signal_Base = NA,
          Phylogenetic_Signal_Expanded = NA,
          Convergence_Base = FALSE,
          Convergence_Expanded = FALSE,
          Note = "Insufficient data",
          stringsAsFactors = FALSE
        )
        
        return(row_data)
      }
      
      # Match tree
      test_tree <- keep.tip(base_tree, test_data$species)
      test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
      
      cat("   N species with", predictor_name, "data:", nrow(test_data), "\n")
      
      # Initialize row data
      row_data <- data.frame(
        Analysis = analysis_name,
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
      row_data$Phylogenetic_Signal_Base <- safe_extract(base_subset_fit$alpha, NA)
      row_data$Convergence_Base <- !is.null(base_subset_fit$convergence) && base_subset_fit$convergence
      
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
      row_data$Phylogenetic_Signal_Expanded <- safe_extract(expanded_fit$alpha, NA)
      row_data$Convergence_Expanded <- !is.null(expanded_fit$convergence) && expanded_fit$convergence
      
      # Likelihood ratio test
      lrt_stat <- 2 * (expanded_fit$logLik - base_subset_fit$logLik)
      lrt_df <- expanded_fit$d - base_subset_fit$d
      lrt_p <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
      
      row_data$LRT_Statistic <- lrt_stat
      row_data$LRT_df <- lrt_df
      row_data$LRT_p_value <- lrt_p
      
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
        row_data$Predictor_z_value <- coef_summary[pred_row, "z.value"]
        row_data$Predictor_p_value <- coef_summary[pred_row, "p.value"]
        
        # Calculate odds ratio
        or <- exp(row_data$Predictor_Coefficient)
        or_ci_lower <- exp(row_data$Predictor_Coefficient - 1.96 * row_data$Predictor_SE)
        or_ci_upper <- exp(row_data$Predictor_Coefficient + 1.96 * row_data$Predictor_SE)
        
        row_data$Predictor_OR <- or
        row_data$Predictor_OR_CI_lower <- or_ci_lower
        row_data$Predictor_OR_CI_upper <- or_ci_upper
      } else {
        # Set defaults if predictor not found
        row_data$Predictor_Coefficient <- NA
        row_data$Predictor_SE <- NA
        row_data$Predictor_z_value <- NA
        row_data$Predictor_p_value <- NA
        row_data$Predictor_OR <- NA
        row_data$Predictor_OR_CI_lower <- NA
        row_data$Predictor_OR_CI_upper <- NA
      }
      
      cat("   Model with", predictor_name, "AIC:", round(new_aic, 2), "\n")
      cat("   TRUE AIC improvement:", round(improvement, 2), "\n")
      
      # Add note about significance (if not already set for categorical)
      if (!("Note" %in% names(row_data)) || is.na(row_data$Note[1])) {
        if (length(row_data$Predictor_p_value) > 0 && !is.na(row_data$Predictor_p_value[1])) {
          pval <- row_data$Predictor_p_value[1]
          if (pval < 0.001) {
            row_data$Note <- "***"
          } else if (pval < 0.01) {
            row_data$Note <- "**"
          } else if (pval < 0.05) {
            row_data$Note <- "*"
          } else {
            row_data$Note <- "ns"
          }
        } else {
          row_data$Note <- ""
        }
      }
      
      return(row_data)
    }
    
    # Test each predictor
    predictors_to_test <- list()
    
    # 1. Territory as numeric (if not already in as Terr3)
    if (!grepl("Terr3", analysis_name)) {
      predictors_to_test[["Territory"]] <- list(
        name = "Territory (1-3)",
        values = as.numeric(full_data$Territory),
        var_name = "Territory_num",
        var_type = "ordinal"
      )
    }
    
    # 2. Migration
    predictors_to_test[["Migration"]] <- list(
      name = "Migration (1-3)",
      values = as.numeric(full_data$Migration_AVONET),
      var_name = "Migration_num",
      var_type = "ordinal"
    )
    
    # 3. Absolute Latitude
    if (!grepl("absLat|Latitude", analysis_name)) {
      predictors_to_test[["Latitude"]] <- list(
        name = "Absolute Latitude",
        values = abs(full_data$Centroid.Latitude_AVONET),
        var_name = "abs_Latitude",
        var_type = "continuous"
      )
    }
    
    # 4. Geographic Region
    if (!grepl("Region", analysis_name)) {
      predictors_to_test[["Region"]] <- list(
        name = "Geographic Region",
        values = full_data$GeographicRegion_Jetz,
        var_name = "GeographicRegion_Jetz",
        var_type = "categorical"
      )
    }
    
    # 5. Wing Dimorphism
    if (!grepl("WingDim", analysis_name)) {
      predictors_to_test[["Wing"]] <- list(
        name = "Wing Dimorphism",
        values = full_data$PercentAbsLogWingDimorphism,
        var_name = "PercentAbsLogWingDimorphism",
        var_type = "continuous"
      )
    }
    
    # 6. Plumage Dimorphism
    if (!grepl("PlumDim", analysis_name)) {
      predictors_to_test[["Plumage"]] <- list(
        name = "Plumage Dimorphism",
        values = full_data$logMaleFemalePlumageDiffAbs,
        var_name = "logMaleFemalePlumageDiffAbs",
        var_type = "continuous"
      )
    }
    
    # 7. Familial Living
    predictors_to_test[["FamilialLiving"]] <- list(
      name = "Familial Living",
      values = full_data$Griesser2017FamilialLiving,
      var_name = "Griesser2017FamilialLiving",
      var_type = "binary"
    )
    
    # Test each predictor
    for (pred in predictors_to_test) {
      row_result <- test_predictor(pred$name, pred$values, pred$var_name, pred$var_type)
      comprehensive_table <- rbind(comprehensive_table, row_result)
    }
  }
  
  # Sort table by analysis and AIC improvement
  comprehensive_table <- comprehensive_table %>%
    arrange(Analysis, desc(AIC_Improvement))
  
  # Write comprehensive table
  write.csv(comprehensive_table, 
            file.path(output_dir, "comprehensive_improvement_comparison_table.csv"),
            row.names = FALSE)
  
  # Create summary statistics
  summary_stats <- comprehensive_table %>%
    group_by(Predictor) %>%
    summarise(
      Mean_AIC_Improvement = mean(AIC_Improvement, na.rm = TRUE),
      SD_AIC_Improvement = sd(AIC_Improvement, na.rm = TRUE),
      Min_AIC_Improvement = min(AIC_Improvement, na.rm = TRUE),
      Max_AIC_Improvement = max(AIC_Improvement, na.rm = TRUE),
      N_Significant = sum(Note %in% c("*", "**", "***"), na.rm = TRUE),
      Mean_N_Species = mean(N_Species, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_AIC_Improvement))
  
  write.csv(summary_stats,
            file.path(output_dir, "predictor_summary_statistics.csv"),
            row.names = FALSE)
  
  # Create visualization of comprehensive results
  plot_data <- comprehensive_table %>%
    filter(!is.na(AIC_Improvement))
  
  if (nrow(plot_data) > 0) {
    p1 <- ggplot(plot_data, 
                 aes(x = reorder(Predictor, AIC_Improvement), 
                     y = AIC_Improvement,
                     fill = Analysis)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, color = "black") +
      geom_text(aes(label = paste0("n=", N_Species)), 
                position = position_dodge(width = 0.9),
                hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Comprehensive Stepwise Model Improvements",
           subtitle = "AIC improvements calculated on same species subsets",
           x = "", y = "AIC Improvement",
           fill = "Analysis") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "comprehensive_improvements_plot.png"),
           p1, width = 12, height = 8, dpi = 300)
    
    # Create plot showing sample size effects
    p2 <- ggplot(plot_data,
                 aes(x = N_Species, y = AIC_Improvement, 
                     color = Predictor, shape = Analysis)) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = aic_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "AIC Improvement vs Sample Size",
           subtitle = "Showing relationship between data availability and model improvement",
           x = "Number of Species", y = "AIC Improvement") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "sample_size_vs_improvement.png"),
           p2, width = 10, height = 6, dpi = 300)
  }
  
  # Create summary text file
  summary_text <- "COMPREHENSIVE STEPWISE EXPANSION RESULTS\n"
  summary_text <- paste0(summary_text, paste(rep("=", 70), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Generated: ", Sys.Date(), "\n")
  summary_text <- paste0(summary_text, "AIC improvement threshold: ", aic_threshold, "\n\n")
  
  summary_text <- paste0(summary_text, "SUMMARY STATISTICS BY PREDICTOR\n")
  summary_text <- paste0(summary_text, paste(rep("-", 40), collapse=""), "\n\n")
  
  for (i in 1:nrow(summary_stats)) {
    summary_text <- paste0(summary_text, summary_stats$Predictor[i], ":\n")
    summary_text <- paste0(summary_text, "  Mean AIC improvement: ", 
                          round(summary_stats$Mean_AIC_Improvement[i], 2), 
                          " (SD = ", round(summary_stats$SD_AIC_Improvement[i], 2), ")\n")
    summary_text <- paste0(summary_text, "  Range: ", 
                          round(summary_stats$Min_AIC_Improvement[i], 2), " to ",
                          round(summary_stats$Max_AIC_Improvement[i], 2), "\n")
    summary_text <- paste0(summary_text, "  Significant in ", 
                          summary_stats$N_Significant[i], " analyses\n")
    summary_text <- paste0(summary_text, "  Mean N species: ", 
                          round(summary_stats$Mean_N_Species[i], 0), "\n\n")
  }
  
  writeLines(summary_text, file.path(output_dir, "comprehensive_summary.txt"))
  
  cat("\nComprehensive results saved to:", output_dir, "\n")
  cat("Main output: comprehensive_improvement_comparison_table.csv\n")
  
  return(comprehensive_table)
}

# Run the comprehensive analysis
results <- run_stepwise_expansion_comprehensive_table()