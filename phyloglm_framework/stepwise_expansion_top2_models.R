# Enhanced stepwise expansion for top 2 models from each analysis
# Includes tracking of main FS<->CB associations

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

# Helper function to extract FS/CB coefficients and CIs
extract_main_effects <- function(model, analysis_direction) {
  coef_summary <- summary(model)$coefficients
  
  # Determine which coefficient to look for based on direction
  if (grepl("^FS_vs_CB", analysis_direction)) {
    # FS is response, CB is predictor
    main_pred <- "HighConfidence_Coop"
  } else {
    # CB is response, FS is predictor
    main_pred <- "FemaleSong_Agg01"
  }
  
  # Find the main effect (could be in interaction model)
  main_rows <- grep(paste0("^", main_pred, "$"), rownames(coef_summary))
  
  if (length(main_rows) > 0) {
    coef <- coef_summary[main_rows[1], "Estimate"]
    se <- coef_summary[main_rows[1], "StdErr"]
    pval <- coef_summary[main_rows[1], "p.value"]
    
    # Calculate CI
    ci_lower <- coef - 1.96 * se
    ci_upper <- coef + 1.96 * se
    
    # For logistic regression, also calculate OR
    or <- exp(coef)
    or_lower <- exp(ci_lower)
    or_upper <- exp(ci_upper)
    
    return(list(
      coefficient = coef,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = pval,
      or = or,
      or_lower = or_lower,
      or_upper = or_upper
    ))
  } else {
    return(list(
      coefficient = NA,
      se = NA,
      ci_lower = NA,
      ci_upper = NA,
      p_value = NA,
      or = NA,
      or_lower = NA,
      or_upper = NA
    ))
  }
}

run_stepwise_expansion_top2 <- function(
  results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/stepwise_top2_models_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data and results...\n")
  all_results <- readRDS(results_path)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  
  # Results storage
  comprehensive_table <- data.frame()
  
  # Analyze each direction
  for (analysis_name in c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass", 
                          "FS_vs_CB_Terr3_Mass", "CB_vs_FS_Terr3_Mass")) {
    
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("Analyzing:", analysis_name, "\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n")
    
    # Get base analysis results
    base_result <- all_results[[analysis_name]]
    if (is.null(base_result)) {
      cat("Analysis", analysis_name, "not found. Skipping.\n")
      next
    }
    
    # Get top 2 models
    model_comparison <- base_result$comparison$comparison
    top_models <- model_comparison[1:2,]
    
    cat("Top 2 models:\n")
    print(top_models)
    cat("\n")
    
    # Store deltaAIC values for plot labels
    if (!exists("model_deltaAIC")) {
      model_deltaAIC <- data.frame()
    }
    for (i in 1:2) {
      model_deltaAIC <- rbind(model_deltaAIC, data.frame(
        Analysis = analysis_name,
        Model_Rank = i,
        Model_Name = top_models$Model[i],
        deltaAIC = top_models$deltaAIC[i]
      ))
    }
    
    # Process each of the top 2 models
    for (model_idx in 1:2) {
      model_name <- top_models$Model[model_idx]
      model_obj <- base_result$models$models[[model_name]]
      original_aic <- top_models$AIC[model_idx]
      model_rank <- model_idx
      
      cat("\n--- Processing model", model_idx, ":", model_name, "---\n")
      
      base_data <- base_result$prepared_data$data
      base_tree <- base_result$prepared_data$tree
      base_formula <- formula(model_obj)
      base_formula_str <- paste(deparse(base_formula), collapse = " ")
      
      # Extract main effects from base model
      base_main_effects <- extract_main_effects(model_obj, analysis_name)
      
      cat("Base AIC:", round(original_aic, 2), "\n")
      cat("Base formula:", base_formula_str, "\n")
      cat("Main association coefficient:", round(base_main_effects$coefficient, 3), 
          " (p=", round(base_main_effects$p_value, 4), ")\n\n")
      
      # Original sample size
      original_n <- nrow(base_data)
      
      # Test each predictor
      cat("Testing predictors:\n")
      
      # Helper function to test a predictor
      test_predictor <- function(predictor_name, predictor_values, var_name, var_type = "continuous") {
        cat("  ", predictor_name, "... ", sep="")
        
        # Add predictor to base data
        test_data <- base_data
        test_data[[var_name]] <- predictor_values[match(test_data$species, full_data$species)]
        
        # Count missing data
        n_missing_predictor <- sum(is.na(test_data[[var_name]]))
        
        # Remove NAs
        test_data <- test_data[!is.na(test_data[[var_name]]), ]
        
        if (nrow(test_data) < 100) {
          cat("insufficient data (n=", nrow(test_data), ")\n")
          return(NULL)
        }
        
        # Match tree
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        # Initialize row data
        row_data <- data.frame(
          Analysis = analysis_name,
          Model_Rank = model_rank,
          Base_Model_Name = model_name,
          Predictor = predictor_name,
          Base_Formula = base_formula_str,
          Base_Original_AIC = original_aic,
          Base_Original_N = original_n,
          N_Species = nrow(test_data),
          N_Missing_Predictor = n_missing_predictor,
          Percent_Retained = round(100 * nrow(test_data) / original_n, 1),
          stringsAsFactors = FALSE
        )
        
        # Add base model main effects
        row_data$Base_Main_Coef <- base_main_effects$coefficient
        row_data$Base_Main_SE <- base_main_effects$se
        row_data$Base_Main_CI_Lower <- base_main_effects$ci_lower
        row_data$Base_Main_CI_Upper <- base_main_effects$ci_upper
        row_data$Base_Main_P <- base_main_effects$p_value
        row_data$Base_Main_OR <- base_main_effects$or
        
        # CRITICAL: Refit base model to this subset
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
          cat("base refit failed\n")
          return(NULL)
        })
        
        if (is.null(base_subset_fit)) {
          return(NULL)
        }
        
        # Extract base subset info
        base_subset_aic <- -2 * base_subset_fit$logLik + 2 * base_subset_fit$d
        base_subset_main <- extract_main_effects(base_subset_fit, analysis_name)
        
        row_data$Base_Subset_AIC <- base_subset_aic
        row_data$Base_Subset_Main_Coef <- base_subset_main$coefficient
        row_data$Base_Subset_Main_SE <- base_subset_main$se
        row_data$Base_Subset_Main_CI_Lower <- base_subset_main$ci_lower
        row_data$Base_Subset_Main_CI_Upper <- base_subset_main$ci_upper
        row_data$Base_Subset_Main_P <- base_subset_main$p_value
        row_data$Base_Subset_Main_OR <- base_subset_main$or
        
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
          cat("expanded fit failed\n")
          return(NULL)
        })
        
        if (is.null(expanded_fit)) {
          return(NULL)
        }
        
        # Extract expanded model info
        new_aic <- -2 * expanded_fit$logLik + 2 * expanded_fit$d
        improvement <- base_subset_aic - new_aic
        expanded_main <- extract_main_effects(expanded_fit, analysis_name)
        
        row_data$Expanded_Model_AIC <- new_aic
        row_data$AIC_Improvement <- improvement
        
        # Add expanded model main effects
        row_data$Expanded_Main_Coef <- expanded_main$coefficient
        row_data$Expanded_Main_SE <- expanded_main$se
        row_data$Expanded_Main_CI_Lower <- expanded_main$ci_lower
        row_data$Expanded_Main_CI_Upper <- expanded_main$ci_upper
        row_data$Expanded_Main_P <- expanded_main$p_value
        row_data$Expanded_Main_OR <- expanded_main$or
        
        # Calculate change in main effect
        if (!is.na(base_subset_main$coefficient) && !is.na(expanded_main$coefficient)) {
          row_data$Main_Coef_Change <- expanded_main$coefficient - base_subset_main$coefficient
          row_data$Main_Coef_Pct_Change <- 100 * row_data$Main_Coef_Change / abs(base_subset_main$coefficient)
        } else {
          row_data$Main_Coef_Change <- NA
          row_data$Main_Coef_Pct_Change <- NA
        }
        
        # Extract predictor coefficient info
        coef_summary <- summary(expanded_fit)$coefficients
        
        # For categorical variables, look for any row that starts with the variable name
        if (var_type == "categorical") {
          pred_rows <- grep(paste0("^", var_name), rownames(coef_summary))
        } else {
          pred_rows <- which(rownames(coef_summary) == var_name)
        }
        
        if (length(pred_rows) > 0) {
          if (length(pred_rows) > 1) {
            # Find most significant p-value
            p_values <- coef_summary[pred_rows, "p.value"]
            pred_row <- pred_rows[which.min(p_values)]
          } else {
            pred_row <- pred_rows[1]
          }
          
          row_data$Predictor_Coefficient <- coef_summary[pred_row, "Estimate"]
          row_data$Predictor_SE <- coef_summary[pred_row, "StdErr"]
          row_data$Predictor_P <- coef_summary[pred_row, "p.value"]
          
          # Significance
          if (row_data$Predictor_P < 0.001) {
            sig <- "***"
          } else if (row_data$Predictor_P < 0.01) {
            sig <- "**"
          } else if (row_data$Predictor_P < 0.05) {
            sig <- "*"
          } else {
            sig <- "ns"
          }
          row_data$Predictor_Sig <- sig
        }
        
        cat("AIC improvement =", round(improvement, 1), "\n")
        
        return(row_data)
      }
      
      # Define predictors to test
      predictors_to_test <- list()
      
      # Avoid adding predictors already in the model
      existing_vars <- all.vars(base_formula)
      
      if (!("Territory_num" %in% existing_vars) && !("Territory_12vs3" %in% existing_vars) && 
          !("Territory_12vs31" %in% existing_vars) && !grepl("Terr3", model_name)) {
        predictors_to_test[["Territory"]] <- list(
          name = "Territory (1-3)",
          values = as.numeric(full_data$Territory),
          var_name = "Territory_num",
          var_type = "ordinal"
        )
      }
      
      if (!("Migration_num" %in% existing_vars)) {
        predictors_to_test[["Migration"]] <- list(
          name = "Migration (1-3)",
          values = as.numeric(full_data$Migration_AVONET),
          var_name = "Migration_num",
          var_type = "ordinal"
        )
      }
      
      if (!("abs_Latitude" %in% existing_vars)) {
        predictors_to_test[["Latitude"]] <- list(
          name = "Absolute Latitude",
          values = abs(full_data$Centroid.Latitude_AVONET),
          var_name = "abs_Latitude",
          var_type = "continuous"
        )
      }
      
      if (!("GeographicRegion_Jetz" %in% existing_vars)) {
        predictors_to_test[["Region"]] <- list(
          name = "Geographic Region",
          values = full_data$GeographicRegion_Jetz,
          var_name = "GeographicRegion_Jetz",
          var_type = "categorical"
        )
      }
      
      if (!("PercentAbsLogWingDimorphism" %in% existing_vars)) {
        predictors_to_test[["Wing"]] <- list(
          name = "Wing Dimorphism",
          values = full_data$PercentAbsLogWingDimorphism,
          var_name = "PercentAbsLogWingDimorphism",
          var_type = "continuous"
        )
      }
      
      if (!("logMaleFemalePlumageDiffAbs" %in% existing_vars)) {
        predictors_to_test[["Plumage"]] <- list(
          name = "Plumage Dimorphism",
          values = full_data$logMaleFemalePlumageDiffAbs,
          var_name = "logMaleFemalePlumageDiffAbs",
          var_type = "continuous"
        )
      }
      
      if (!("Griesser2017FamilialLiving" %in% existing_vars)) {
        predictors_to_test[["FamilialLiving"]] <- list(
          name = "Familial Living",
          values = full_data$Griesser2017FamilialLiving,
          var_name = "Griesser2017FamilialLiving",
          var_type = "binary"
        )
      }
      
      # Test each predictor
      for (pred in predictors_to_test) {
        row_result <- test_predictor(pred$name, pred$values, pred$var_name, pred$var_type)
        if (!is.null(row_result)) {
          comprehensive_table <- rbind(comprehensive_table, row_result)
        }
      }
    }
  }
  
  # Sort table
  comprehensive_table <- comprehensive_table %>%
    arrange(Analysis, Model_Rank, desc(AIC_Improvement))
  
  # Write comprehensive table
  write.csv(comprehensive_table, 
            file.path(output_dir, "top2_models_comprehensive_table.csv"),
            row.names = FALSE)
  
  # Create summary by analysis and model rank
  summary_stats <- comprehensive_table %>%
    group_by(Analysis, Model_Rank, Base_Model_Name) %>%
    summarise(
      N_Predictors_Tested = n(),
      N_Improvements = sum(AIC_Improvement > aic_threshold, na.rm = TRUE),
      Best_Predictor = Predictor[which.max(AIC_Improvement)],
      Best_Improvement = max(AIC_Improvement, na.rm = TRUE),
      Mean_Main_Coef_Change = mean(Main_Coef_Pct_Change, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(summary_stats,
            file.path(output_dir, "top2_models_summary.csv"),
            row.names = FALSE)
  
  # Create visualization of main effect changes
  plot_data <- comprehensive_table %>%
    filter(!is.na(Main_Coef_Change)) %>%
    left_join(model_deltaAIC, by = c("Analysis", "Model_Rank" = "Model_Rank")) %>%
    mutate(
      Direction = ifelse(grepl("^FS_vs_CB", Analysis), "FS→CB", "CB→FS"),
      Model_Label = paste0("Model ", Model_Rank, ": ", Base_Model_Name)
    )
  
  if (nrow(plot_data) > 0) {
    p1 <- ggplot(plot_data, 
                 aes(x = reorder(Predictor, AIC_Improvement), 
                     y = Main_Coef_Pct_Change,
                     fill = factor(Model_Rank))) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      facet_wrap(~ Direction, scales = "free_y") +
      coord_flip() +
      scale_fill_brewer(palette = "Set1", name = "Model Rank") +
      labs(title = "Change in Main Association When Adding Predictors",
           subtitle = "Percentage change in FS↔CB coefficient",
           x = "", y = "% Change in Main Effect") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "main_effect_changes.png"),
           p1, width = 12, height = 8, dpi = 300)
    
    # Create plot showing confidence intervals with improved layout
    plot_data_p2 <- comprehensive_table %>%
      filter(!is.na(Base_Subset_Main_Coef)) %>%
      left_join(model_deltaAIC, by = c("Analysis", "Model_Rank" = "Model_Rank")) %>%
      mutate(
        Direction = ifelse(grepl("^FS_vs_CB", Analysis), "CB→FS", "FS→CB"),
        Model_Label = paste0(Analysis, " - ", Base_Model_Name, " (ΔAIC=", round(deltaAIC, 1), ")"),
        Predictor_Bold = AIC_Improvement > aic_threshold,
        Analysis_Order = factor(Analysis, levels = c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass",
                                                    "FS_vs_CB_Terr3_Mass", "CB_vs_FS_Terr3_Mass"))
      ) %>%
      arrange(Analysis_Order, desc(Predictor))
    
    # Create separate plots for each model rank
    p2_rank1 <- plot_data_p2 %>%
      filter(Model_Rank == 1) %>%
      ggplot(aes(x = reorder(Predictor, desc(Predictor)))) +
      geom_point(aes(y = Base_Subset_Main_Coef), shape = 1, size = 3) +
      geom_errorbar(aes(ymin = Base_Subset_Main_CI_Lower, 
                       ymax = Base_Subset_Main_CI_Upper),
                   width = 0.2, alpha = 0.5) +
      geom_point(aes(y = Expanded_Main_Coef), shape = 16, size = 3) +
      geom_errorbar(aes(ymin = Expanded_Main_CI_Lower, 
                       ymax = Expanded_Main_CI_Upper),
                   width = 0.2) +
      facet_wrap(~ Model_Label, scales = "free", ncol = 1, 
                 labeller = label_wrap_gen(width = 50)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      coord_flip() +
      labs(x = "", y = "Coefficient", title = "Best Models (Model Rank = 1)") +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 9),
        axis.text.y = element_text(face = c("plain", "bold")[as.numeric(plot_data_p2 %>% 
                                                                        filter(Model_Rank == 1) %>% 
                                                                        arrange(Analysis_Order, desc(Predictor)) %>% 
                                                                        pull(Predictor_Bold)) + 1]),
        panel.spacing.y = unit(0.5, "lines")
      )
    
    p2_rank2 <- plot_data_p2 %>%
      filter(Model_Rank == 2) %>%
      ggplot(aes(x = reorder(Predictor, desc(Predictor)))) +
      geom_point(aes(y = Base_Subset_Main_Coef), shape = 1, size = 3) +
      geom_errorbar(aes(ymin = Base_Subset_Main_CI_Lower, 
                       ymax = Base_Subset_Main_CI_Upper),
                   width = 0.2, alpha = 0.5) +
      geom_point(aes(y = Expanded_Main_Coef), shape = 16, size = 3) +
      geom_errorbar(aes(ymin = Expanded_Main_CI_Lower, 
                       ymax = Expanded_Main_CI_Upper),
                   width = 0.2) +
      facet_wrap(~ Model_Label, scales = "free", ncol = 1,
                 labeller = label_wrap_gen(width = 50)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      coord_flip() +
      labs(x = "", y = "Coefficient", title = "Second-Best Models (Model Rank = 2)") +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 9),
        axis.text.y = element_text(face = c("plain", "bold")[as.numeric(plot_data_p2 %>% 
                                                                        filter(Model_Rank == 2) %>% 
                                                                        arrange(Analysis_Order, desc(Predictor)) %>% 
                                                                        pull(Predictor_Bold)) + 1]),
        panel.spacing.y = unit(0.5, "lines")
      )
    
    # Function to add custom x-axis labels to each facet
    add_custom_labels <- function(plot_obj, data_subset) {
      plot_build <- ggplot_build(plot_obj)
      plot_table <- ggplot_gtable(plot_build)
      
      # Get panel info
      panels <- data_subset %>%
        select(Model_Label, Direction) %>%
        distinct() %>%
        arrange(Model_Label)
      
      # Find axis labels and modify
      for (i in seq_len(nrow(panels))) {
        # Add custom label
        axis_label <- grid::textGrob(panels$Direction[i], 
                                   gp = grid::gpar(fontsize = 9, col = "gray40"))
        
        # Find the correct position for this panel
        panel_name <- paste0("panel-", i, "-1")
        panel_pos <- which(plot_table$layout$name == panel_name)
        
        if (length(panel_pos) > 0) {
          # Get coordinates
          l <- plot_table$layout$l[panel_pos]
          r <- plot_table$layout$r[panel_pos]
          t <- max(plot_table$layout$t) + 1
          
          # Add the label
          plot_table <- gtable::gtable_add_rows(plot_table, unit(0.5, "lines"), t-1)
          plot_table <- gtable::gtable_add_grob(plot_table, axis_label, t = t, l = l, r = r)
        }
      }
      
      return(plot_table)
    }
    
    # Apply custom labels
    p2_rank1_labeled <- add_custom_labels(p2_rank1, filter(plot_data_p2, Model_Rank == 1))
    p2_rank2_labeled <- add_custom_labels(p2_rank2, filter(plot_data_p2, Model_Rank == 2))
    
    # Combine plots using gridExtra
    library(gridExtra)
    p2_combined <- arrangeGrob(
      grobs = list(p2_rank1_labeled, p2_rank2_labeled),
      ncol = 2,
      top = grid::textGrob("Main Association Coefficients: Base (○) vs Expanded (●) Models", 
                          gp = grid::gpar(fontsize = 16, fontface = "bold")),
      bottom = grid::textGrob("With 95% confidence intervals. Bold parameter names indicate AIC improvement > 2", 
                             gp = grid::gpar(fontsize = 10)),
      left = grid::textGrob("Parameter added to expanded model", rot = 90, 
                           gp = grid::gpar(fontsize = 12))
    )
    
    ggsave(file.path(output_dir, "main_coefficients_comparison.png"),
           p2_combined, width = 14, height = 10, dpi = 300)
    
    # Create comprehensive improvements plot (like the one you showed)
    p3 <- comprehensive_table %>%
      filter(!is.na(AIC_Improvement)) %>%
      mutate(
        Analysis_Clean = case_when(
          Analysis == "CB_vs_FS_Terr3_Mass" ~ "CB_vs_FS_Terr3_Mass",
          Analysis == "CB_vs_FS_TerrWS_Mass" ~ "CB_vs_FS_TerrWS_Mass", 
          Analysis == "FS_vs_CB_Terr3_Mass" ~ "FS_vs_CB_Terr3_Mass",
          Analysis == "FS_vs_CB_TerrWS_Mass" ~ "FS_vs_CB_TerrWS_Mass"
        ),
        Model_Type = paste0("Model ", Model_Rank)
      ) %>%
      ggplot(aes(x = AIC_Improvement, 
                 y = reorder(Predictor, AIC_Improvement),
                 fill = Analysis_Clean)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_vline(xintercept = aic_threshold, linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, color = "black") +
      geom_text(aes(label = paste0("n=", N_Species),
                    hjust = ifelse(AIC_Improvement > 0, -0.1, 1.1)), 
                position = position_dodge(width = 0.9),
                size = 3) +
      scale_fill_manual(values = c("CB_vs_FS_Terr3_Mass" = "#1b9e77",
                                   "CB_vs_FS_TerrWS_Mass" = "#d95f02",
                                   "FS_vs_CB_Terr3_Mass" = "#7570b3",
                                   "FS_vs_CB_TerrWS_Mass" = "#e7298a")) +
      labs(title = "Comprehensive Stepwise Model Improvements",
           subtitle = "AIC improvements calculated on same species subsets",
           x = "AIC Improvement", y = "",
           fill = "Analysis") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      xlim(-60, 65)
    
    ggsave(file.path(output_dir, "comprehensive_improvements_plot.png"),
           p3, width = 12, height = 8, dpi = 300)
  }
  
  # Create summary text file
  summary_text <- "TOP 2 MODELS STEPWISE EXPANSION RESULTS\n"
  summary_text <- paste0(summary_text, paste(rep("=", 70), collapse=""), "\n")
  summary_text <- paste0(summary_text, "Generated: ", Sys.Date(), "\n")
  summary_text <- paste0(summary_text, "AIC improvement threshold: ", aic_threshold, "\n\n")
  
  for (analysis in unique(comprehensive_table$Analysis)) {
    summary_text <- paste0(summary_text, "\n", analysis, "\n")
    summary_text <- paste0(summary_text, paste(rep("-", 40), collapse=""), "\n")
    
    analysis_data <- comprehensive_table[comprehensive_table$Analysis == analysis,]
    
    for (rank in 1:2) {
      rank_data <- analysis_data[analysis_data$Model_Rank == rank,]
      if (nrow(rank_data) > 0) {
        model_name <- rank_data$Base_Model_Name[1]
        summary_text <- paste0(summary_text, "\nModel ", rank, ": ", model_name, "\n")
        summary_text <- paste0(summary_text, "Base AIC: ", 
                              round(rank_data$Base_Original_AIC[1], 2), "\n")
        
        # Show predictors with improvements > threshold
        good_preds <- rank_data[rank_data$AIC_Improvement > aic_threshold,]
        if (nrow(good_preds) > 0) {
          summary_text <- paste0(summary_text, "Significant improvements:\n")
          for (i in 1:nrow(good_preds)) {
            summary_text <- paste0(summary_text, "  - ", good_preds$Predictor[i], 
                                  ": AIC improvement = ", 
                                  round(good_preds$AIC_Improvement[i], 2),
                                  ", Main effect change = ",
                                  round(good_preds$Main_Coef_Pct_Change[i], 1), "%\n")
          }
        } else {
          summary_text <- paste0(summary_text, "No predictors improved AIC > ", 
                                aic_threshold, "\n")
        }
      }
    }
  }
  
  writeLines(summary_text, file.path(output_dir, "top2_models_summary.txt"))
  
  cat("\n\nResults saved to:", output_dir, "\n")
  cat("Main output: top2_models_comprehensive_table.csv\n")
  
  return(comprehensive_table)
}

# Run the analysis
results <- run_stepwise_expansion_top2()