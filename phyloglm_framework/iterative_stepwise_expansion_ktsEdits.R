# Iterative stepwise forward selection
# Builds on initial stepwise results to create optimal models
# KTS note: some very incorrect-seeming loops in run_iterative_stepwise().

library(phylolm)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("phyloglm_framework/batch_runner_helpers.R")

# Helper function to fit a model and extract key info
fit_and_extract <- function(formula, data, tree, method = "logistic_MPLE") {
  tryCatch({
    fit <- phyloglm(
      formula = formula,
      data = data,
      phy = tree,
      method = method,
      btol = 50,
      log.alpha.bound = 4
    )
    
    # Extract key metrics
    aic <- -2 * fit$logLik + 2 * fit$d
    
    # Get main effect coefficient if applicable
    coef_summary <- summary(fit)$coefficients
    
    # Determine main predictor based on response variable
    response_var <- all.vars(formula)[1]
    if (response_var == "FemaleSong_Agg01") {
      main_pred <- "HighConfidence_Coop"
    } else {
      main_pred <- "FemaleSong_Agg01"
    }
    
    # Extract main effect
    main_rows <- grep(paste0("^", main_pred, "$"), rownames(coef_summary))
    if (length(main_rows) > 0) {
      main_coef <- coef_summary[main_rows[1], "Estimate"]
      main_se <- coef_summary[main_rows[1], "StdErr"]
      main_p <- coef_summary[main_rows[1], "p.value"]
    } else {
      main_coef <- NA
      main_se <- NA
      main_p <- NA
    }
    
    return(list(
      success = TRUE,
      aic = aic,
      logLik = fit$logLik,
      alpha = fit$alpha,
      main_coef = main_coef,
      main_se = main_se,
      main_p = main_p,
      coefficients = coef_summary,
      fit = fit
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message
    ))
  })
}

# Main function for iterative stepwise expansion
run_iterative_stepwise <- function(
  initial_results_file = "Outputs/PhyloglmResults/stepwise_top2_models_20250715/top2_models_comprehensive_table1.csv",
  data_path = "Data_R_2025-06-09.csv",
  tree_path = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  all_results_path = "Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds",
  output_dir = NULL,
  aic_threshold = 2
) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("Outputs/PhyloglmResults/iterative_stepwise_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data...\n")
  initial_results <- read.csv(initial_results_file)
  full_data <- read.csv(data_path)
  tree <- read.nexus(tree_path)
  all_results <- readRDS(all_results_path)
  
  # Identify starting points (models with at least one significant improvement)
  starting_points <- initial_results %>%
    filter(AIC_Improvement > aic_threshold & Model_Rank == 1 & str_detect(Analysis, "Terr3", negate = T)) %>%
    group_by(Analysis, Model_Rank, Base_Model_Name) %>%
    summarise(
      Best_Predictor = Predictor[which.max(AIC_Improvement)],
      Best_Variable = Predictor_Variable_Name[which.max(AIC_Improvement)],
      Best_Improvement = max(AIC_Improvement),
      N_Significant = sum(AIC_Improvement > aic_threshold),
      Significant_Predictors = paste(Predictor[AIC_Improvement > aic_threshold], collapse = ", "),
      .groups = "drop"
    )
  
  cat("\nIdentified", nrow(starting_points), "starting points for iterative expansion\n\n")
  
  # Store all iterative results
  all_iterations <- list()
  final_models <- list()
  
  # Process each starting point
  for (i in 1:nrow(starting_points)) {
    sp <- starting_points[i,]
    
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat("Processing:", sp$Analysis, "- Model", sp$Model_Rank, "(", sp$Base_Model_Name, ")\n")
    cat("Starting with:", sp$Best_Predictor, "(AIC improvement:", round(sp$Best_Improvement, 2), ")\n")
    cat("Other candidates:", sp$Significant_Predictors, "\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n")
    
    # Get base model info
    base_result <- all_results[[sp$Analysis]]
    base_model <- base_result$models$models[[sp$Base_Model_Name]]
    base_data <- base_result$prepared_data$data
    base_tree <- base_result$prepared_data$tree
    base_formula <- formula(base_model)
    
    # Get list of candidate predictors
    candidates <- initial_results %>%
      filter(Analysis == sp$Analysis, 
             Model_Rank == sp$Model_Rank,
             AIC_Improvement > aic_threshold) %>%
      arrange(desc(AIC_Improvement)) %>%
      select(Predictor, Predictor_Variable_Name, AIC_Improvement)
    
    # Initialize tracking
    current_formula <- base_formula
    current_predictors <- character(0)
    iteration_results <- list()
    iter <- 0
    
    # Iterative addition
    continue_adding <- FALSE
    
    while (continue_adding && nrow(candidates) > 0) {
      iter <- iter + 1
      cat("\n--- Iteration", iter, "---\n")
      
      # Test each remaining candidate
      test_results <- data.frame()
      
      for (j in 1:nrow(candidates)) {
        cand <- candidates[j,]
        cat("Testing addition of", cand$Predictor, "... ")
        
        # Add all current predictors plus the candidate to data
        test_data <- base_data
        
        # Add all predictors we've already selected
        if (length(current_predictors) > 0) {
          for (k in 1:length(current_predictors)) {
            # Find the variable name from the iteration results
            prev_var <- iteration_results[[k]]$added_variable
            if (!is.null(prev_var) && prev_var %in% names(full_data)) {
              test_data[[prev_var]] <- full_data[[prev_var]][match(test_data$species, full_data$species)]
            }
          }
        }
        
        # Add the candidate predictor
        pred_values <- full_data[[cand$Predictor_Variable_Name]][match(test_data$species, full_data$species)]
        test_data[[cand$Predictor_Variable_Name]] <- pred_values
        
        # Remove NAs from all predictors
        # complete_mask <- complete.cases(test_data[, c(names(base_data), 
        #                sapply(iteration_results, function(x) x$added_variable), cand$Predictor_Variable_Name)]) # this was removing all data
        if (!is.null(iteration_results)) {
          complete_mask <- complete.cases(test_data[, c(sapply(iteration_results, function(x) x$added_variable), cand$Predictor_Variable_Name)])
        } else {
          complete_mask <- complete.cases(test_data[, c(cand$Predictor_Variable_Name)])
        }
        test_data <- test_data[complete_mask, ]
        
        if (nrow(test_data) < 100) {
          cat("insufficient data (n=", nrow(test_data), ")\n")
          next
        }
        
        # Match tree
        test_tree <- keep.tip(base_tree, test_data$species)
        test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
        
        # Update formula
        test_formula <- update(current_formula, paste("~ . +", cand$Predictor_Variable_Name))
        
        # Fit current model on this subset (for fair comparison)
        current_fit <- fit_and_extract(current_formula, test_data, test_tree)
        
        if (!current_fit$success) {
          cat("current model failed\n")
          next
        }
        
        # Fit expanded model
        expanded_fit <- fit_and_extract(test_formula, test_data, test_tree)
        
        if (!expanded_fit$success) {
          cat("expanded model failed\n")
          next
        }
        
        # Calculate improvement
        improvement <- current_fit$aic - expanded_fit$aic
        
        if (improvement > 2) {
          current_predictors <- c(current_predictors, cand$Predictor_Variable_Name)
        }
        
        cat("AIC improvement =", round(improvement, 2), "\n")
        
        # Store results
        test_results <- rbind(test_results, data.frame(
          Predictor = cand$Predictor,
          Variable = cand$Predictor_Variable_Name,
          N_Species = nrow(test_data),
          Current_AIC = current_fit$aic,
          Expanded_AIC = expanded_fit$aic,
          AIC_Improvement = improvement,
          Current_Main_Coef = current_fit$main_coef,
          Expanded_Main_Coef = expanded_fit$main_coef,
          Main_Coef_Change = expanded_fit$main_coef - current_fit$main_coef
        ))
      }
      
      # Check if any predictor improves model
      if (nrow(test_results) == 0 || max(test_results$AIC_Improvement, na.rm = TRUE) <= aic_threshold) {
        cat("\nNo predictor improves model by >", aic_threshold, "AIC units. Stopping.\n")
        continue_adding <- FALSE
      } else {
        # Select best predictor
        best_idx <- which.max(test_results$AIC_Improvement)
        best_pred <- test_results[best_idx,]
        
        cat("\nBest predictor:", best_pred$Predictor, 
            "(AIC improvement:", round(best_pred$AIC_Improvement, 2), ")\n")
        cat("Adding to model...\n")
        
        # Update current model
        current_formula <- update(current_formula, paste("~ . +", best_pred$Variable))
        current_predictors <- c(current_predictors, best_pred$Predictor)
        
        # Remove selected predictor from candidates
        candidates <- candidates[candidates$Predictor != best_pred$Predictor,]
        
        # Store iteration results
        iteration_results[[iter]] <- list(
          iteration = iter,
          added_predictor = best_pred$Predictor,
          added_variable = best_pred$Variable,
          aic_improvement = best_pred$AIC_Improvement,
          cumulative_predictors = current_predictors,
          test_results = test_results,
          selected = best_pred
        )
      }
    }
    
    # Store final model info
    model_key <- paste0(sp$Analysis, "_Model", sp$Model_Rank)
    
    all_iterations[[model_key]] <- iteration_results
    
    final_models[[model_key]] <- list(
      analysis = sp$Analysis,
      model_rank = sp$Model_Rank,
      base_model_name = sp$Base_Model_Name,
      base_formula = base_formula,
      final_formula = current_formula,
      predictors_added = current_predictors,
      n_iterations = iter,
      iteration_details = iteration_results
    )
    
    cat("\nFinal model includes:", length(current_predictors), "additional predictors\n")
    cat("Order of addition:", paste(current_predictors, collapse = " → "), "\n")
  }
  
  # Create summary table
  summary_table <- data.frame()
  
  for (model_key in names(final_models)) {
    fm <- final_models[[model_key]]
    
    # Create summary row
    summary_row <- data.frame(
      Analysis = fm$analysis,
      Model_Rank = fm$model_rank,
      Base_Model = fm$base_model_name,
      N_Predictors_Added = length(fm$predictors_added),
      Predictors_Added = paste(fm$predictors_added, collapse = " + "),
      Addition_Order = paste(fm$predictors_added, collapse = " → "),
      stringsAsFactors = FALSE
    )
    
    # Add cumulative AIC improvements
    if (fm$n_iterations > 0) {
      cumulative_improvements <- sapply(fm$iteration_details, function(x) x$aic_improvement)
      summary_row$Total_AIC_Improvement <- sum(cumulative_improvements)
      summary_row$Final_Formula <- paste(deparse(fm$final_formula), collapse = " ")
    } else {
      summary_row$Total_AIC_Improvement <- 0
      summary_row$Final_Formula <- paste(deparse(fm$base_formula), collapse = " ")
    }
    
    summary_table <- rbind(summary_table, summary_row)
  }
  
  # Save results
  write.csv(summary_table, 
            file.path(output_dir, "iterative_expansion_summary.csv"),
            row.names = FALSE)
  
  saveRDS(list(
    iterations = all_iterations,
    final_models = final_models,
    summary = summary_table
  ), file.path(output_dir, "iterative_expansion_results.rds"))
  
  # Create visualizations
  create_iterative_plots(all_iterations, output_dir)
  
  cat("\n\nResults saved to:", output_dir, "\n")
  
  return(list(
    summary = summary_table,
    final_models = final_models,
    iterations = all_iterations
  ))
}

# Function to create visualizations
create_iterative_plots <- function(all_iterations, output_dir) {
  
  # 1. Create cumulative improvement plot
  plot_data <- data.frame()
  
  for (model_key in names(all_iterations)) {
    iterations <- all_iterations[[model_key]]
    
    if (length(iterations) > 0) {
      # Build cumulative data
      cumulative_improvement <- 0
      predictors <- character(0)
      
      for (i in 1:length(iterations)) {
        iter <- iterations[[i]]
        cumulative_improvement <- cumulative_improvement + iter$aic_improvement
        predictors <- c(predictors, iter$added_predictor)
        
        plot_data <- rbind(plot_data, data.frame(
          Model = model_key,
          Iteration = i,
          Predictor_Added = iter$added_predictor,
          Step_Improvement = iter$aic_improvement,
          Cumulative_Improvement = cumulative_improvement,
          Predictors = paste(predictors, collapse = "\n+ "),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(plot_data) > 0) {
    # Cumulative improvement plot
    p1 <- ggplot(plot_data, aes(x = Iteration, y = Cumulative_Improvement, 
                                 color = Model, group = Model)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      geom_text(aes(label = Predictor_Added), 
                hjust = -0.1, vjust = -0.5, size = 3, angle = 45) +
      scale_x_continuous(breaks = 1:max(plot_data$Iteration)) +
      labs(title = "Cumulative AIC Improvement Through Iterative Addition",
           x = "Iteration", y = "Cumulative AIC Improvement") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "cumulative_improvements.png"),
           p1, width = 12, height = 8, dpi = 300)
    
    # Step-wise improvements
    p2 <- ggplot(plot_data, aes(x = Iteration, y = Step_Improvement, 
                                 fill = Model)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      geom_text(aes(label = round(Step_Improvement, 1)), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3) +
      scale_x_continuous(breaks = 1:max(plot_data$Iteration)) +
      labs(title = "Step-wise AIC Improvements",
           subtitle = "Red line indicates threshold of 2 AIC units",
           x = "Iteration", y = "AIC Improvement") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_dir, "stepwise_improvements.png"),
           p2, width = 12, height = 8, dpi = 300)
  }
  
  # 2. Create predictor frequency plot
  all_predictors <- unlist(lapply(all_iterations, function(x) {
    sapply(x, function(iter) iter$added_predictor)
  }))
  
  if (length(all_predictors) > 0) {
    predictor_freq <- as.data.frame(table(all_predictors))
    names(predictor_freq) <- c("Predictor", "Frequency")
    
    p3 <- ggplot(predictor_freq, aes(x = reorder(Predictor, Frequency), 
                                      y = Frequency)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Frequency of Predictors in Final Models",
           x = "", y = "Number of Models") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "predictor_frequency.png"),
           p3, width = 10, height = 6, dpi = 300)
  }
}

# Run the analysis
results <- run_iterative_stepwise()