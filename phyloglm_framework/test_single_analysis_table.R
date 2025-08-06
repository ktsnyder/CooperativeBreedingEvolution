# Test comprehensive table with just one analysis to debug

source("phyloglm_framework/stepwise_expansion_comprehensive_table.R")

# Override the function to just do one analysis
run_single_test <- function() {
  output_dir <- "Outputs/PhyloglmResults/test_comprehensive_table"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  cat("Loading data...\n")
  all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")
  full_data <- read.csv("Data_R_2025-06-09.csv")
  tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  
  # Just test FS_vs_CB_TerrWS_Mass
  analysis_name <- "FS_vs_CB_TerrWS_Mass"
  
  # Get base model
  base_result <- all_results[[analysis_name]]
  best_model_name <- base_result$comparison$comparison$Model[1]
  best_model <- base_result$models$models[[best_model_name]]
  original_base_aic <- base_result$comparison$comparison$AIC[1]
  base_data <- base_result$prepared_data$data
  base_tree <- base_result$prepared_data$tree
  base_formula <- formula(best_model)
  
  cat("Testing Territory predictor only...\n")
  
  # Just test Territory
  test_data <- base_data
  test_data$Territory_num <- as.numeric(full_data$Territory[match(test_data$species, full_data$species)])
  test_data <- test_data[!is.na(test_data$Territory_num), ]
  test_tree <- keep.tip(base_tree, test_data$species)
  test_data <- test_data[match(test_tree$tip.label, test_data$species), ]
  
  # Refit base model
  base_fit <- phyloglm(
    formula = base_formula,
    data = test_data,
    phy = test_tree,
    method = "logistic_MPLE",
    btol = 50,
    log.alpha.bound = 4
  )
  
  # Fit expanded model
  new_formula <- update(base_formula, ~ . + Territory_num)
  expanded_fit <- phyloglm(
    formula = new_formula,
    data = test_data,
    phy = test_tree,
    method = "logistic_MPLE",
    btol = 50,
    log.alpha.bound = 4
  )
  
  # Create simple output
  output <- data.frame(
    Analysis = analysis_name,
    Base_Model = best_model_name,
    Base_Formula = deparse(base_formula),
    N_Species = nrow(test_data),
    Base_AIC = -2 * base_fit$logLik + 2 * base_fit$d,
    Expanded_AIC = -2 * expanded_fit$logLik + 2 * expanded_fit$d,
    AIC_Improvement = (-2 * base_fit$logLik + 2 * base_fit$d) - (-2 * expanded_fit$logLik + 2 * expanded_fit$d),
    Territory_Coef = coef(expanded_fit)["Territory_num"],
    Territory_p = summary(expanded_fit)$coefficients["Territory_num", "p.value"]
  )
  
  write.csv(output, file.path(output_dir, "test_output.csv"), row.names = FALSE)
  
  cat("\nTest completed. Output saved to:", output_dir, "\n")
  print(output)
}

# Run the test
run_single_test()