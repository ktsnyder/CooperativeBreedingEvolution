# Debug interaction plot issue

library(ggplot2)
library(dplyr)
library(patchwork)

# Load data
batch_results <- readRDS('./Outputs/PhyloglmResults/PhyloGLM_Batch_20250710_235142_boot100/all_results.rds')

# Get a specific result
result <- batch_results[["FS_CB_TerrWS_Mass"]]

# Find interaction model
int_models <- result$comparison$comparison %>%
  filter(grepl("x", Model)) %>%
  arrange(AIC)

cat("Interaction models found:\n")
print(int_models)

if (nrow(int_models) > 0) {
  best_model_name <- int_models$Model[1]
  cat("\nBest interaction model:", best_model_name, "\n")
  
  model <- result$models$models[[best_model_name]]
  
  # Check model formula
  cat("\nModel formula:\n")
  print(formula(model))
  
  # Create prediction data
  pred_data <- expand.grid(
    HighConfidence_Coop = c(0, 1),
    TerritorialityWeakVsStrong = c(0, 1)
  )
  
  cat("\nPrediction data:\n")
  print(pred_data)
  
  # Add control variable
  control_var <- result$config$controls[1]
  cat("\nControl variable:", control_var, "\n")
  
  if (!is.null(control_var) && control_var %in% names(result$prepared_data$data)) {
    pred_data[[control_var]] <- mean(result$prepared_data$data[[control_var]], na.rm = TRUE)
    cat("Added control variable with mean:", pred_data[[control_var]][1], "\n")
  }
  
  # Try to create model matrix
  cat("\nTrying to create model matrix...\n")
  tryCatch({
    X <- model.matrix(formula(model), data = pred_data)
    cat("Model matrix created successfully!\n")
    cat("Dimensions:", dim(X), "\n")
    
    # Get predictions
    coefs <- coef(model)
    pred_data$predicted <- plogis(X %*% coefs)
    
    cat("\nPredictions:\n")
    print(pred_data)
    
  }, error = function(e) {
    cat("Error creating model matrix:", e$message, "\n")
    
    # Try alternative approach
    cat("\nTrying alternative approach...\n")
    
    # Get the formula terms
    form <- formula(model)
    cat("Formula terms:", as.character(form), "\n")
    
    # Check what variables the model expects
    cat("\nModel terms:\n")
    print(terms(model))
    
    # Look at the original data used
    cat("\nVariables in original data:\n")
    print(names(model$data)[1:10])
  })
}