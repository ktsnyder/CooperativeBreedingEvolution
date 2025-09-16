## Test phyloglm output plots
## 

## Forest plots for coefficient comparison ----
library(ggplot2)
library(dplyr)

# Extract coefficients from multiple models
coef_data <- data.frame()
for (key in names(all_model_results)) {
  # Get best model for each analysis
  best_model_name <- all_model_results[[key]]$comparison$Model[1]
  best_model <- all_model_results[[key]]$models[[best_model_name]]
  
  # Extract model summary
  model_summary <- summary(best_model)
  coefs <- model_summary$coefficients
  
  # Create data frame for this model
  model_df <- data.frame(
    Model = key,
    Parameter = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    StdErr = coefs[, "StdErr"],
    LowerCI = coefs[, "Estimate"] - 1.96 * coefs[, "StdErr"],
    UpperCI = coefs[, "Estimate"] + 1.96 * coefs[, "StdErr"],
    pvalue = coefs[, "p.value"]
  )
  
  coef_data <- rbind(coef_data, model_df)
}

# Create forest plot for a specific parameter
ggplot(subset(coef_data, Parameter == "HighConfidence_Coop"), 
       aes(x = Model, y = Estimate)) +
  geom_point(aes(size = -log10(pvalue))) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(title = "Effect of Cooperative Breeding Across Models",
       y = "Coefficient Estimate (log-odds scale)",
       x = "Territoriality Variable")


## Faceted interaction plots ----
library(ggplot2)
library(dplyr)

# Function to manually calculate predictions for multiple models
create_prediction_grid_manual <- function(model_results_list) {
  # Create empty data frame to store all predictions
  all_predictions <- data.frame()
  
  # Loop through each analysis
  for (key in names(model_results_list)) {
    # Skip if models is NULL or empty
    if (is.null(model_results_list[[key]]$models) || length(model_results_list[[key]]$models) == 0) {
      message(paste("Skipping", key, "- no models available"))
      next
    }
    
    # Get parts from key
    parts <- strsplit(key, "_")[[1]]
    resp_var <- parts[1]
    terr_var <- paste(parts[-1], collapse = "_")
    
    # Get best model
    best_model_name <- model_results_list[[key]]$comparison$Model[1]
    best_model <- model_results_list[[key]]$models[[best_model_name]]
    
    # Get data
    data <- model_results_list[[key]]$data
    
    # Mean mass for predictions
    mean_mass <- mean(data$logMass_AVONET, na.rm = TRUE)
    
    # Create prediction grid
    grid <- expand.grid(
      HighConfidence_Coop = c(0, 1),
      Territoriality = c(0, 1),
      logMass_AVONET = mean_mass
    )
    
    # Create a cleaner model key for display
    clean_key <- gsub("_", " + ", key)
    clean_key <- gsub("FemaleSong_Agg01", "Female Song", clean_key)
    clean_key <- gsub("HighConfidence_Coop", "Cooperative Breeding", clean_key)
    clean_key <- gsub("Territoriality", "Territoriality: ", clean_key)
    
    # Add model info
    grid$ModelKey <- key
    grid$Model <- clean_key
    grid$Response <- resp_var
    grid$TerrVar <- terr_var
    
    # Rename territoriality column to match model
    names(grid)[names(grid) == "Territoriality"] <- "TerrValue"
    
    # Manual calculation using coefficients
    coefs <- coef(best_model)
    
    # Start with intercept
    grid$linear_pred <- coefs["(Intercept)"]
    
    # Add effect of cooperative breeding if present
    if ("HighConfidence_Coop" %in% names(coefs)) {
      grid$linear_pred <- grid$linear_pred + coefs["HighConfidence_Coop"] * grid$HighConfidence_Coop
    }
    
    # Add effect of territoriality if present
    if (terr_var %in% names(coefs)) {
      grid$linear_pred <- grid$linear_pred + coefs[terr_var] * grid$TerrValue
    }
    
    # Add effect of log mass if present
    if ("logMass_AVONET" %in% names(coefs)) {
      grid$linear_pred <- grid$linear_pred + coefs["logMass_AVONET"] * grid$logMass_AVONET
    }
    
    # Add interaction if present
    interaction_term <- paste0("HighConfidence_Coop:", terr_var)
    if (interaction_term %in% names(coefs)) {
      grid$linear_pred <- grid$linear_pred + 
        coefs[interaction_term] * grid$HighConfidence_Coop * grid$TerrValue
    }
    
    # Calculate predicted probabilities
    grid$pred_prob <- plogis(grid$linear_pred)
    
    # Add to overall predictions
    all_predictions <- rbind(all_predictions, grid)
  }
  
  return(all_predictions)
}

# Create predictions
all_preds <- create_prediction_grid_manual(all_model_results)

# Create fancy labels
all_preds$CoopLabel <- ifelse(all_preds$HighConfidence_Coop == 0, 
                              "Cooperative Breeding: Absent", 
                              "Cooperative Breeding: Present")
all_preds$TerrLabel <- ifelse(all_preds$TerrValue == 0, 
                              "Territoriality: 0", 
                              "Territoriality: 1")

# Create a combined scenario label
all_preds$Scenario <- paste(all_preds$CoopLabel, all_preds$TerrLabel, sep = ", ")

# Filter to just female song predictions for plotting
female_song_preds <- all_preds %>%
  filter(Response == "FemaleSong_Agg01")

# Create faceted plot by territoriality variable
ggplot(female_song_preds, 
       aes(x = factor(HighConfidence_Coop, labels = c("Absent", "Present")), 
           y = pred_prob, 
           color = factor(TerrValue, labels = c("0", "1")), 
           group = factor(TerrValue))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ TerrVar, scales = "free") +
  theme_bw() +
  labs(x = "Cooperative Breeding",
       y = "Probability of Female Song",
       color = "Territoriality Value",
       title = "Effects of Cooperative Breeding Across Territoriality Variables") +
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 8))



## Coefficient heat map ----
library(ggplot2)
library(dplyr)
library(reshape2)

# Extract coefficients from multiple models
coef_data <- data.frame()
for (key in names(all_model_results)) {
  # Skip if models is NULL or empty
  if (is.null(all_model_results[[key]]$models) || length(all_model_results[[key]]$models) == 0) {
    next
  }
  
  # Get best model for each analysis
  best_model_name <- all_model_results[[key]]$comparison$Model[1]
  best_model <- all_model_results[[key]]$models[[best_model_name]]
  
  # Extract model coefficients
  model_coefs <- coef(best_model)
  model_summary <- summary(best_model)
  
  # Get p-values
  p_values <- model_summary$coefficients[, "p.value"]
  
  # For each coefficient in this model
  for (coef_name in names(model_coefs)) {
    # Add a row to the data frame
    coef_data <- rbind(coef_data, data.frame(
      Model = key,
      Parameter = coef_name,
      Estimate = model_coefs[coef_name],
      p_value = p_values[coef_name],
      stringsAsFactors = FALSE
    ))
  }
}

# Add significance indicators
coef_data$Significance <- ifelse(coef_data$p_value < 0.001, "***",
                                 ifelse(coef_data$p_value < 0.01, "**",
                                        ifelse(coef_data$p_value < 0.05, "*",
                                               ifelse(coef_data$p_value < 0.1, ".", ""))))

# Create a plot-ready data frame
plot_data <- coef_data

# Create a simplified heatmap
ggplot(plot_data, aes(x = Model, y = Parameter, fill = Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Coefficient") +
  geom_text(aes(label = Significance), color = "black", size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coefficients and Significance Across Models",
       x = "Model", y = "Parameter")



## Predicted probability plots ----
# Reshape data for plotting
pred_comparison <- all_preds %>%
  filter(Response == "FemaleSong_Agg01") %>%
  mutate(Scenario = paste0(
    "Coop: ", ifelse(HighConfidence_Coop == 0, "Absent", "Present"),
    ", Terr: ", ifelse(get(Territoriality) == 0, "0", "1")
  ))

# Create comparison plot
ggplot(pred_comparison, 
       aes(x = Territoriality, y = predicted_prob, 
           fill = Scenario)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = round(predicted_prob, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predicted Probabilities Across Territoriality Variables",
       x = "Territoriality Variable", 
       y = "Probability of Female Song")