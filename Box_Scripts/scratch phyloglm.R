### phyloglm
### 

library(phylolm)


boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

longevitydata = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/justlongevitydata.csv")

dfIn1 = read.csv(newdata)
dfIn = merge(dfIn1, longevitydata, by = "species", all.x = T)
dflongsubset = merge(dfIn1, longevitydata, by = "species")
tree = read.nexus(treefile)
rownames(dfIn) = dfIn$species

complete_vars <- c("TerritorialityWeakVsStrong", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("Territory", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
  
dfIn_clean <- dfIn[complete.cases(dfIn[,complete_vars]),]
dfIn_clean <- dfIn_clean[which(!duplicated(dfIn_clean$species)),]
dfIn_clean$logMass_AVONET = log(dfIn_clean$Mass_AVONET)
dfIn_clean$Territory = as.character(dfIn_clean$Territory)

tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfIn_clean$species))
rownames(dfIn_clean) = dfIn_clean$species

fit0 = phyloglm(FemaleSong_Agg01 ~ 1, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit1 = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit2 = phyloglm(FemaleSong_Agg01 ~ Territory, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit2b = phyloglm(FemaleSong_Agg01 ~ logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit3ni = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit3i = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * Territory, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit1wMassNI = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit1wMassI = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit2wMassNI = phyloglm(FemaleSong_Agg01 ~ Territory + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit2wMassI = phyloglm(FemaleSong_Agg01 ~ Territory * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit3niWmassNI = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fit3iWmassI = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * Territory * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fitMassPlusTerrICoop = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * Territory + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fitCoopPlusTerrIMass = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)
fitTerrPlusCoopIMass = phyloglm(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE, save = FALSE)


AIC_list <- list(
  m0 = AIC(fit0),
  m1 = AIC(fit1), 
  m2 = AIC(fit2),
  m3ni = AIC(fit3ni),
  m3i = AIC(fit3i),
  m1wMassNI = AIC(fit1wMassNI),
  m1wMassI = AIC(fit1wMassI), 
  m2wMassNI = AIC(fit2wMassNI),
  m2wMassI = AIC(fit2wMassI),
  m3niWmassNI = AIC(fit3niWmassNI),
  m3iWmassI = AIC(fit3iWmassI),
  fitMassPlusTerrICoop = AIC(fitMassPlusTerrICoop),
  fitCoopPlusTerrIMass = AIC(fitCoopPlusTerrIMass),
  fitTerrPlusCoopIMass = AIC(fitTerrPlusCoopIMass)
)
model_comparison <- data.frame(
  Model = names(AIC_list),
  AIC = unlist(AIC_list)
)

# Sort by AIC (ascending order)
model_comparison <- model_comparison[order(model_comparison$AIC), ]

# Calculate delta AIC
model_comparison$deltaAIC <- model_comparison$AIC - min(model_comparison$AIC)

# Calculate Akaike weights
model_comparison$weight <- exp(-0.5 * model_comparison$deltaAIC) / 
  sum(exp(-0.5 * model_comparison$deltaAIC))

# Add a column for cumulative weights (useful for determining 95% confidence set)
model_comparison$cum_weight <- cumsum(model_comparison$weight)

# Format the numeric columns to have fewer decimal places
model_comparison$AIC <- round(model_comparison$AIC, 2)
model_comparison$deltaAIC <- round(model_comparison$deltaAIC, 2)
model_comparison$weight <- round(model_comparison$weight, 3)
model_comparison$cum_weight <- round(model_comparison$cum_weight, 3)

# Display the comparison table
print(model_comparison)


## Increase bootstrapping for top few models
fitMassPlusTerrICoop = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)
fitCoopPlusTerrIMass = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)
fitTerrPlusCoopIMass = phyloglm(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop * logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)
fit3niWmassNI = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)
fit3i = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)
fit3ni = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = FALSE)

model_comparison <- data.frame(
  Model = c("fit3ni", "fit3i",
            "fit3niWmassNI", "fitMassPlusTerrICoop", "fitCoopPlusTerrIMass", "fitTerrPlusCoopIMass" ),
  AIC = c(AIC(fit3ni), AIC(fit3i),
          AIC(fit3niWmassNI), AIC(fitMassPlusTerrICoop), AIC(fitCoopPlusTerrIMass), AIC(fitTerrPlusCoopIMass))
)

# Sort by AIC (ascending order)
model_comparison <- model_comparison[order(model_comparison$AIC), ]

# Calculate delta AIC
model_comparison$deltaAIC <- model_comparison$AIC - min(model_comparison$AIC)

# Calculate Akaike weights
model_comparison$weight <- exp(-0.5 * model_comparison$deltaAIC) / 
  sum(exp(-0.5 * model_comparison$deltaAIC))

# Add a column for cumulative weights (useful for determining 95% confidence set)
model_comparison$cum_weight <- cumsum(model_comparison$weight)

# Format the numeric columns to have fewer decimal places
model_comparison$AIC <- round(model_comparison$AIC, 2)
model_comparison$deltaAIC <- round(model_comparison$deltaAIC, 2)
model_comparison$weight <- round(model_comparison$weight, 3)
model_comparison$cum_weight <- round(model_comparison$cum_weight, 3)

# Display the comparison table
print(model_comparison)

summary(fitMassPlusTerrICoop)
summary(fit3niWmassNI)
summary(fitCoopPlusTerrIMass)


# Calculate Akaike weights
model_comparison <- data.frame(
  Model = c("fitMassPlusTerrICoop", "fit3niWmassNI", "fitCoopPlusTerrIMass"),
  AIC = c(AIC(fitMassPlusTerrICoop), AIC(fit3niWmassNI), AIC(fitCoopPlusTerrIMass))
)
model_comparison$deltaAIC <- model_comparison$AIC - min(model_comparison$AIC)
model_comparison$weight <- exp(-0.5 * model_comparison$deltaAIC) / 
  sum(exp(-0.5 * model_comparison$deltaAIC))

# Extract coefficients from each model
coef_model1 <- coef(fitMassPlusTerrICoop)
coef_model2 <- coef(fit3niWmassNI)
coef_model3 <- coef(fitCoopPlusTerrIMass)

# Create a common set of parameter names
all_params <- unique(c(names(coef_model1), names(coef_model2), names(coef_model3)))

# Create matrix to store coefficients
coef_matrix <- matrix(0, nrow = length(all_params), ncol = 3)
rownames(coef_matrix) <- all_params
colnames(coef_matrix) <- c("Model1", "Model2", "Model3")

# Fill in the coefficient matrix (parameters not in a model remain 0)
for (param in intersect(all_params, names(coef_model1)))
  coef_matrix[param, "Model1"] <- coef_model1[param]
for (param in intersect(all_params, names(coef_model2)))
  coef_matrix[param, "Model2"] <- coef_model2[param]
for (param in intersect(all_params, names(coef_model3)))
  coef_matrix[param, "Model3"] <- coef_model3[param]

# Calculate weighted average coefficients
avg_coefs <- coef_matrix %*% model_comparison$weight
names(avg_coefs) <- all_params

# Create a data frame with the results
avg_coef_df <- data.frame(
  Parameter = all_params,
  Estimate = avg_coefs,
  Model1 = coef_matrix[, "Model1"],
  Model2 = coef_matrix[, "Model2"],
  Model3 = coef_matrix[, "Model3"],
  Weight1 = model_comparison$weight[1],
  Weight2 = model_comparison$weight[2],
  Weight3 = model_comparison$weight[3]
)

# Print the results
print(avg_coef_df)

#### plot averaged models ----
# Create a grid of predictors
new_data <- expand.grid(
  HighConfidence_Coop = c(0, 1),
  TerritorialityWeakVsStrong = c(0, 1),
  logMass_AVONET = mean(dfIn_clean$logMass_AVONET)
)

# Calculate predictions from each model
pred1 <- get_model_predictions(fitMassPlusTerrICoop, new_data)
pred2 <- get_model_predictions(fit3niWmassNI, new_data)
pred3 <- get_model_predictions(fitCoopPlusTerrIMass, new_data)

# Weight the predictions
weighted_pred <- pred1 * model_comparison$weight[1] + 
  pred2 * model_comparison$weight[2] + 
  pred3 * model_comparison$weight[3]

# Add to the data frame
new_data$avg_pred <- weighted_pred

# Create the plot
ggplot(new_data, aes(x = factor(HighConfidence_Coop), y = avg_pred, 
                     group = factor(TerritorialityWeakVsStrong), 
                     color = factor(TerritorialityWeakVsStrong))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(x = "Cooperative Breeding", 
       y = "Model-Averaged Probability of Female Song",
       color = "Territoriality") +
  scale_color_manual(values = c("darkblue", "darkred"),
                     labels = c("Weak", "Strong")) +
  scale_x_discrete(labels = c("Absent", "Present")) +
  theme_bw() +
  ggtitle("Model-Averaged Effect of Cooperative Breeding and Territoriality") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Helper function for model predictions
get_model_predictions <- function(model, newdata) {
  # Extract model coefficients
  coefs <- coef(model)
  
  # Calculate linear predictor based on model formula
  if (identical(model, fitMassPlusTerrICoop)) {
    lp <- coefs[1] + 
      coefs[2] * newdata$HighConfidence_Coop + 
      coefs[3] * newdata$TerritorialityWeakVsStrong + 
      coefs[4] * newdata$logMass_AVONET + 
      coefs[5] * newdata$HighConfidence_Coop * newdata$TerritorialityWeakVsStrong
  } else if (identical(model, fit3niWmassNI)) {
    lp <- coefs[1] + 
      coefs[2] * newdata$HighConfidence_Coop + 
      coefs[3] * newdata$TerritorialityWeakVsStrong + 
      coefs[4] * newdata$logMass_AVONET
  } else if (identical(model, fitCoopPlusTerrIMass)) {
    lp <- coefs[1] + 
      coefs[2] * newdata$HighConfidence_Coop + 
      coefs[3] * newdata$TerritorialityWeakVsStrong + 
      coefs[4] * newdata$logMass_AVONET + 
      coefs[5] * newdata$TerritorialityWeakVsStrong * newdata$logMass_AVONET
  }
  
  # Convert to probabilities
  return(plogis(lp))
}


fitMassPlusTerrICoop_with_boot = phyloglm(FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_AVONET, data = dfIn_clean, phy = tree_clean, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL, boot = 100, full.matrix = TRUE, save = TRUE)

## Visualize interaction
# Set up a grid of values for prediction
new_data <- expand.grid(
  HighConfidence_Coop = c(0, 1),
  TerritorialityWeakVsStrong = c(0, 1),
  logMass_AVONET = mean(dfIn_clean$logMass_AVONET)  # Hold mass at its mean
)

# Extract the original coefficients
beta0 <- coef(fitMassPlusTerrICoop_with_boot)[1]  # Intercept
beta1 <- coef(fitMassPlusTerrICoop_with_boot)[2]  # HighConfidence_Coop
beta2 <- coef(fitMassPlusTerrICoop_with_boot)[3]  # TerritorialityWeakVsStrong
beta3 <- coef(fitMassPlusTerrICoop_with_boot)[4]  # logMass_AVONET
beta4 <- coef(fitMassPlusTerrICoop_with_boot)[5]  # Interaction term

# Calculate predictions for the original model
new_data$linear_pred <- beta0 + 
  beta1 * new_data$HighConfidence_Coop + 
  beta2 * new_data$TerritorialityWeakVsStrong + 
  beta3 * new_data$logMass_AVONET + 
  beta4 * new_data$HighConfidence_Coop * new_data$TerritorialityWeakVsStrong
new_data$pred_prob <- plogis(new_data$linear_pred)

# Initialize matrices to store bootstrap predictions
n_boot <- nrow(fitMassPlusTerrICoop_with_boot$bootstrap)
boot_preds <- matrix(NA, nrow = nrow(new_data), ncol = n_boot)

# Access bootstrap samples - these should be in the 'bootstrap' element of your model
if (exists("bootstrap", where = fitMassPlusTerrICoop_with_boot)) {
  
  # Assuming coefficients are in columns of the bootstrap matrix
  boot_coefs <- fitMassPlusTerrICoop_with_boot$bootstrap[i, ]
  
  # For each bootstrap sample
  for (i in 1:n_boot) {
    # Extract coefficients from this bootstrap sample
    # Assuming coefficients are in columns of the bootstrap matrix
    boot_coefs <- fitMassPlusTerrICoop_with_boot$bootstrap[i, ]
    
    # Calculate linear predictor for this sample
    boot_linear_pred <- boot_coefs[1] + 
      boot_coefs[2] * new_data$HighConfidence_Coop + 
      boot_coefs[3] * new_data$TerritorialityWeakVsStrong + 
      boot_coefs[4] * new_data$logMass_AVONET + 
      boot_coefs[5] * new_data$HighConfidence_Coop * new_data$TerritorialityWeakVsStrong
    
    # Convert to probabilities
    boot_preds[, i] <- plogis(boot_linear_pred)
  }
  
  # Calculate quantiles for confidence intervals (2.5% and 97.5% for 95% CI)
  new_data$lower_ci <- apply(boot_preds, 1, quantile, probs = 0.025)
  new_data$upper_ci <- apply(boot_preds, 1, quantile, probs = 0.975)
} else {
  # If bootstrap samples aren't available, create a warning
  warning("Bootstrap samples not available. Run model with save=TRUE to get confidence intervals.")
  # Set dummy CIs equal to the point estimate (will be invisible in plot)
  new_data$lower_ci <- new_data$pred_prob
  new_data$upper_ci <- new_data$pred_prob
}

# Create plot with confidence intervals
library(ggplot2)

#Figure: Predicted probability of female song occurrence in relation to cooperative breeding status, moderated by territoriality level. This figure illustrates the significant interaction between cooperative breeding and territoriality on the likelihood of female song in birds, while controlling for body mass. Lines represent model-predicted probabilities derived from the best-fitting phylogenetic logistic regression model. Points represent predicted values at mean body mass. Shaded regions indicate 95% confidence intervals based on 100 parametric bootstrap replicates.
ggplot(new_data, aes(x = factor(HighConfidence_Coop), y = pred_prob, 
                     group = factor(TerritorialityWeakVsStrong), 
                     color = factor(TerritorialityWeakVsStrong),
                     fill = factor(TerritorialityWeakVsStrong))) +
  # Add confidence interval ribbons
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  # Add lines and points
  geom_line(size = 1) +
  geom_point(size = 3) +
  # Labels and styling
  labs(x = "Cooperative Breeding", 
       y = "Probability of Female Song",
       color = "Territoriality",
       fill = "Territoriality") +
  scale_color_manual(values = c("darkblue", "darkred"),
                     labels = c("Weak", "Strong")) +
  scale_fill_manual(values = c("darkblue", "darkred"),
                    labels = c("Weak", "Strong")) +
  scale_x_discrete(labels = c("Absent", "Present")) +
  ylim(0, 1) +  # Ensure y-axis covers full probability range
  theme_bw() +
  ggtitle("Interaction between Cooperative Breeding and Territoriality") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")
