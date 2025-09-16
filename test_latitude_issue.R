library(phylolm)
library(ape)

# Load data
data <- read.csv("Data_R_2025-06-09.csv")
tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Prepare a simple test case
test_data <- data[c("species", "HighConfidence_Coop", "FemaleSong_Agg01")]
na_mask <- is.na(test_data$HighConfidence_Coop) | is.na(test_data$FemaleSong_Agg01)
test_data <- test_data[\!na_mask,]
test_data$abs_Latitude <- abs(data$Centroid.Latitude_AVONET[match(test_data$species, data$species)])
test_data <- test_data[\!is.na(test_data$abs_Latitude),]

# Match tree
test_tree <- keep.tip(tree, test_data$species)
test_data <- test_data[match(test_tree$tip.label, test_data$species),]

cat("N species:", nrow(test_data), "\n")

# Fit base model
base_fit <- phyloglm(HighConfidence_Coop ~ FemaleSong_Agg01, 
                     data = test_data, phy = test_tree, 
                     method = "logistic_MPLE", btol = 50)
cat("Base logLik:", base_fit$logLik, "\n")
cat("Base AIC:", -2 * base_fit$logLik + 2 * base_fit$d, "\n")

# Fit model with latitude
lat_fit <- phyloglm(HighConfidence_Coop ~ FemaleSong_Agg01 + abs_Latitude, 
                    data = test_data, phy = test_tree, 
                    method = "logistic_MPLE", btol = 50)
cat("\nWith latitude logLik:", lat_fit$logLik, "\n")
cat("With latitude AIC:", -2 * lat_fit$logLik + 2 * lat_fit$d, "\n")
cat("LogLik difference:", lat_fit$logLik - base_fit$logLik, "\n")
cat("AIC improvement:", (-2 * base_fit$logLik + 2 * base_fit$d) - (-2 * lat_fit$logLik + 2 * lat_fit$d), "\n")

# Check coefficient
coef_info <- summary(lat_fit)$coefficients["abs_Latitude",]
cat("\nLatitude coefficient:", coef_info["Estimate"], "\n")
cat("Latitude p-value:", coef_info["p.value"], "\n")

# Check optimization details
cat("\nBase alpha:", base_fit$alpha, "\n")
cat("Latitude alpha:", lat_fit$alpha, "\n")
cat("Base convergence:", base_fit$convergence, "\n")
cat("Latitude convergence:", lat_fit$convergence, "\n")

# Test with standard glm for comparison
cat("\n--- Standard GLM comparison ---\n")
base_glm <- glm(HighConfidence_Coop ~ FemaleSong_Agg01, 
                data = test_data, family = binomial)
lat_glm <- glm(HighConfidence_Coop ~ FemaleSong_Agg01 + abs_Latitude, 
               data = test_data, family = binomial)
cat("Base GLM logLik:", logLik(base_glm), "\n")
cat("Latitude GLM logLik:", logLik(lat_glm), "\n")
cat("GLM LogLik difference:", as.numeric(logLik(lat_glm) - logLik(base_glm)), "\n")
