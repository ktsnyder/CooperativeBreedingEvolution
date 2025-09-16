# PhyloPath scratch

library(phylopath)
library(phytools)


boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

#setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

source("subsettreedata.R")

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

tree = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

complete_vars <- c("Territory", "FemaleSong_Agg01", "HighConfidence_Coop")
dfIn_clean <- dfIn[complete.cases(dfIn[,complete_vars]),]

#dfIn_clean$HighConfidence_Coop = as.factor(dfIn_clean$HighConfidence_Coop)
#dfIn_clean$FemaleSong_Agg01 = as.factor(dfIn_clean$FemaleSong_Agg01)
dfIn_clean$Territory <- as.character(dfIn_clean$Territory)
dfIn_clean$Territory <- factor(dfIn_clean$Territory, levels = c("1", "2", "3"))
unique(dfIn_clean$HighConfidence_Coop)
unique(dfIn_clean$FemaleSong_Agg01)
unique(dfIn_clean$Territory)

# Create split binary variables for Territory
Territory.matrix <- model.matrix(~ Territory - 1, data = dfIn_clean)
dfIn_clean$Territory1 <- Territory.matrix[,1]
dfIn_clean$Territory2 <- Territory.matrix[,2]
dfIn_clean$Territory3 <- Territory.matrix[,3]

dfIn_clean$HighConfidence_Coop <- as.factor(dfIn_clean$HighConfidence_Coop)
dfIn_clean$Territory1 <- as.factor(dfIn_clean$Territory1)
dfIn_clean$Territory2 <- as.factor(dfIn_clean$Territory2)
dfIn_clean$Territory3 <- as.factor(dfIn_clean$Territory3)
dfIn_clean$FemaleSong_Agg01 <- as.factor(dfIn_clean$FemaleSong_Agg01)

tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfIn_clean$species))
rownames(dfIn_clean) <- dfIn_clean$species

# models <- define_model_set(
#   m0 = c(FemaleSong_Agg01 ~ Territory),
#   m1 = c(FemaleSong_Agg01 ~ Territory, FemaleSong_Agg01 ~ HighConfidence_Coop), 
#   m2 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
#   m3 = c(FemaleSong_Agg01 ~ Territory, FemaleSong_Agg01 ~ HighConfidence_Coop, HighConfidence_Coop ~ Territory), 
#   m4 = c(FemaleSong_Agg01 ~ Territory, HighConfidence_Coop ~ Territory), 
#   m5 = c(FemaleSong_Agg01 ~ HighConfidence_Coop, HighConfidence_Coop ~ Territory),
#   m6 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
#   m7 = c(HighConfidence_Coop ~ FemaleSong_Agg01, FemaleSong_Agg01 ~ Territory),
#   .common = NULL
# )
# Define models using dummy variables
models <- define_model_set(
  m0 = c(FemaleSong_Agg01 ~ Territory2 + Territory3),
  m0b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory2),
  m0c = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3),
  m0d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory2),
  m1b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory2,
          HighConfidence_Coop ~ Territory2 + Territory3),
  m1c = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3,
          HighConfidence_Coop ~ Territory2 + Territory3),
  m1 = c(FemaleSong_Agg01 ~ Territory2 + Territory3 + HighConfidence_Coop, 
         HighConfidence_Coop ~ Territory2 + Territory3),
  m2 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ Territory2 + Territory3),
  m3 = c(FemaleSong_Agg01 ~ Territory2 + Territory3,
         HighConfidence_Coop ~ Territory2 + Territory3),
  m4 = c(FemaleSong_Agg01 ~ Territory2 + Territory3,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  m5 = c(FemaleSong_Agg01 ~ Territory2 + Territory3,
         HighConfidence_Coop ~ FemaleSong_Agg01 + Territory2 + Territory3),
  m6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory2,
         HighConfidence_Coop ~ Territory2 + Territory3),
  m7 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3,
         HighConfidence_Coop ~ Territory2 + Territory3),
  m8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory2,
         HighConfidence_Coop ~ Territory2),
  m9 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory2,
         HighConfidence_Coop ~ Territory3)
)

models$m0
plot(models$m0)
plot(models$m1)
plot(models$m2)
plot(models$m4)
plot(models$m5)
plot(models$m6)
plot(models$m3)

pdf(paste0("plot phylopath models_",Sys.Date(),".pdf"), height = 15, width = 15)
plot_model_set(models)
dev.off()


result <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'lambda',
  na.rm = TRUE
)
result

(s <- summary(result))

plot(s)

(best_model <- best(result))

plot(best_model)




#### Testing alternative: phylolm ----

tree = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

complete_vars <- c("FemaleSong_Agg01", "HighConfidence_Coop", "Territory")
df <- dfIn[complete.cases(dfIn[,complete_vars]),]

# Check tree and data match
df <- df[which(df$species %in% tree$tip.label),]
tree <- drop.tip(tree, setdiff(tree$tip.label, df$species))

# Convert Territory to numeric contrasts
df$Territory <- factor(df$Territory, levels = c("1", "2", "3"))
Territory.matrix <- model.matrix(~ Territory - 1, data = df)
df$Territory2 <- Territory.matrix[,2]
df$Territory3 <- Territory.matrix[,3]

# Convert response variables to numeric
df$HighConfidence_Coop <- as.numeric(df$HighConfidence_Coop) - 1
df$FemaleSong_Agg01 <- as.numeric(df$FemaleSong_Agg01) - 1

# Territory -> HighConfidence_Coop
m1 <- phylolm(HighConfidence_Coop ~ Territory, phy = tree, model = "lambda", data = df)

# Territory + HighConfidence_Coop -> FemaleSong_Agg01  
m2 <- phylolm(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop, phy = tree, model = "lambda", data = df)

# Territory -> FemaleSong_Agg01
m3 <- phylolm(FemaleSong_Agg01 ~ Territory, phy = tree, model = "lambda", data = df)

# HighConfidence_Coop -> FemaleSong_Agg01
m4 <- phylolm(FemaleSong_Agg01 ~ HighConfidence_Coop, phy = tree, model = "lambda", data = df)



#### phyloglm ----
library(phylolm)

# Convert factors to numeric (0/1) for binary responses
df$HighConfidence_Coop_num = as.numeric(df$HighConfidence_Coop) - 1
df$FemaleSong_num = as.numeric(df$FemaleSong_Agg01) - 1

# # Territory -> HighConfidence_Coop
# m1 <- phyloglm(HighConfidence_Coop_num ~ Territory, phy = tree, method = "logistic_MPLE", data = df)
# 
# # Territory + HighConfidence_Coop -> FemaleSong_Agg01  
# m2 <- phyloglm(FemaleSong_num ~ Territory + HighConfidence_Coop, phy = tree, method = "logistic_MPLE", data = df)
# 
# # Territory -> FemaleSong_Agg01
# m3 <- phyloglm(FemaleSong_num ~ Territory, phy = tree, method = "logistic_MPLE", data = df)
# 
# # HighConfidence_Coop -> FemaleSong_Agg01
# m4 <- phyloglm(FemaleSong_num ~ HighConfidence_Coop, phy = tree, method = "logistic_MPLE", data = df)


# Territory -> Cooperative Breeding
m1 <- phylolm(HighConfidence_Coop ~ Territory2 + Territory3, 
              phy = tree, 
              model = "lambda", 
              data = df)

# Territory + Cooperative Breeding -> Female Song
m2 <- phylolm(FemaleSong_Agg01 ~ Territory2 + Territory3 + HighConfidence_Coop, 
              phy = tree,
              model = "lambda", 
              data = df)

# Territory -> Female Song
m3 <- phylolm(FemaleSong_Agg01 ~ Territory2 + Territory3,
              phy = tree,
              model = "lambda",
              data = df)

# Cooperative Breeding -> Female Song  
m4 <- phylolm(FemaleSong_Agg01 ~ HighConfidence_Coop,
              phy = tree,
              model = "lambda",
              data = df)

# Territory2 + Cooperative Breeding -> Female Song
m5 <- phylolm(FemaleSong_Agg01 ~ Territory2 + HighConfidence_Coop, 
              phy = tree,
              model = "lambda", 
              data = df)

phylolm(FemaleSong_Agg01 ~ Territory2 , 
        phy = tree,
        model = "lambda", 
        data = df)

# Territory2 + Cooperative Breeding -> Female Song
m11 <- phylolm(FemaleSong_Agg01 ~ Territory3 + HighConfidence_Coop, 
              phy = tree,
              model = "lambda", 
              data = df)

# Effects on CB
m6 <- phylolm(HighConfidence_Coop ~ Territory2 + Territory3 + FemaleSong_Agg01, 
              phy = tree,
              model = "lambda", 
              data = df)

m7 <- phylolm(HighConfidence_Coop ~ Territory2 + FemaleSong_Agg01, 
              phy = tree,
              model = "lambda", 
              data = df)

phylolm(HighConfidence_Coop ~ Territory2 , 
        phy = tree,
        model = "lambda", 
        data = df)

m8 <- phylolm(HighConfidence_Coop ~ Territory3, 
              phy = tree,
              model = "lambda", 
              data = df)

m9 <- phylolm(HighConfidence_Coop ~ FemaleSong_Agg01, 
              phy = tree,
              model = "lambda", 
              data = df)

m10 <- phylolm(HighConfidence_Coop ~ Territory3 + FemaleSong_Agg01, 
              phy = tree,
              model = "lambda", 
              data = df)


# Compare models
AIC_list <- list(
  m1 = AIC(m1),
  m2 = AIC(m2), 
  m3 = AIC(m3),
  m4 = AIC(m4),
  m5 = AIC(m5),
  m6 = AIC(m6),
  m7 = AIC(m7), 
  m8 = AIC(m8),
  m9 = AIC(m9),
  m10 = AIC(m10)
)

# Compare models using AIC
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
AICs <- sapply(models, AIC)

models
