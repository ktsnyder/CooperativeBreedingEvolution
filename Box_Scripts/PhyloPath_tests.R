# PhyloPath scratch

library(phylopath)
library(phytools)


boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')
newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

longevitydata = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/justlongevitydata.csv")

dfIn1 = read.csv(newdata)
dfIn = merge(dfIn1, longevitydata, by = "species", all.x = T)
dflongsubset = merge(dfIn1, longevitydata, by = "species")
tree = read.nexus(treefile)
rownames(dfIn) = dfIn$species

complete_vars <- c("Territory", "FemaleSong_Agg01", "HighConfidence_Coop", "WeakStrongTerr.Rough")
complete_vars <- c("TerritorialityWeakVsStrong", "FemaleSong_Agg01", "HighConfidence_Coop", "log_Longevity")
complete_vars <- c("TerritorialityWeakVsStrong", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("TerritorialityWeakVsStrongHighConf", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("TerritorialityPermissiveExclusive", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("TerritorialityPermissiveColonialCoopVsExclusive", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("TerritorialityPermissiveExclusiveHighConf", "FemaleSong_Agg01", "HighConfidence_Coop", "Mass_AVONET")
complete_vars <- c("Territory_12vs3", "FemaleSong_Agg01", "HighConfidence_Coop", "log_Longevity")
dfIn_clean <- dfIn[complete.cases(dfIn[,complete_vars]),]
dfIn_clean <- dfIn_clean[which(!duplicated(dfIn_clean$species)),]
dfIn_clean$logMass_AVONET = log(dfIn_clean$Mass_AVONET)

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

dfIn_clean$FemaleSong_Agg01 <- as.factor(dfIn_clean$FemaleSong_Agg01)
dfIn_clean$HighConfidence_Coop <- as.factor(dfIn_clean$HighConfidence_Coop)
dfIn_clean$Territory1 <- as.factor(dfIn_clean$Territory1)
dfIn_clean$Territory2 <- as.factor(dfIn_clean$Territory2)
dfIn_clean$Territory3 <- as.factor(dfIn_clean$Territory3)


dfIn_clean$TerritorialityWeakVsStrong <- as.character(dfIn_clean$TerritorialityWeakVsStrong)
dfIn_clean$TerritorialityWeakVsStrongHighConf <- as.character(dfIn_clean$TerritorialityWeakVsStrongHighConf)
dfIn_clean$TerritorialityPermissiveExclusive <- as.character(dfIn_clean$TerritorialityPermissiveExclusive)
dfIn_clean$TerritorialityPermissiveExclusiveHighConf <- as.character(dfIn_clean$TerritorialityPermissiveExclusiveHighConf)
dfIn_clean$TerritorialityPermissiveColonialCoopVsExclusive <- as.character(dfIn_clean$TerritorialityPermissiveColonialCoopVsExclusive)

tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfIn_clean$species))
rownames(dfIn_clean) <- dfIn_clean$species

Territoriality_sym <- sym("TerritorialityWeakVsStrong")
col2_sym <- sym(columns[2])
col3_sym <- sym(columns[3])

# Get counts for this group
dfIn_clean %>% 
  group_by(!!Territoriality_sym) %>% 
  count() 

#### simple models ----
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
# 

models <- define_model_set(
  r1 = c(FemaleSong_Agg01 ~ !!Territoriality_sym),
  r2 = c(HighConfidence_Coop ~ !!Territoriality_sym)#,
  # r3 = c(FemaleSong_Agg01 ~ Territoriality_sym,
  #        HighConfidence_Coop ~ Territoriality_sym),
  # r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territoriality_sym),
  # r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + Territoriality_sym),
  # r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
  #        HighConfidence_Coop ~ Territoriality_sym),
  # r7 = c(FemaleSong_Agg01 ~ Territoriality_sym,
  #        HighConfidence_Coop ~ FemaleSong_Agg01),
  # r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  # r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01)#,
  # r10 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop,
  #         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  # r11 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
  #        HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong)
)
result <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'lambda',
  na.rm = TRUE
)

#### uses TerrWeakVsStrong and logMass ----
models <- define_model_set(
  rn1 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong),
  rn2 = c(HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  rn3 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r7 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r11 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong),
  rn1b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET),
  rn2b = c(HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  rn3b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong + logMass_AVONET),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + logMass_AVONET),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r7b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
         HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  rn7b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
          HighConfidence_Coop ~ logMass_AVONET),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  rn10b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r11b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + logMass_AVONET),
  rn11b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  rn3c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r7c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  r11c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + logMass_AVONET),
  rn3d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong,
          TerritorialityWeakVsStrong ~ logMass_AVONET),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong,
          TerritorialityWeakVsStrong ~ logMass_AVONET),
  r7d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          TerritorialityWeakVsStrong ~ logMass_AVONET),
  rn7d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
          TerritorialityWeakVsStrong ~ logMass_AVONET),
  r10d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong,
           TerritorialityWeakVsStrong ~ logMass_AVONET),
  rn10d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong,
           TerritorialityWeakVsStrong ~ logMass_AVONET),
  r11d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong,
           TerritorialityWeakVsStrong ~ logMass_AVONET),
  rn11d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong,
           TerritorialityWeakVsStrong ~ logMass_AVONET),
  n1 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong + logMass_AVONET),
  n2 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  n3 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ logMass_AVONET),
  n4 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  n5 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ logMass_AVONET)
)

#### uses TerrWeakVsStrongHighConf and logMass ----
models <- define_model_set(
  rn1 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf),
  rn2 = c(HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  rn3 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
          HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrongHighConf),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  r7 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  r11 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf),
  rn1b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn2b = c(HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn3b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r7b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  rn7b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
           HighConfidence_Coop ~ logMass_AVONET),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn10b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r11b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn11b = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn3c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r7c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  r11c = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  rn3d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf,
           TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf,
          TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  r7d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  rn7d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
           TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  r10d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf,
           TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  rn10d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf,
            TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  r11d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrongHighConf,
           TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  rn11d = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf,
            TerritorialityWeakVsStrongHighConf ~ logMass_AVONET),
  n1 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf + logMass_AVONET),
  n2 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  n3 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
         HighConfidence_Coop ~ logMass_AVONET),
  n4 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrongHighConf,
         HighConfidence_Coop ~ TerritorialityWeakVsStrongHighConf),
  n5 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ logMass_AVONET)
)


#### uses TerritorialityPermissiveExclusive and logMass ----
models <- define_model_set(
  rn1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive),
  rn2 = c(HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  rn3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveExclusive),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  r7 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  r11 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive),
  rn1b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  rn2b = c(HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  rn3b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveExclusive + logMass_AVONET),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive + logMass_AVONET),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  rn7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
           HighConfidence_Coop ~ logMass_AVONET),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  rn10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive + logMass_AVONET),
  rn11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  rn3c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r7c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  r11c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive + logMass_AVONET),
  rn3d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive,
           TerritorialityPermissiveExclusive ~ logMass_AVONET),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusive,
          TerritorialityPermissiveExclusive ~ logMass_AVONET),
  r7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          TerritorialityPermissiveExclusive ~ logMass_AVONET),
  rn7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
           TerritorialityPermissiveExclusive ~ logMass_AVONET),
  r10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusive,
           TerritorialityPermissiveExclusive ~ logMass_AVONET),
  rn10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusive,
            TerritorialityPermissiveExclusive ~ logMass_AVONET),
  r11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusive,
           TerritorialityPermissiveExclusive ~ logMass_AVONET),
  rn11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusive,
            TerritorialityPermissiveExclusive ~ logMass_AVONET),
  n1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusive + logMass_AVONET),
  n2 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  n3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
         HighConfidence_Coop ~ logMass_AVONET),
  n4 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusive,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusive),
  n5 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ logMass_AVONET)
)

#### uses TerritorialityPermissiveExclusiveHighConf and logMass ----
models <- define_model_set(
  rn1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf),
  rn2 = c(HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  rn3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveExclusiveHighConf),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  r7 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  r11 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf),
  rn1b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn2b = c(HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn3b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  rn7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
           HighConfidence_Coop ~ logMass_AVONET),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn3c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r7c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  r11c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  rn3d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf,
           TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf,
          TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  r7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  rn7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
           TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  r10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf,
           TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  rn10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf,
            TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  r11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveExclusiveHighConf,
           TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  rn11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf,
            TerritorialityPermissiveExclusiveHighConf ~ logMass_AVONET),
  n1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf + logMass_AVONET),
  n2 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  n3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
         HighConfidence_Coop ~ logMass_AVONET),
  n4 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveExclusiveHighConf,
         HighConfidence_Coop ~ TerritorialityPermissiveExclusiveHighConf),
  n5 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ logMass_AVONET)
)

#### uses TerritorialityPermissiveColonialCoopVsExclusive and logMass ----
models <- define_model_set(
  rn1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive),
  rn2 = c(HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  rn3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
          HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveColonialCoopVsExclusive),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  r7 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  r11 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive),
  rn1b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn2b = c(HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn3b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  rn7b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
           HighConfidence_Coop ~ logMass_AVONET),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn10b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn11b = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn3c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r7c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
          HighConfidence_Coop ~ FemaleSong_Agg01 + logMass_AVONET),
  r10c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + HighConfidence_Coop,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  r11c = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  rn3d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive,
           TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + logMass_AVONET,
          HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive,
          TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  r7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  rn7d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
           TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  r10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + HighConfidence_Coop + logMass_AVONET,
           HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive,
           TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  rn10d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive,
            TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  r11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
           HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityPermissiveColonialCoopVsExclusive,
           TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  rn11d = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
            HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive,
            TerritorialityPermissiveColonialCoopVsExclusive ~ logMass_AVONET),
  n1 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive + logMass_AVONET),
  n2 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  n3 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
         HighConfidence_Coop ~ logMass_AVONET),
  n4 = c(FemaleSong_Agg01 ~ TerritorialityPermissiveColonialCoopVsExclusive,
         HighConfidence_Coop ~ TerritorialityPermissiveColonialCoopVsExclusive),
  n5 = c(FemaleSong_Agg01 ~ logMass_AVONET,
         HighConfidence_Coop ~ logMass_AVONET)
)

#### uses longevity and territory 12vs3 ----
models <- define_model_set(
  r1 = c(FemaleSong_Agg01 ~ Territory_12vs3),
  r2 = c(HighConfidence_Coop ~ Territory_12vs3),
  r3 = c(FemaleSong_Agg01 ~ Territory_12vs3,
         HighConfidence_Coop ~ Territory_12vs3),
  r4 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3),
  r5 = c(HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3),
  r6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ Territory_12vs3),
  r7 = c(FemaleSong_Agg01 ~ Territory_12vs3,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  r9 = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  r10 = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop,
          HighConfidence_Coop ~ Territory_12vs3),
  r11 = c(FemaleSong_Agg01 ~ Territory_12vs3,
          HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3),
  r1b = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity),
  r2b = c(HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r3b = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity,
          HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r4b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory_12vs3 + log_Longevity),
  r5b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3 + log_Longevity),
  r6b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + log_Longevity,
          HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r7b = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity,
          HighConfidence_Coop ~ FemaleSong_Agg01 + log_Longevity),
  r8b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + log_Longevity),
  r9b = c(HighConfidence_Coop ~ FemaleSong_Agg01 + log_Longevity),
  r10b = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop + log_Longevity,
           HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r11b = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity,
           HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3 + log_Longevity),
  r3c = c(FemaleSong_Agg01 ~ Territory_12vs3,
          HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r6c = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
          HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r7c = c(FemaleSong_Agg01 ~ Territory_12vs3,
          HighConfidence_Coop ~ FemaleSong_Agg01 + log_Longevity),
  r10c = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop,
           HighConfidence_Coop ~ Territory_12vs3 + log_Longevity),
  r11c = c(FemaleSong_Agg01 ~ Territory_12vs3,
           HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3 + log_Longevity),
  r3d = c(FemaleSong_Agg01 ~ Territory_12vs3,
          HighConfidence_Coop ~ Territory_12vs3,
          Territory_12vs3 ~ log_Longevity),
  r6d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + log_Longevity,
          HighConfidence_Coop ~ Territory_12vs3,
          Territory_12vs3 ~ log_Longevity),
  r7d = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity,
          HighConfidence_Coop ~ FemaleSong_Agg01,
          Territory_12vs3 ~ log_Longevity),
  r10d = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop + log_Longevity,
           HighConfidence_Coop ~ Territory_12vs3,
           Territory_12vs3 ~ log_Longevity),
  r11d = c(FemaleSong_Agg01 ~ Territory_12vs3 + log_Longevity,
           HighConfidence_Coop ~ FemaleSong_Agg01 + Territory_12vs3,
           Territory_12vs3 ~ log_Longevity)
)
dfIn_clean$Territory_12vs3 = as.factor(dfIn_clean$Territory_12vs3)
result <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'lambda',
  na.rm = TRUE
)


#### uses Terr1, Terr3, and TerrWeakVsStrong ----
models <- define_model_set(
  m0 = c(FemaleSong_Agg01 ~ Territory1 + Territory3),
  m0b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory1),
  m0c = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3),
  m0d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory1),
  m1b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory1,
          HighConfidence_Coop ~ Territory1 + Territory3),
  m1c = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3,
          HighConfidence_Coop ~ Territory1 + Territory3),
  m1 = c(FemaleSong_Agg01 ~ Territory1 + Territory3 + HighConfidence_Coop, 
         HighConfidence_Coop ~ Territory1 + Territory3),
  m2 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ Territory1 + Territory3),
  m3 = c(FemaleSong_Agg01 ~ Territory1 + Territory3,
         HighConfidence_Coop ~ Territory1 + Territory3),
  m4 = c(FemaleSong_Agg01 ~ Territory1 + Territory3,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  m5 = c(FemaleSong_Agg01 ~ Territory1 + Territory3,
         HighConfidence_Coop ~ FemaleSong_Agg01 + Territory1 + Territory3),
  m6 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory1,
         HighConfidence_Coop ~ Territory1 + Territory3),
  m7 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3,
         HighConfidence_Coop ~ Territory1 + Territory3),
  m8 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory1,
         HighConfidence_Coop ~ Territory1),
  m9 = c(FemaleSong_Agg01 ~ HighConfidence_Coop + Territory3 + Territory1,
         HighConfidence_Coop ~ Territory3),
  r0 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong),
  r0d = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong),
  r1b = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong,
          HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r1 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r2 = c(FemaleSong_Agg01 ~ HighConfidence_Coop,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r3 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ TerritorialityWeakVsStrong),
  r4 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ FemaleSong_Agg01),
  r5 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong,
         HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong),
  mr5 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1 + Territory3,
         HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong),
  mr6 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1 + Territory3,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + Territory1 + Territory3),
  mr7 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory3,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + Territory1),
  mr8 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + Territory3),
  mr9 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + Territory1),
  mr10 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory3,
          HighConfidence_Coop ~ FemaleSong_Agg01 + TerritorialityWeakVsStrong + Territory3),
  mr11 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory3,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory3),
  mr12 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory3 + Territory1,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory3 + Territory1),
  mr13 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory1),
  mr14 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory3,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory1),
  mr15 = c(FemaleSong_Agg01 ~ TerritorialityWeakVsStrong + Territory1,
           HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory3),
  mr10r = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong + Territory3,
    HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory3),
  mr6r = c(FemaleSong_Agg01 ~ HighConfidence_Coop + TerritorialityWeakVsStrong + Territory1 + Territory3,
            HighConfidence_Coop ~ TerritorialityWeakVsStrong + Territory1 + Territory3)
)

#### plots, analyses ----
models$m0
plot(models$m0)
plot(models$m1)
plot(models$m2)
plot(models$m4)
plot(models$m5)
plot(models$mr6)
plot(models$m3)

pdf(paste0("plot phylopath models TerrWS_MassAVONET_",Sys.Date(),".pdf"), height = 55, width = 55)
plot_model_set(models, text_size = 4)
dev.off()


result <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'lambda',
  na.rm = TRUE
)
result

result_logisticMPLE <- phylo_path(
  models, 
  data = dfIn_clean, 
  tree = tree_clean, 
  model = 'logistic_MPLE',
  na.rm = TRUE
)
result_logisticMPLE

(s <- summary(result))
(sLog <- summary(result_logisticMPLE))

plot(s)
plot(sLog)


(best_model <- best(result))
plot(average(result), text_size = 3)

(best_model <- best(result_logisticMPLE))
plot(average(result_logisticMPLE), text_size = 3)

plot(choice(result, choice = "r11d"), text_size = 3)
plot(choice(result, choice = "r10d"), text_size = 3)
plot(choice(result, choice = "r10b"), text_size = 3)
plot(choice(result, choice = "r4b"), text_size = 3)
plot(choice(result, choice = "r11b"), text_size = 3)

plot(best_model, text_size = 3)


plot_model_set(models[c("r7b", "r7d", "r10d", "r10b", "r11b", "r11d")], text_size = 3)
plot_model_set(models[c("r10b", "r11b", "r10d", "r11d", "rn3b", "rn10b", "rn11b", "n1", "rn10d", "rn11d")], text_size = 3)

est_DAG(models["r11b"], data = dfIn_clean, tree = tree_clean, model = "lambda")
DAG(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop + log_Longevity,
    HighConfidence_Coop ~ Territory_12vs3,
    Territory_12vs3 ~ log_Longevity)

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
