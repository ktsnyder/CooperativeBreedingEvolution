# Merge columns from R data into supplementary information data

library(readxl)

submission1_table <- read_excel('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - RESUBMISSION1/Supplementary Information - Dataset 2 - Cooperative Breeding and Sociality Data.xlsx', sheet = 1)
r_data = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Data_R_2025-07-23_wNormCaretakers_SocBonds_GroupNorm.csv")
colnames(r_data)[which(colnames(r_data) == "Territory")] <- "Territory_Tobias2016"

merge(submission1_table, r_data[,c("species", "Territory_Tobias2016", "Mass_AVONET", "Centroid.Latitude_AVONET", "Migration_AVONET", "log_Wing.Length_mean_F_AVONET", "log_Wing.Length_mean_M_AVONET", "Realm_Jetz2011", "Female_plumage_score_Dale2015", "Male_plumage_score_Dale2015")])