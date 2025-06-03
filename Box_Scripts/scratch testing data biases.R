# scratch quantify biases

library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)


boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-04-29.csv')
newdata = "Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-07.csv"
newdata = "Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv"
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')
tree = read.nexus(treefile)

allcols = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/2024-8-4AllCol_Duet_Female_AllTraits_fromJiaying2025-03-04.csv")
AVONET = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/AVONET Supplementary dataset 1_Sheet_Raw_Data.csv")
unique(AVONET$Sex)
table(AVONET$Sex)
table(AVONET$Data.type)
AVONET %>% group_by()


df = read.csv(newdata)
df = merge(df, allcols[,c("BirdtreeSpecies", "Female_plumage_score_Dale2015", "Male_plumage_score_Dale2015", "Realm_Jetz2011", "Region_Cockburn2006", "FemalePC1sum_Dunn2015", "FemalePC2sum_Dunn2015", "MalePC1sum_Dunn2015", "MalePC2sum_Dunn2015", "sumDiffPC1_Dunn2015", "sumDiffPC2_Dunn2015" )], by.x = "species", by.y = "BirdtreeSpecies")

# Subset to just Oscine species
dfOs = df[which(df$species %in% tree$tip.label),]
dfOs$HaveFSData = !is.na(dfOs$FemaleSong_Agg01)
dfOs$HaveCBData = !is.na(dfOs$HighConfidence_Coop)
dfOs$MassDiff_AVONET = abs(dfOs$Male_AVONET - dfOs$Female_AVONET)+0.01
logMassDiff_AVONET = log(dfOs$MassDiff_AVONET)

# Realm (from Jetz and Rubenstein 2011) - some regions have been studied more intensively. If certain regions are understudied, species would be less likely to have a female song and/or cooperative breeding classification there
# Realm; AT - Afrotropics, PA - Palearctic, NA - Nearctic, NT - Neotropics, IM - Indomalay, AA - Australasia, OC - Oceania). Holarctic species are significantly more likely to have evidence for or against female song than Tropical species, but they are not significantly more likely to have evidence for or against cooperative breeding. 
dfOs$GeographicRegion_Jetz = NA
dfOs$GeographicRegion_Jetz[which(dfOs$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
dfOs$GeographicRegion_Jetz[which(dfOs$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
dfOs %>% group_by(GeographicRegion_Jetz) %>% count

# Region (from Cockburn 2006) - "Nearctic", "Holarctic", "Palearctic" = "Holarctic"; "Indomalayan" "Australia"   "Neotropical" "Africa"  "Widespread"  = Tropical
dfOs$GeographicRegion_Cockburn = NA
dfOs$GeographicRegion_Cockburn[which(dfOs$Region_Cockburn2006 %in% c("Nearctic", "Holarctic", "Palearctic"))] <- "Holarctic"
dfOs$GeographicRegion_Cockburn[which(dfOs$Region_Cockburn2006 %in% c("Indomalayan", "Australia",   "Neotropical", "Africa"))] <- "Tropical"

#### Biases we checked  ----
# Cooperative breeding - if birds that are cooperative are more likely to be studied, we should see more FS classifications in cooperatively breeding birds. We do not. 
table_data <- dfOs %>%
  filter(!is.na(HighConfidence_Coop)) %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$HighConfidence_Coop
mat
chisq.test(mat)

# Colonial - if birds that group together are generally easier to study, we should see more FS and CB classifications in Colonial birds. We do not, for either.
# Female song data presence
dfOs %>% filter(!is.na(colonial_Griesser2023)) %>% group_by(colonial_Griesser2023, HaveFSData) %>% count
table_data <- dfOs %>%
  filter(!is.na(colonial_Griesser2023)) %>%
  group_by(colonial_Griesser2023, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$colonial_Griesser2023
mat
chisq.test(mat)

# Cooperative breeding data presence
dfOs %>% filter(!is.na(colonial_Griesser2023)) %>% group_by(colonial_Griesser2023, HaveCBData) %>% count
table_data <- dfOs %>%
  filter(!is.na(colonial_Griesser2023)) %>%
  group_by(colonial_Griesser2023, HaveCBData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveCBData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$colonial_Griesser2023
mat
chisq.test(mat)


# Familial living - if birds that group together are generally easier to study, we should see more FS and CB classifications in familial-living birds. In fact, familial-living birds are significantly less likely to have Female Song classifications, and also less likely to have Cooperative breeding classifications
# Female song data presence
dfOs %>% filter(!is.na(Griesser2017FamilialLiving)) %>% group_by(Griesser2017FamilialLiving, HaveFSData) %>% count
table_data <- dfOs %>%
  filter(!is.na(Griesser2017FamilialLiving)) %>%
  group_by(Griesser2017FamilialLiving, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Griesser2017FamilialLiving
mat
chisq.test(mat)

# Cooperative breeding data presence
dfOs %>% filter(!is.na(Griesser2017FamilialLiving)) %>% group_by(Griesser2017FamilialLiving, HaveCBData) %>% count
table_data <- dfOs %>%
  filter(!is.na(Griesser2017FamilialLiving)) %>%
  group_by(Griesser2017FamilialLiving, HaveCBData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveCBData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Griesser2017FamilialLiving
mat
chisq.test(mat)


# Breeding system - if birds that have interesting breeding systems are more likely to be studied intensely enough to have FS or CB classifications, we should see more classifications in polygynous species. However, there is no significant difference in FS and CB data availability between socially monogamous species and polygynous species
# Female song data presence
dfOs %>% filter(!is.na(Final.polygyny)) %>% group_by(Final.polygyny, HaveFSData) %>% count
table_data <- dfOs %>%
  filter(!is.na(Final.polygyny)) %>%
  group_by(Final.polygyny, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Final.polygyny
mat
chisq.test(mat)

# Female song data presence
dfOs %>% filter(!is.na(Final.polygyny)) %>% group_by(Final.polygyny, HaveCBData) %>% count
table_data <- dfOs %>%
  filter(!is.na(Final.polygyny)) %>%
  group_by(Final.polygyny, HaveCBData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveCBData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Final.polygyny
mat
chisq.test(mat)

#### Biases we checked that appear to be present in our study ----
# Year-round territoriality from Tobias et al 2016 - if birds that are year-round territorial are more likely to be studied, we should see more FS and CB classifications in year-round territorial birds. Indeed we do. 
table_data <- dfOs %>%
  filter(!is.na(Territory_12vs3)) %>%
  group_by(Territory_12vs3, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Territory_12vs3
mat
chisq.test(mat)

table_data <- dfOs %>%
  filter(!is.na(Territory_12vs3)) %>%
  group_by(Territory_12vs3, HaveCBData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveCBData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$Territory_12vs3
mat
chisq.test(mat)


# Geographic region
table_data <- dfOs %>%
  filter(!is.na(GeographicRegion)) %>%
  group_by(GeographicRegion, HaveFSData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$GeographicRegion
mat
chisq.test(mat)

table_data <- dfOs %>%
  filter(!is.na(GeographicRegion)) %>%
  group_by(GeographicRegion, HaveCBData) %>%
  count() %>%
  pivot_wider(
    names_from = HaveCBData,
    values_from = n,
    values_fill = 0
  )
# Convert to matrix for chisq.test
mat <- as.matrix(table_data[, -1])
rownames(mat) <- table_data$GeographicRegion
mat
chisq.test(mat)


# Plumage dimorphism (Dale 2015) - test whether sexually monochromatic species are less likely to have female song data or more likely to have a classification of female song absent
# tests
dfOs$MaleFemalePlumageDiff = dfOs$Male_plumage_score_Dale2015-dfOs$Female_plumage_score_Dale2015
dfOs$MaleFemalePlumageDiffAbs = abs(dfOs$MaleFemalePlumageDiff) # absolute value of difference Male-Female

# Species with female song classifications do have significantly higher differences between male and female plumage coloration (both not-absolute-value and absolute value)
boxplot(formula = MaleFemalePlumageDiff ~ HaveFSData, data = dfOs)
t.test(formula = MaleFemalePlumageDiff ~ HaveFSData, data = dfOs)

boxplot(formula = MaleFemalePlumageDiffAbs ~ HaveFSData, data = dfOs)
t.test(formula = MaleFemalePlumageDiffAbs ~ HaveFSData, data = dfOs)

# Species with a classification of Female Song Absent do have significantly higher differences between male and female plumage coloration (both not-absolute-value and absolute value)
boxplot(formula = MaleFemalePlumageDiff ~ FemaleSong_Agg01, data = dfOs)
t.test(formula = MaleFemalePlumageDiff ~ FemaleSong_Agg01, data = dfOs)

boxplot(formula = MaleFemalePlumageDiffAbs ~ FemaleSong_Agg01, data = dfOs)
t.test(formula = MaleFemalePlumageDiffAbs ~ FemaleSong_Agg01, data = dfOs)

# Repeat for cooperative breeding:
boxplot(formula = MaleFemalePlumageDiffAbs ~ HaveCBData, data = dfOs)
t.test(formula = MaleFemalePlumageDiffAbs ~ HaveCBData, data = dfOs)


# plots - might not use these
dfOs$FemaleSong_Agg01 <- as.character(dfOs$FemaleSong_Agg01)
ggplot(dfOs, mapping = aes(x = Female_plumage_score_Dale2015, y = Male_plumage_score_Dale2015)) +
  geom_point(aes(col = HaveFSData), size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

ggplot(dfOs %>% filter(HaveFSData == T), mapping = aes(x = Female_plumage_score_Dale2015, y = Male_plumage_score_Dale2015)) +
  geom_point(aes(col = HaveFSData), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

ggplot(dfOs %>% filter(HaveFSData == F), mapping = aes(x = Female_plumage_score_Dale2015, y = Male_plumage_score_Dale2015)) +
  geom_point(aes(col = HaveFSData), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

ggplot(dfOs %>% filter(HaveFSData == T), mapping = aes(x = Female_plumage_score_Dale2015, y = Male_plumage_score_Dale2015)) +
  geom_point(aes(col = FemaleSong_Agg01), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

ggplot(dfOs, mapping = aes(x = Female_plumage_score_Dale2015, y = Male_plumage_score_Dale2015)) +
  geom_point(aes(col = FemaleSong_Agg01), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line



## Repeat tests for Dunn et al (2015) dichromatism measures. PC1 corresponds with sex difference in brightness; PC2 corresponds with sex difference in hue
boxplot(formula = sumDiffPC1_Dunn2015 ~ HaveFSData, data = dfOs)
t.test(formula = sumDiffPC1_Dunn2015 ~ HaveFSData, data = dfOs)

boxplot(formula = sumDiffPC2_Dunn2015 ~ HaveFSData, data = dfOs)
t.test(formula = sumDiffPC2_Dunn2015 ~ HaveFSData, data = dfOs)


dfOs$FemaleSong_Agg01 <- as.character(dfOs$FemaleSong_Agg01)
ggplot(dfOs, mapping = aes(x = MalePC1sum_Dunn2015, y = FemalePC1sum_Dunn2015)) +
  geom_point(aes(col = HaveFSData), size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

dfOs$FemaleSong_Agg01 <- as.character(dfOs$FemaleSong_Agg01)
ggplot(dfOs, mapping = aes(x = MalePC2sum_Dunn2015, y = FemalePC2sum_Dunn2015)) +
  geom_point(aes(col = HaveFSData), size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

dfOs$FemaleSong_Agg01 <- as.character(dfOs$FemaleSong_Agg01)
ggplot(dfOs %>% filter(HaveFSData == T), mapping = aes(x = MalePC1sum_Dunn2015, y = FemalePC1sum_Dunn2015)) +
  geom_point(aes(col = FemaleSong_Agg01), size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line

dfOs$FemaleSong_Agg01 <- as.character(dfOs$FemaleSong_Agg01)
ggplot(dfOs %>% filter(HaveFSData == T), mapping = aes(x = MalePC2sum_Dunn2015, y = FemalePC2sum_Dunn2015)) +
  geom_point(aes(col = FemaleSong_Agg01), size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # y = x line



## Size dimorphism
boxplot(formula = logMassDiff_AVONET ~ HaveFSData, data = dfOs)
t.test(formula = logMassDiff_AVONET ~ HaveFSData, data = dfOs)

boxplot(formula = logMassDiff_AVONET ~ HaveCBData, data = dfOs)
t.test(formula = logMassDiff_AVONET ~ HaveCBData, data = dfOs)



# Create output folder
dir.create("BiasTestMatrices", showWarnings = FALSE)

collect_chisq <- function(data, trait_col, data_col) {
  # Build contingency table
  table_data <- data %>%
    filter(!is.na(.data[[trait_col]])) %>%
    group_by(across(all_of(c(trait_col, data_col)))) %>%
    count() %>%
    pivot_wider(names_from = all_of(data_col), values_from = n, values_fill = 0)
  
  mat <- as.matrix(table_data[,-1])
  rownames(mat) <- table_data[[trait_col]]
  
  # Choose appropriate test
  if (any(mat < 5)) {
    test <- fisher.test(mat)
    test_type <- "Fisher"
  } else {
    test <- chisq.test(mat)
    test_type <- "Chi-squared"
  }
  test <- chisq.test(mat)
  test_type <- "Chi-squared"
  
  print(test)
  
  # Prepare labeled matrix for CSV output
  df_mat <- as.data.frame(mat)
  df_mat <- tibble::rownames_to_column(df_mat, var = trait_col)
  
  # Create header row with column label above column names
  col_header <- c(trait_col, rep(data_col, ncol(df_mat) - 1))
  names(df_mat) <- colnames(df_mat)
  final_mat <- rbind(col_header, colnames(df_mat), df_mat)
  
  # Save to CSV
  matrix_filename <- paste0("BiasTestMatrices/", trait_col, "_", data_col, "_matrix.csv")
  write.table(final_mat, file = matrix_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Return test summary
  tibble(
    Trait = trait_col,
    DataType = data_col,
    Test = test_type,
    p_value = round(test$p.value, 4),
    X2 = ifelse(test_type == "Chi-squared", round(test$statistic, 2), NA),
    df = ifelse(test_type == "Chi-squared", test$parameter, NA)
  )
}



# Geographic region transformation
dfOs$GeographicRegion <- NA
dfOs$GeographicRegion[dfOs$Realm_Jetz2011 %in% c("PA", "NeA")] <- "Holarctic"
dfOs$GeographicRegion[dfOs$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC")] <- "Tropical"


## Interaction tests
summary(glm(HaveFSData ~ HighConfidence_Coop * GeographicRegion, data = dfOs, family = binomial))
summary(glm(HaveFSData ~ HighConfidence_Coop * Griesser2017FamilialLiving, data = dfOs, family = binomial))
summary(glm(HaveFSData ~ HighConfidence_Coop * Territory_12vs3, data = dfOs, family = binomial))
summary(glm(HaveCBData ~ FemaleSong_Agg01 * GeographicRegion, data = dfOs, family = binomial))
summary(glm(HaveCBData ~ FemaleSong_Agg01 * Griesser2017FamilialLiving, data = dfOs, family = binomial))
summary(glm(HaveCBData ~ FemaleSong_Agg01 * Territory_12vs3, data = dfOs, family = binomial))



# summary(loglm( ~ HaveFSData * HighConfidence_Coop * Territory_12vs3, data = dfOs))

dfOs$FemaleSong_Agg01 <- as.numeric(dfOs$FemaleSong_Agg01)
summary(glm(FemaleSong_Agg01 ~ HighConfidence_Coop * GeographicRegion, data = dfOs %>% filter(HaveFSData == T), family = binomial))
summary(glm(FemaleSong_Agg01 ~ HighConfidence_Coop * Griesser2017FamilialLiving, data = dfOs %>% filter(HaveFSData == T), family = binomial))
summary(glm(FemaleSong_Agg01 ~ HighConfidence_Coop * Territory_12vs3, data = dfOs %>% filter(HaveFSData == T), family = binomial))


results <- bind_rows(
  collect_chisq(dfOs, "HighConfidence_Coop", "HaveFSData"),
  collect_chisq(dfOs, "colonial_Griesser2023", "HaveFSData"),
  collect_chisq(dfOs, "colonial_Griesser2023", "HaveCBData"),
  collect_chisq(dfOs, "Griesser2017FamilialLiving", "HaveFSData"),
  collect_chisq(dfOs, "Griesser2017FamilialLiving", "HaveCBData"),
  collect_chisq(dfOs, "Final.polygyny", "HaveFSData"),
  collect_chisq(dfOs, "Final.polygyny", "HaveCBData"),
  collect_chisq(dfOs, "Territory_12vs3", "HaveFSData"),
  collect_chisq(dfOs, "Territory_12vs3", "HaveCBData"),
  collect_chisq(dfOs, "TerritorialityWeakVsStrong", "HaveFSData"),
  collect_chisq(dfOs, "TerritorialityWeakVsStrong", "HaveCBData"),
  collect_chisq(dfOs, "TerritorialityWeakVsStrongHighConf", "HaveFSData"),
  collect_chisq(dfOs, "TerritorialityWeakVsStrongHighConf", "HaveCBData"),
  collect_chisq(dfOs, "TerritorialityPermissiveColonialCoopVsExclusive", "HaveFSData"),
  collect_chisq(dfOs, "TerritorialityPermissiveColonialCoopVsExclusive", "HaveCBData"),
  collect_chisq(dfOs, "GeographicRegion_Jetz", "HaveFSData"),
  collect_chisq(dfOs, "GeographicRegion_Cockburn", "HaveCBData")
)

write.csv(results, paste(Sys.Date(), "Bias_ChiSquare_Results.csv"), row.names = FALSE)
