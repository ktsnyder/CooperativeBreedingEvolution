# testing data biases - organized outputs
# 

library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(broom)
library(gt)
library(knitr)
library(kableExtra)



# Update these paths to your local directory structure
basepath = getwd()  # Assumes running from project root

newdata = file.path(basepath, 'Box_Scripts/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')
treefile = file.path(basepath, 'ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')
tree = read.nexus(treefile)

# NOTE: The allcols file needs to be obtained from Box separately if needed
# allcols = read.csv("path/to/2024-8-4AllCol_Duet_Female_AllTraits_fromJiaying2025-03-04.csv")

df = read.csv(newdata)
# NOTE: Merge with allcols commented out - uncomment when allcols file is available
# df = merge(df, allcols[,c("BirdtreeSpecies", "Female_plumage_score_Dale2015", "Male_plumage_score_Dale2015", "Realm_Jetz2011", "FemalePC1sum_Dunn2015", "FemalePC2sum_Dunn2015", "MalePC1sum_Dunn2015", "MalePC2sum_Dunn2015", "sumDiffPC1_Dunn2015", "sumDiffPC2_Dunn2015", "Tropical_life_history_ppca_Dale2015", "Region_Cockburn2006", "location_Downing2015")], by.x = "species", by.y = "BirdtreeSpecies")

# Subset to just Oscine species
dfOs = df[which(df$species %in% tree$tip.label),]
dfOs$HaveFSData = !is.na(dfOs$FemaleSong_Agg01)
dfOs$HaveCBData = !is.na(dfOs$HighConfidence_Coop)

dfOs %>% filter(!is.na(HighConfidence_Coop) & !is.na(FemaleSong_Agg01)) %>% group_by(Realm_Jetz2011, Region_Cockburn2006) %>% count %>% print(n=40)

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
#dfOs$GeographicRegion_Cockburn[which(dfOs$Region_Cockburn2006 == "Widespread")] <- "Widespread"
dfOs %>% group_by(GeographicRegion_Cockburn) %>% count
dfOs %>% filter(!is.na(HighConfidence_Coop)) %>% group_by(HaveFSData, GeographicRegion_Cockburn, GeographicRegion_Jetz) %>% count

dfOs %>% filter(GeographicRegion_Cockburn != GeographicRegion_Jetz) %>% group_by(GeographicRegion_Jetz,GeographicRegion_Cockburn, Realm_Jetz2011, Region_Cockburn2006) %>% count %>% print(n=40)

#dfOs$GeographicRegion = "Earth"


df_FS <- dfOs %>% filter(HaveFSData == TRUE)

dfOs %>% filter(HaveCBData) %>% group_by(HighConfidence_Coop, GeographicRegion_Cockburn, HaveFSData) %>% count

dfOs_char = dfOs
dfOs_char$HighConfidence_Coop = as.character(dfOs_char$HighConfidence_Coop)

#### calculate downsampling - geographic regions ----
# Create Holarctic matrix
holarctic_matrix <- dfOs_char %>%
  #filter(HaveCBData, GeographicRegion_Cockburn == "Holarctic") %>%
  filter(HaveCBData, GeographicRegion_Jetz == "Holarctic") %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    names_prefix = "HaveFSData_"
  ) %>%
  mutate(Total = HaveFSData_FALSE + HaveFSData_TRUE) %>%
  bind_rows(
    summarize(., 
              HighConfidence_Coop = "Total", 
              HaveFSData_FALSE = sum(HaveFSData_FALSE), 
              HaveFSData_TRUE = sum(HaveFSData_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSData_FALSE = HaveFSData_FALSE/Total, pct_HaveFSData_TRUE = HaveFSData_TRUE/Total) 

# Get totals for percentage calculations
total_row <- holarctic_matrix %>% filter(HighConfidence_Coop == "Total")

# Add percentage rows in one step
pct_rows <- holarctic_matrix %>%
  filter(HighConfidence_Coop %in% c("0", "1")) %>%
  mutate(
    HaveFSData_FALSE = HaveFSData_FALSE / total_row$HaveFSData_FALSE,
    HaveFSData_TRUE = HaveFSData_TRUE / total_row$HaveFSData_TRUE,
    Total = Total / total_row$Total,
    HighConfidence_Coop = paste0("pct_HighConfCoop_", HighConfidence_Coop),
    pct_HaveFSData_FALSE = NA,
    pct_HaveFSData_TRUE = NA
  )

# Combine everything
final_holarctic <- bind_rows(holarctic_matrix, pct_rows)

# Create Tropical matrix
tropical_matrix <- dfOs_char %>%
  #filter(HaveCBData, GeographicRegion_Cockburn == "Tropical") %>%
  filter(HaveCBData, GeographicRegion_Jetz == "Tropical") %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    names_prefix = "HaveFSData_"
  ) %>%
  mutate(Total = HaveFSData_FALSE + HaveFSData_TRUE) %>%
  bind_rows(
    summarize(., 
              HighConfidence_Coop = "Total", 
              HaveFSData_FALSE = sum(HaveFSData_FALSE), 
              HaveFSData_TRUE = sum(HaveFSData_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSData_FALSE = HaveFSData_FALSE/Total, pct_HaveFSData_TRUE = HaveFSData_TRUE/Total) 

# Get totals for percentage calculations
total_row <- tropical_matrix %>% filter(HighConfidence_Coop == "Total")

# Add percentage rows in one step
pct_rows <- tropical_matrix %>%
  filter(HighConfidence_Coop %in% c("0", "1")) %>%
  mutate(
    HaveFSData_FALSE = HaveFSData_FALSE / total_row$HaveFSData_FALSE,
    HaveFSData_TRUE = HaveFSData_TRUE / total_row$HaveFSData_TRUE,
    Total = Total / total_row$Total,
    HighConfidence_Coop = paste0("pct_HighConfCoop_", HighConfidence_Coop),
    pct_HaveFSData_FALSE = NA,
    pct_HaveFSData_TRUE = NA
  )

# Combine everything
final_tropical <- bind_rows(tropical_matrix, pct_rows)

## Matrix without splitting by Holarctic and Tropical
# Create Holarctic matrix
global_matrix <- dfOs_char %>%
  filter(HaveCBData) %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    names_prefix = "HaveFSData_"
  ) %>%
  mutate(Total = HaveFSData_FALSE + HaveFSData_TRUE) %>%
  bind_rows(
    summarize(., 
              HighConfidence_Coop = "Total", 
              HaveFSData_FALSE = sum(HaveFSData_FALSE), 
              HaveFSData_TRUE = sum(HaveFSData_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSData_FALSE = HaveFSData_FALSE/Total, pct_HaveFSData_TRUE = HaveFSData_TRUE/Total) 

# Get totals for percentage calculations
total_row <- global_matrix %>% filter(HighConfidence_Coop == "Total")

# Add percentage rows in one step
pct_rows <- global_matrix %>%
  filter(HighConfidence_Coop %in% c("0", "1")) %>%
  mutate(
    HaveFSData_FALSE = HaveFSData_FALSE / total_row$HaveFSData_FALSE,
    HaveFSData_TRUE = HaveFSData_TRUE / total_row$HaveFSData_TRUE,
    Total = Total / total_row$Total,
    HighConfidence_Coop = paste0("pct_HighConfCoop_", HighConfidence_Coop),
    pct_HaveFSData_FALSE = NA,
    pct_HaveFSData_TRUE = NA
  )

# Combine everything
final_global <- bind_rows(global_matrix, pct_rows)


# Print the matrices
print("Holarctic Matrix:")
print(final_holarctic)

print("Tropical Matrix:")
print(final_tropical)

print("Global Matrix:")
print(final_global)

## Downsample Holarctic (Jetz) Non-cooperative breeders
# Holarctic:
#   Numerator: [HighConfidence_Coop0, HaveFSData_TRUE] = 209
# Denominator: [HighConfidence_Coop0, Total] = 592
# 209/592 = 0.353
# Tropical:
#   Numerator: [HighConfidence_Coop0, HaveFSData_TRUE] = 692
# Denominator: [HighConfidence_Coop0, Total] = 3251
# 692/3251 = 0.213
# so more holarctic non-cooperative birds have FS data than tropical non cooperative birds. to downsample we want to include in our analysis x/592=0.213 so 126 birds, thus remove 209-126 = 83 birds

## Downsample Holarctic (Cockburn) Non-cooperative breeders
# Holarctic: 
#   Numerator: [HighConfidence_Coop0, HaveFSData_TRUE] = 234
#   Denominator: [HighConfidence_Coop0, Total] = 693
#   234/693 = 0.3377
# Tropical:
#   Numerator: [HighConfidence_Coop0, HaveFSData_TRUE] = 647
#   Denominator: [HighConfidence_Coop0, Total] = 3061
#   647/3061 = 0.2114
# x/693=0.2114 so x = 146.49 birds, thus remove 234-146 = 88 species


## Downsample Tropical Cooperative breeders - Jetz geographic region
# In this next case, we want to mitigate the overrepresentation of HaveFSData_TRUE in Tropical+Coop1 species.
# Because the denominator is the total Tropical species with HaveFSData_TRUE, the denominator is going to change if we’re removing the excess Coop=1 species with HaveFSData_TRUE.
# Tropical+Coop1+HaveFSData_TRUE = 112 species
# Tropical+HaveFSData_TRUE = 804 species
# Proportion of tropical species that are Cooperative = 0.113
#   (112-x)/(804-x) = 0.113
#   which works out to x = 23.8 (=24) species to remove (i.e. move from HaveFSData_TRUE to HaveFSData_FALSE). 

## Downsample Tropical Cooperative breeders - Cockburn geographic region
# In this next case, we want to mitigate the overrepresentation of HaveFSData_TRUE in Tropical+Coop1 species.
# Tropical+Coop1+HaveFSData_TRUE = 113 species
# Tropical+HaveFSData_TRUE = 760 species
# Proportion of tropical species that are Cooperative = 0.119
#   (113-x)/(760-x) = 0.119
#   which works out to x = 25.61 (=26) species to remove (i.e. move from HaveFSData_TRUE to HaveFSData_FALSE). 
#

## Downsample Global Cooperative breeders
# Want to mitigate the overrepresentation of HaveFSData_TRUE in Coop1 species.
# Global rate of Coop1 = 0.107
# Rate of Coop1 where HaveFSData_TRUE = 0.120
# Coop=1 + HaveFSData_TRUE = 125 species
# HaveFSData = TRUE (total) = 1041 species
# (125-x)/(1041-x) = 0.107
# which works out to x = 15.244 (=15) species to remove.


#### Calculate downsampling - Territory12vs3 ----
## Assumes that highly or year-round territorial birds are easier to study, and thus more likely to have data for FS and Coop
terr12_matrix <- dfOs_char %>%
  filter(HaveCBData, Territory_12vs3==0) %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    names_prefix = "HaveFSData_"
  ) %>%
  mutate(Total = HaveFSData_FALSE + HaveFSData_TRUE) %>%
  bind_rows(
    summarize(., 
              HighConfidence_Coop = "Total", 
              HaveFSData_FALSE = sum(HaveFSData_FALSE), 
              HaveFSData_TRUE = sum(HaveFSData_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSData_FALSE = HaveFSData_FALSE/Total, pct_HaveFSData_TRUE = HaveFSData_TRUE/Total) 

# Get totals for percentage calculations
total_row <- terr12_matrix %>% filter(HighConfidence_Coop == "Total")

# Add percentage rows in one step
pct_rows <- terr12_matrix %>%
  filter(HighConfidence_Coop %in% c("0", "1")) %>%
  mutate(
    HaveFSData_FALSE = HaveFSData_FALSE / total_row$HaveFSData_FALSE,
    HaveFSData_TRUE = HaveFSData_TRUE / total_row$HaveFSData_TRUE,
    Total = Total / total_row$Total,
    HighConfidence_Coop = paste0("pct_HighConfCoop_", HighConfidence_Coop),
    pct_HaveFSData_FALSE = NA,
    pct_HaveFSData_TRUE = NA
  )

# Combine everything
final_terr12 <- bind_rows(terr12_matrix, pct_rows)

## Territory=3 subset matrix
terr3_matrix <- dfOs_char %>%
  filter(HaveCBData, Territory_12vs3==1) %>%
  group_by(HighConfidence_Coop, HaveFSData) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSData,
    values_from = n,
    names_prefix = "HaveFSData_"
  ) %>%
  mutate(Total = HaveFSData_FALSE + HaveFSData_TRUE) %>%
  bind_rows(
    summarize(., 
              HighConfidence_Coop = "Total", 
              HaveFSData_FALSE = sum(HaveFSData_FALSE), 
              HaveFSData_TRUE = sum(HaveFSData_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSData_FALSE = HaveFSData_FALSE/Total, pct_HaveFSData_TRUE = HaveFSData_TRUE/Total) 

# Get totals for percentage calculations
total_row <- terr3_matrix %>% filter(HighConfidence_Coop == "Total")

# Add percentage rows in one step
pct_rows <- terr3_matrix %>%
  filter(HighConfidence_Coop %in% c("0", "1")) %>%
  mutate(
    HaveFSData_FALSE = HaveFSData_FALSE / total_row$HaveFSData_FALSE,
    HaveFSData_TRUE = HaveFSData_TRUE / total_row$HaveFSData_TRUE,
    Total = Total / total_row$Total,
    HighConfidence_Coop = paste0("pct_HighConfCoop_", HighConfidence_Coop),
    pct_HaveFSData_FALSE = NA,
    pct_HaveFSData_TRUE = NA
  )

# Combine everything
final_terr3 <- bind_rows(terr3_matrix, pct_rows)

## Alt approach - downsample based on data presence for both FS and CB
dfOs_char$HaveFSCBdata = !is.na(dfOs_char$HighConfidence_Coop) & !is.na(dfOs_char$FemaleSong_Agg01)
dfOs_char$Territory_12vs3 = as.character(dfOs_char$Territory_12vs3)

terr12v3_matrix <- dfOs_char %>%
  filter(!is.na(Territory_12vs3)) %>%
  group_by(Territory_12vs3, HaveFSCBdata) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = HaveFSCBdata,
    values_from = n,
    names_prefix = "HaveFSCBdata_"
  ) %>%
  mutate(Total = HaveFSCBdata_FALSE + HaveFSCBdata_TRUE) %>%
  bind_rows(
    summarize(., 
              Territory_12vs3 = "Total", 
              HaveFSCBdata_FALSE = sum(HaveFSCBdata_FALSE), 
              HaveFSCBdata_TRUE = sum(HaveFSCBdata_TRUE),
              Total = sum(Total))
  ) %>%
  mutate(pct_HaveFSCBdata_FALSE = HaveFSCBdata_FALSE/Total, pct_HaveFSCBdata_TRUE = HaveFSCBdata_TRUE/Total) 

# Get totals for percentage calculations
total_row <- terr12v3_matrix %>% filter(Territory_12vs3 == "Total")

# Add percentage rows in one step
pct_rows <- terr12v3_matrix %>%
  filter(Territory_12vs3 %in% c("0", "1")) %>%
  mutate(
    HaveFSCBdata_FALSE = HaveFSCBdata_FALSE / total_row$HaveFSCBdata_FALSE,
    HaveFSCBdata_TRUE = HaveFSCBdata_TRUE / total_row$HaveFSCBdata_TRUE,
    Total = Total / total_row$Total,
    Territory_12vs3 = paste0("pct_Terr12v3_", Territory_12vs3),
    pct_HaveFSCBdata_FALSE = NA,
    pct_HaveFSCBdata_TRUE = NA
  )

# Combine everything
final_terr12v3 <- bind_rows(terr12v3_matrix, pct_rows)


# Print the matrices
print("Territory12 Matrix:")
print(final_terr12)

print("Territory3 Matrix:")
print(final_terr3)

print("Territory12v3 Matrix:")
print(final_terr12v3)

## Downsample based on final_terr12v3
# it does seem that year-round territorial species are more likely to have both FS and CB data 
# (32.4% of Territory3 species have FS and CB classifications; 20.0% of Territory1or2 species have FS and CB classifications).
# species Territory_12vs3=1 (aka Territory=3) with FS/CB data: 404
# Total species with Territory_12vs3=1 (aka Territory=3): 1246
# 404/1246 = 0.324 (want 0.200)
# (404-x)/1246 = 0.200
# x = 154.8 (--> remove 155 year-round territorial species with FS/CB data from analyses)
# 


#### Calculate and run downsampling on multiple Territoriality measures ----
TerritorialityCols = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "TerritorialityPermissiveColonialCoopVsExclusive") 
for (i in 1:length(TerritorialityCols) ) {
  tempterr = TerritorialityCols[i]
  print(tempterr)
  dfOs_char[,tempterr] = as.character(dfOs_char[,tempterr])
  dfOs_char$HaveFSCBdata = !is.na(dfOs_char$FemaleSong_Agg01) & !is.na(dfOs_char$HighConfidence_Coop)
  downsampleout <- calculate_downsampling(dfOs_char, territoriality_col = tempterr, territory_value_high = "1", territory_value_low = "0", data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop"))
  print(downsampleout$report)
  print(downsampleout$filter_statement)
  numSpeciesToRemove = downsampleout$n_to_remove
  valueToDownsample = downsampleout$report$downsample_info$territory_value
  
  
  results <- run_multiple_phylopath(
    dfIn = dfIn,
    tree = tree,
    downsample_columns = c(tempterr, "HaveFSCBdata"),
    downsample_values = c(valueToDownsample, TRUE),
    numToRemove = numSpeciesToRemove,
    n_iterations = 500,
    female_song_var = "FemaleSong_Agg01",
    coop_breeding_var = "HighConfidence_Coop",
    territoriality_var = tempterr,
    mass_var = "logMass_AVONET",
    plotlabel = paste0(tempterr, " Remove", numSpeciesToRemove,"Territory",  valueToDownsample)
  )
  
}


# Function to calculate number of species to downsample
# Based on a specified territoriality column and data presence columns

calculate_downsampling <- function(df, 
                                   territoriality_col,
                                   territory_value_high = 1,  # Value for high territoriality (e.g., year-round)
                                   territory_value_low = 0,   # Value for low territoriality (e.g., seasonal) 
                                   data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop"),
                                   target_proportion = NULL) {
  
  # Validate inputs
  if (!territoriality_col %in% colnames(df)) {
    stop(paste("Territoriality column", territoriality_col, "not found in dataframe"))
  }
  
  for (col in data_cols) {
    if (!col %in% colnames(df)) {
      stop(paste("Data column", col, "not found in dataframe"))
    }
  }
  
  print(territoriality_col)
  
  df[,territoriality_col] <- as.character(df[,territoriality_col])
  
  # Create a variable for having complete data across specified columns
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  # Filter to rows with non-NA territory values
  df_filtered <- df[!is.na(df[[territoriality_col]]), ]
  
  # Calculate data availability for high and low territoriality
  terr_high <- df_filtered[df_filtered[[territoriality_col]] == territory_value_high, ]
  terr_low <- df_filtered[df_filtered[[territoriality_col]] == territory_value_low, ]
  
  # Count species with complete data for each territoriality type
  n_complete_high <- sum(terr_high$has_complete_data)
  n_complete_low <- sum(terr_low$has_complete_data)
  
  # Total species counts for each territoriality type
  n_total_high <- nrow(terr_high)
  n_total_low <- nrow(terr_low)
  
  # Calculate current proportions
  prop_high <- n_complete_high / n_total_high
  prop_low <- n_complete_low / n_total_low
  
  # If target proportion is not specified, use the lower of the two
  if (is.null(target_proportion)) {
    target_proportion <- min(prop_high, prop_low)
  }
  
  # Calculate how many species need to be removed from the group with higher proportion
  if (prop_high > prop_low) {
    # Need to remove species from the high territoriality group
    target_n_complete_high <- target_proportion * n_total_high
    n_to_remove <- n_complete_high - target_n_complete_high
  } else if (prop_low > prop_high) {
    # Need to remove species from the low territoriality group
    target_n_complete_low <- target_proportion * n_total_low
    n_to_remove <- n_complete_low - target_n_complete_low
  } else {
    # Proportions are already equal
    n_to_remove <- 0
  }
  
  samplingMatrixOut <- create_downsampling_matrix(df = df, territoriality_col = territoriality_col, territory_value_high = territory_value_high, territory_value_low = territory_value_low, data_cols = data_cols)
  
  # Create a detailed report
  report <- list(
    proportions = data.frame(
      territoriality = c("High", "Low"),
      species_with_data = c(n_complete_high, n_complete_low),
      total_species = c(n_total_high, n_total_low),
      proportion = c(prop_high, prop_low)
    ),
    target_proportion = target_proportion,
    downsample_info = if (prop_high > prop_low) {
      list(
        group_to_downsample = "High territoriality",
        territory_value = territory_value_high,
        n_to_remove = round(n_to_remove),
        original_proportion = prop_high,
        new_proportion = target_proportion
      )
    } else if (prop_low > prop_high) {
      list(
        group_to_downsample = "Low territoriality",
        territory_value = territory_value_low, 
        n_to_remove = round(n_to_remove),
        original_proportion = prop_low,
        new_proportion = target_proportion
      )
    } else {
      list(
        group_to_downsample = "None - already balanced",
        n_to_remove = 0
      )
    }
  )
  
  # Return the number to remove and the detailed report
  return(list(
    n_to_remove = ceiling(n_to_remove),
    report = report,
    filter_statement = if (n_to_remove > 0) {
      if (prop_high > prop_low) {
        paste0("To balance, remove ", ceiling(n_to_remove), " species with ", 
               territoriality_col, " = ", territory_value_high, 
               " and complete data for ", paste(data_cols, collapse = " and "))
      } else {
        paste0("To balance, remove ", ceiling(n_to_remove), " species with ", 
               territoriality_col, " = ", territory_value_low, 
               " and complete data for ", paste(data_cols, collapse = " and "))
      }
    } else {
      "No downsampling needed - proportions are already balanced"
    }
  ))
}

# Function to create a downsampling matrix table (similar to the original code)
create_downsampling_matrix <- function(df, 
                                       territoriality_col,
                                       territory_value_high = 1, 
                                       territory_value_low = 0,
                                       data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")) {
  
  # Create a variable for having complete data
  df$has_complete_data <- !apply(df[, data_cols, drop = FALSE], 1, function(row) any(is.na(row)))
  
  df[,territoriality_col] <- as.character(df[,territoriality_col])
  
  # Filter to rows with non-NA territory values
  df_filtered <- df[!is.na(df[[territoriality_col]]), ]
  
  # Create a matrix similar to the original code
  matrix <- df_filtered %>%
    group_by(.data[[territoriality_col]], has_complete_data) %>%
    summarize(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = has_complete_data,
      values_from = n,
      names_prefix = "has_data_"
    ) %>%
    mutate(Total = has_data_FALSE + has_data_TRUE) %>%
    bind_rows(
      summarize(.,
                !!territoriality_col := "Total",
                has_data_FALSE = sum(has_data_FALSE),
                has_data_TRUE = sum(has_data_TRUE),
                Total = sum(Total))
    ) %>%
    mutate(
      pct_has_data_FALSE = has_data_FALSE/Total, 
      pct_has_data_TRUE = has_data_TRUE/Total
    )
  
  # Get totals for percentage calculations
  total_row <- matrix %>% filter(.data[[territoriality_col]] == "Total")
  
  # Add percentage rows 
  pct_rows <- matrix %>%
    filter(.data[[territoriality_col]] %in% c(territory_value_high, territory_value_low)) %>%
    mutate(
      has_data_FALSE = has_data_FALSE / total_row$has_data_FALSE,
      has_data_TRUE = has_data_TRUE / total_row$has_data_TRUE,
      Total = Total / total_row$Total,
      !!territoriality_col := paste0("pct_", territoriality_col, "_", .data[[territoriality_col]]),
      pct_has_data_FALSE = NA,
      pct_has_data_TRUE = NA
    )
  
  # Combine everything
  final_matrix <- bind_rows(matrix, pct_rows)
  print(final_matrix)
  
  return(final_matrix)
}

#Example usage:

# Calculate downsampling amount
result <- calculate_downsampling(
  df = dfIn,
  territoriality_col = "Territory_12vs3",
  territory_value_high = 1, # Territory=3 (year-round)
  territory_value_low = 0,  # Territory=1or2 (seasonal)
  data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")
)

# Print results
print(result$n_to_remove)
print(result$filter_statement)
print(result$report$proportions)

# Create the matrix for visualization
matrix <- create_downsampling_matrix(
  df = dfIn,
  territoriality_col = "Territory_12vs3"
)
print(matrix)

#### Analyses ----
# Binary predictors
predictors <- c(
  "HighConfidence_Coop", "Griesser2017FamilialLiving", "colonial_Griesser2023",
  "Final.polygyny", "Territory_12vs3", "TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "PermissiveExclusive", "GeographicRegion"
)

# Run logistic regression and collect summaries
regression_results <- map_dfr(predictors, function(var) {
  model <- glm(as.formula(paste0("FemaleSong_Agg01 ~ ", var)), data = df_FS, family = binomial)
  tidy(model) %>%
    mutate(Predictor = var)
})

# Format nicely and save
regression_results %>%
  select(Predictor, term, estimate, std.error, statistic, p.value) %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  write.csv("All_LogisticRegressions_FemaleSong_Agg01.csv", row.names = FALSE)



## Summary by geographic region
library(dplyr)

# Binary and continuous trait columns
binary_cols <- c("HaveFSData", "HaveCBData", "FemaleSong_Agg01", "HighConfidence_Coop",
                 "Griesser2017FamilialLiving", "colonial_Griesser2023", "Final.polygyny",
                 "Territory_12vs3", "TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "PermissiveExclusive")

continuous_cols <- c("MaleFemalePlumageDiffAbs", "sumDiffPC1_Dunn2015", "sumDiffPC2_Dunn2015")

# Generate summary table by GeographicRegion
summary_table <- dfOs %>%
  group_by(GeographicRegion) %>%
  summarise(
    N_species = n(),
    
    # % WITH DATA columns
    FS_data_pct     = mean(!is.na(FemaleSong_Agg01)) * 100,
    CB_data_pct     = mean(!is.na(HighConfidence_Coop)) * 100,
    
    # % TRAIT PRESENT among those with data
    FS_present_pct  = mean(FemaleSong_Agg01 == 1, na.rm = TRUE) * 100,
    CB_present_pct  = mean(HighConfidence_Coop == 1, na.rm = TRUE) * 100,
    
    FamLiving_pct   = mean(Griesser2017FamilialLiving == 1, na.rm = TRUE) * 100,
    Colonial_pct    = mean(colonial_Griesser2023 == "colonial", na.rm = TRUE) * 100,
    Polygynous_pct  = mean(Final.polygyny == 1, na.rm = TRUE) * 100,
    YearRoundTerr_pct = mean(Territory_12vs3 == 1, na.rm = TRUE) * 100,
    
    # Continuous trait means
    #PlumageDiffAbs_mean = mean(MaleFemalePlumageDiffAbs, na.rm = TRUE),
    PC1Diff_mean         = mean(sumDiffPC1_Dunn2015, na.rm = TRUE),
    PC2Diff_mean         = mean(sumDiffPC2_Dunn2015, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))  # round numeric columns

# View the result
print(summary_table)
View(summary_table)

# Optionally export to CSV
write.csv(summary_table, paste(Sys.Date(),"Summary_By_GeographicRegion.csv"), row.names = FALSE)


## full stratefied summary

library(dplyr)

# Binary columns to cross with GeographicRegion
binary_traits <- c("HighConfidence_Coop", "Griesser2017FamilialLiving", 
                   "colonial_Griesser2023", "Final.polygyny", "Territory_12vs3", "TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "TerritorialityPermissiveExclusive", "TerritorialityPermissiveExclusiveHighConf")

# Initialize empty list to hold each result
all_summaries <- list()

# Loop through each trait to compute stratified summaries
for (trait in binary_traits) {
  df_temp <- dfOs %>%
    filter(!is.na(.data[[trait]]), !is.na(GeographicRegion)) %>%
    group_by(GeographicRegion, !!sym(trait)) %>%
    summarise(
      N_species = n(),
      HaveFSData_N = sum(HaveFSData),
      HaveFSData_pct = mean(HaveFSData, na.rm = TRUE) * 100,
      FS_present_N = sum(FemaleSong_Agg01, na.rm = TRUE),
      FS_present_pct = mean(FemaleSong_Agg01 == 1, na.rm = TRUE) * 100,
      HaveCBData_N = sum(HaveCBData),
      HaveCBData_pct = mean(HaveCBData, na.rm = TRUE) * 100,
      CB_present_N = sum(HighConfidence_Coop, na.rm = TRUE),
      CB_present_pct = mean(HighConfidence_Coop == 1, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(Trait = trait, Trait_Value = as.character(!!sym(trait))) %>%
    select(Trait, Trait_Value, GeographicRegion, N_species, HaveFSData_N, HaveFSData_pct, FS_present_N, FS_present_pct, HaveCBData_N, HaveCBData_pct, CB_present_N, CB_present_pct)
  
  all_summaries[[trait]] <- df_temp
}

# Combine all results
final_summary <- bind_rows(all_summaries)

# Round numeric columns
final_summary <- final_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# View or export
print(final_summary)
write.csv(final_summary, paste(Sys.Date(),"FS_CB_Presence_StratifiedSummary.csv"), row.names = FALSE)

#### phylANOVA cooperative breeding and body mass ----
newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-09-2.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')
allcols = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/2024-8-4AllCol_Duet_Female_AllTraits_fromJiaying2025-03-04.csv")
longevitydata = read.csv("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/justlongevitydata.csv")

# dfIn1 = read.csv(newdata)
# dfIn = merge(dfIn1, longevitydata, by = "species", all.x = T)
# dflongsubset = merge(dfIn1, longevitydata, by = "species")

tree = read.nexus(treefile)

df = read.csv(newdata)
rownames(df) <- df$species
df = df[which(df$species %in% tree$tip.label),]
df$logMass_AVONET = log(df$Mass_AVONET)
df$absCentroid.Latitude.AVONET = abs(df$Centroid.Latitude_AVONET)

dflongsubset = merge(df, longevitydata, by = "species")


layout(mat = matrix(1:6,2,3, byrow = T))

dfsub <- df[complete.cases(df[,c("HighConfidence_Coop", "logMass_AVONET")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$logMass_AVONET
discvec = dfsub$HighConfidence_Coop
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$logMass_AVONET ~ dfsub$HighConfidence_Coop, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))

dfsub <- df[complete.cases(df[,c("HighConfidence_Coop", "absCentroid.Latitude.AVONET")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$absCentroid.Latitude.AVONET
discvec = dfsub$HighConfidence_Coop
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$absCentroid.Latitude.AVONET ~ dfsub$HighConfidence_Coop, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))

dfsub <- dflongsubset[complete.cases(dflongsubset[,c("HighConfidence_Coop", "log_Longevity")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$log_Longevity
discvec = dfsub$HighConfidence_Coop
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$log_Longevity ~ dfsub$HighConfidence_Coop, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))

dfsub <- df[complete.cases(df[,c("FemaleSong_Agg01", "logMass_AVONET")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$logMass_AVONET
discvec = dfsub$FemaleSong_Agg01
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$logMass_AVONET ~ dfsub$FemaleSong_Agg01, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))

dfsub <- df[complete.cases(df[,c("FemaleSong_Agg01", "absCentroid.Latitude.AVONET")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$absCentroid.Latitude.AVONET
discvec = dfsub$FemaleSong_Agg01
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$absCentroid.Latitude.AVONET ~ dfsub$FemaleSong_Agg01, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))


dfsub <- dflongsubset[complete.cases(dflongsubset[,c("FemaleSong_Agg01", "log_Longevity")]),]
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfsub$species))
contvec = dfsub$log_Longevity
discvec = dfsub$FemaleSong_Agg01
names(contvec) <- names(discvec) <- dfsub$species
phylanovaout <- phylANOVA(tree_clean, x = discvec, y = contvec)
boxplot(dfsub$log_Longevity ~ dfsub$FemaleSong_Agg01, main = paste("phylanova p =", phylanovaout$Pf, "\nNspecies =", length(contvec)))

layout(matrix(1:3,1,3, byrow = T))
plot(dflongsubset$Max_longevity, dflongsubset$Mass_AVONET)
plot(dflongsubset$log_Longevity, dflongsubset$logMass_AVONET)

plot(dflongsubset$Max_longevity, dflongsubset$logMass_AVONET)
