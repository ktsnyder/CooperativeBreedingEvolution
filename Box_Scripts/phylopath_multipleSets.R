# multi-stage model testing - phylopath
# NO LONGER USED
# However, this is the file where "Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04.csv" was created (also no longer using)

library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(phytools)
library(phylopath)

boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

jiayingdata = read.csv('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/2024-8-4AllCol_Duet_Female_AllTraits_fromJiaying2025-03-04.csv')
passertree = read.nexus('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex')
sum(jiayingdata$BirdtreeSpecies %in% passertree$tip.label)


# Read your data
bird_data <- dfIn <- read.csv(newdata)
bird_data <- dfIn
bird_tree = read.nexus(treefile)

jiayingdata = jiayingdata[which(!is.na(jiayingdata$NF_Solo_Duet)), c("BirdtreeSpecies", "Family3_BirdtreeMatchSpecies2_AVONET", "Order_AVONET", "NF_Solo_Duet")]
jiayingdata[!jiayingdata$BirdtreeSpecies %in% bird_data$species, "Family3_BirdtreeMatchSpecies2_AVONET"]
sum(!jiayingdata$BirdtreeSpecies[which(!is.na(jiayingdata$NF_Solo_Duet))] %in% bird_data$species)
jiayingdata[which(!jiayingdata$BirdtreeSpecies %in% bird_data$species),] %>% group_by(NF_Solo_Duet) %>% count
jiayingdata$JiayingSolo = ifelse(jiayingdata$NF_Solo_Duet == 1, "present", "absent")
jiayingdata$JiayingDuet = ifelse(jiayingdata$NF_Solo_Duet == 2, "present", "absent")
jiayingdata$JiayingSolo[which(is.na(jiayingdata$NF_Solo_Duet))] = NA
jiayingdata$JiayingDuet[which(is.na(jiayingdata$NF_Solo_Duet))] = NA
bird_data = merge(bird_data, jiayingdata[,c("BirdtreeSpecies", "JiayingSolo", "JiayingDuet", "NF_Solo_Duet", "Order_AVONET")], by.x = "species", by.y = "BirdtreeSpecies", all = T, suffixes = c("", "_J"))
#write.csv(bird_data, "Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04.csv", row.names = F)



## make binary 0/1 data into character factors
binary_vars = c("FemaleSong_Agg01", "HighConfidence_Coop", "Duet", "Chorus", "Final.polygyny", "Final.EPP", "Griesser2023.Colonial01", "Griesser2017FamilialLiving", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.LongSocialBonds", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.GroupsLargerThanPair", "Griesser2023.Asocial0VsSocial1", "Griesser2023.LargestGroupSizes", "Communal", "Territory_12vs3", "Territory_1vs23")
for(var in binary_vars) {
  # Convert to 0/1 numeric first
  bird_data[[var]] <- as.numeric(as.character(bird_data[[var]]))
  
  # Then explicitly convert to binary factor with labels
  bird_data[[var]] <- factor(bird_data[[var]], levels = c(0, 1), labels = c("absent", "present"))
  
  # Print confirmation
  cat("Converted", var, "to binary factor with levels:", paste(levels(bird_data[[var]]), collapse=", "), "\n")
}

## make other binary columns into factors
binary_non_integer = c("sedentariness_Griesser2023", "insularity_Griesser2023")
bird_data[["sedentariness_Griesser2023"]] <- as.factor(bird_data[["sedentariness_Griesser2023"]])
bird_data[["insularity_Griesser2023"]] <- as.factor(bird_data[["insularity_Griesser2023"]])
bird_data[["JiayingSolo"]] <- as.factor(bird_data[["JiayingSolo"]])
bird_data[["JiayingDuet"]] <- as.factor(bird_data[["JiayingDuet"]])



## Data columns that are ordered 
can_be_made_ordered_cat = c("social_bonds_Griesser2023", "grouping_Griesser2023", "Trophic.Level_AVONET")
# Convert Trophic.Level_AVONET to ordered integer
bird_data$Trophic.Level_AVONET_int <- NA
bird_data$Trophic.Level_AVONET_int[bird_data$Trophic.Level_AVONET == "Herbivore"] <- 1
bird_data$Trophic.Level_AVONET_int[bird_data$Trophic.Level_AVONET == "Omnivore"] <- 2
bird_data$Trophic.Level_AVONET_int[bird_data$Trophic.Level_AVONET == "Carnivore"] <- 3

# Convert social_bonds_Griesser2023 to ordered integer
bird_data$social_bonds_Griesser2023_int <- NA
bird_data$social_bonds_Griesser2023_int[bird_data$social_bonds_Griesser2023 == "a-short"] <- 1
bird_data$social_bonds_Griesser2023_int[bird_data$social_bonds_Griesser2023 == "b-season"] <- 2
bird_data$social_bonds_Griesser2023_int[bird_data$social_bonds_Griesser2023 == "c-long"] <- 3

# Convert grouping_Griesser2023 to ordered integer
bird_data$grouping_Griesser2023_int <- NA
bird_data$grouping_Griesser2023_int[bird_data$grouping_Griesser2023 == "asocial"] <- 1
bird_data$grouping_Griesser2023_int[bird_data$grouping_Griesser2023 == "pair"] <- 2
bird_data$grouping_Griesser2023_int[bird_data$grouping_Griesser2023 == "small_groups"] <- 3
bird_data$grouping_Griesser2023_int[bird_data$grouping_Griesser2023 == "large_groups"] <- 4


# 2. Scale all ordered categorical variables - directly without any complex logic
ordered_cat_vars <- c("Territory", "Habitat.Density_AVONET", "Migration_AVONET", 
                      "Trophic.Level_AVONET_int", "social_bonds_Griesser2023_int", 
                      "grouping_Griesser2023_int")

for(level_var in ordered_cat_vars) {
  if(sum(!is.na(bird_data[,level_var])) > 0) {
    bird_data[,paste0(level_var,"_scaled")] <- scale(bird_data[,level_var])
    cat("Created scaled variable:", paste0(level_var,"_scaled"), "\n")
  } else {
    warning(paste("Variable", level_var, "contains only NAs and cannot be scaled"))
  }
}
  
## scale and center continuous data
continuous_vars = c("Percentage.of.extra.group.paternity_Cornwallis2017", "EPP1_Biagolini2017", "EPB1_Biagolini2017", "egg_mass_Griesser2023", "weight_Griesser2023", "brain_Griesser2023", "time_fed_Griesser2023", "clutch_size_Griesser2023", "caretakers_Griesser2023", "longevity_Griesser2023", "food_energy_Griesser2023", "food_h_level_Griesser2023", "fibres_Griesser2023", "Syllable.rep.final", "Song.rep.final", "Duration.final", "Beak.Length_Culmen_AVONET", "Beak.Length_Nares_AVONET", "Beak.Width_AVONET", "Beak.Depth_AVONET", "Tarsus.Length_AVONET", "Wing.Length_AVONET", "Mass_AVONET", "Min.Latitude_AVONET", "Max.Latitude_AVONET", "Centroid.Latitude_AVONET", "Centroid.Longitude_AVONET", "Range.Size_AVONET")
library(car)

# Function to check normality and apply appropriate transformation
process_continuous_var <- function(data, var_name) {
  # Extract non-NA values
  valid_values <- data[!is.na(data[,var_name]), var_name]
  
  if(length(valid_values) < 3) {
    warning(paste("Variable", var_name, "has too few non-NA values for normality testing"))
    return("insufficient_data")
  }
  
  # Create scaled version (always)
  data[,paste0(var_name, "_scaled")] <- scale(data[,var_name])
  
  # Test for normality using Shapiro-Wilk test (for samples < 5000)
  if(length(valid_values) < 5000) {
    tryCatch({
      normality_test <- shapiro.test(valid_values)
      is_normal <- normality_test$p.value > 0.05
    }, error = function(e) {
      warning(paste("Error in normality test for", var_name, ":", e$message))
      is_normal <- FALSE
    })
  } else {
    # For larger samples, use QQ plot assessment via Kolmogorov-Smirnov test
    tryCatch({
      ks_test <- ks.test(scale(valid_values), "pnorm")
      is_normal <- ks_test$p.value > 0.05
    }, error = function(e) {
      warning(paste("Error in KS test for", var_name, ":", e$message))
      is_normal <- FALSE
    })
  }
  
  # Special handling for latitude/longitude values which can be negative
  if(grepl("Latitude|Longitude", var_name)) {
    # For geographic coordinates, scaling is usually sufficient
    return("scaled_only")
  }
  
  # Check if log transformation improves normality
  if(!is_normal && min(valid_values, na.rm=TRUE) > 0) {
    # For strictly positive values, use natural log
    log_values <- log(valid_values)
    
    tryCatch({
      if(length(log_values) < 5000) {
        log_normality_test <- shapiro.test(log_values)
        log_is_better <- log_normality_test$p.value > normality_test$p.value
      } else {
        log_ks_test <- ks.test(scale(log_values), "pnorm")
        log_is_better <- log_ks_test$p.value > ks_test$p.value
      }
      
      if(log_is_better) {
        data[,paste0(var_name, "_log")] <- log(data[,var_name])
        data[,paste0(var_name, "_log_scaled")] <- scale(log(data[,var_name]))
        return("log_transformed")
      }
    }, error = function(e) {
      warning(paste("Error in log normality test for", var_name, ":", e$message))
    })
  } else if(!is_normal && min(valid_values, na.rm=TRUE) >= 0) {
    # For values including zeros, use log1p
    log1p_values <- log1p(valid_values)
    
    tryCatch({
      if(length(log1p_values) < 5000) {
        log1p_normality_test <- shapiro.test(log1p_values)
        log1p_is_better <- log1p_normality_test$p.value > normality_test$p.value
      } else {
        log1p_ks_test <- ks.test(scale(log1p_values), "pnorm")
        log1p_is_better <- log1p_ks_test$p.value > ks_test$p.value
      }
      
      if(log1p_is_better) {
        data[,paste0(var_name, "_log1p")] <- log1p(data[,var_name])
        data[,paste0(var_name, "_log1p_scaled")] <- scale(log1p(data[,var_name]))
        return("log1p_transformed")
      }
    }, error = function(e) {
      warning(paste("Error in log1p normality test for", var_name, ":", e$message))
    })
  }
  
  # If we reach here, no transformation improved normality, or there was an error
  return("scaled_only")
}

# scale Centroid.Latitude_AVONET_scaled without normality testing
bird_data$Centroid.Latitude_AVONET_scaled <- scale(bird_data$Centroid.Latitude_AVONET)

# Process each continuous variable, with improved error handling
transformation_results <- data.frame(variable = continuous_vars, 
                                     transformation = NA,
                                     stringsAsFactors = FALSE)
for(i in 1:length(continuous_vars)) {
  var_name <- continuous_vars[i]
  cat("Processing variable:", var_name, "...\n")
  
  # Skip the variable if it doesn't exist in the data
  if(!var_name %in% colnames(bird_data)) {
    warning(paste("Variable", var_name, "not found in dataset"))
    transformation_results$transformation[i] <- "not_found"
    next
  }
  
  # Skip latitude/longitude variables - handle these separately
  if(var_name %in% c("Min.Latitude_AVONET", "Max.Latitude_AVONET", 
                     "Centroid.Latitude_AVONET", "Centroid.Longitude_AVONET")) {
    # Just scale these variables
    bird_data[,paste0(var_name, "_scaled")] <- scale(bird_data[,var_name])
    transformation_results$transformation[i] <- "scaled_only"
    next
  }
  
  # Apply transformations to other variables
  result <- tryCatch({
    process_continuous_var(bird_data, var_name)
  }, error = function(e) {
    warning(paste("Error processing variable", var_name, ":", e$message))
    # Still create a scaled version even if transformation analysis fails
    bird_data[,paste0(var_name, "_scaled")] <- scale(bird_data[,var_name])
    return("scaled_only_with_error")
  })
  
  transformation_results$transformation[i] <- result
}



# Prepare core dataset
#core_vars <- c("FemaleSong_Agg01", "HighConfidence_Coop", "Territory", "Social.bond", "JiayingDuet", "JiayingSolo")
core_vars <- c("HighConfidence_Coop", "Territory", "Social.bond", "JiayingDuet", "JiayingSolo")
core_data <- bird_data %>% 
  select(all_of(core_vars), species) %>% 
  filter(complete.cases(.))
rownames(core_data) <- core_data$species
core_data <- core_data %>% select(-species)

# Define territory variants
core_data$Territory1 <- as.factor(ifelse(core_data$Territory == 1, 1, 0))
core_data$Territory2 <- as.factor(ifelse(core_data$Territory == 2, 1, 0))
core_data$Territory3 <- as.factor(ifelse(core_data$Territory == 3, 1, 0))
core_data$Territory_12vs3 <- as.factor(ifelse(core_data$Territory == 3, 1, 0))
core_data$Territory_1vs23 <- as.factor(ifelse(core_data$Territory == 1, 0, 1))

binary_terr_vars = c("Territory1", "Territory2", "Territory3", "Territory_12vs3", "Territory_1vs23")
for(var in binary_terr_vars) {
  # Convert to 0/1 numeric first
#  core_data[[var]] <- as.numeric(as.character(core_data[[var]]))
  
  # Then explicitly convert to binary factor with labels
  core_data[[var]] <- factor(core_data[[var]], levels = c(0, 1), labels = c("absent", "present"))
  
  # Print confirmation
  cat("Converted", var, "to binary factor with levels:", paste(levels(core_data[[var]]), collapse=", "), "\n")
}

# Match tree to data
core_tree <- ape::drop.tip(bird_tree, setdiff(bird_tree$tip.label, rownames(core_data)))

# Define core models
# Models exploring female song evolution
# Models with Territory1, Territory2, Territory3
fs_models_individual <- define_model_set(
 # fs_terr1 = c(FemaleSong_Agg01 ~ Territory1),
  fs_terr2 = c(FemaleSong_Agg01 ~ Territory2),
  fs_terr3 = c(FemaleSong_Agg01 ~ Territory3),
  fs_coop = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  fs_social = c(FemaleSong_Agg01 ~ Social.bond),
  fs_duet = c(FemaleSong_Agg01 ~ Duet),
  fs_chorus = c(FemaleSong_Agg01 ~ Chorus),
 # fs_terr1_coop = c(FemaleSong_Agg01 ~ Territory1 + HighConfidence_Coop),
  fs_terr2_coop = c(FemaleSong_Agg01 ~ Territory2 + HighConfidence_Coop),
  fs_terr3_coop = c(FemaleSong_Agg01 ~ Territory3 + HighConfidence_Coop),
  fs_terr2_social = c(FemaleSong_Agg01 ~ Territory2 + Social.bond),
  fs_terr3_social = c(FemaleSong_Agg01 ~ Territory3 + Social.bond),
  fs_terr3_social_chorus = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + Chorus),
  fs_terr3_social_duet = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + Duet),
  fs_terr3_social_coop = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + HighConfidence_Coop),
  fs_terr2_social_chorus = c(FemaleSong_Agg01 ~ Territory2 + Social.bond + Chorus),
  fs_terr2_social_duet = c(FemaleSong_Agg01 ~ Territory2 + Social.bond + Duet),
  fs_terr2_social_coop = c(FemaleSong_Agg01 ~ Territory2 + Social.bond + HighConfidence_Coop),
  fs_terr3_social_chorus_coop = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + Chorus + HighConfidence_Coop),
  fs_terr3_social_duet_coop = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + Duet + HighConfidence_Coop),
  fs_terr2_social_duet_coop = c(FemaleSong_Agg01 ~ Territory2 + Social.bond + Duet + HighConfidence_Coop),
  fs_terr2_social_chorus_coop = c(FemaleSong_Agg01 ~ Territory2 + Social.bond + Chorus + HighConfidence_Coop),
  # duet_fs = c(Duet ~ FemaleSong_Agg01),
  # duet_fs_terr2_social = c(Duet ~ FemaleSong_Agg01 + Territory2 + Social.bond),
  # duet_fs_terr3_social = c(Duet ~ FemaleSong_Agg01 + Territory3 + Social.bond),
  # duet_fs_terr2_social_coop = c(Duet ~ FemaleSong_Agg01 + Territory2 + Social.bond + HighConfidence_Coop),
  # duet_fs_terr3_social_coop = c(Duet ~ FemaleSong_Agg01 + Territory3 + Social.bond + HighConfidence_Coop),
 # full_chain = c(FemaleSong_Agg01 ~ Duet, 
 #                Duet ~ Territory3 + Social.bond + HighConfidence_Coop),
 # 
 # partial_chain1 = c(FemaleSong_Agg01 ~ Duet, 
 #                    Duet ~ Territory3 + Social.bond),
 # 
 # partial_chain2 = c(FemaleSong_Agg01 ~ Duet, 
 #                    Duet ~ Social.bond + HighConfidence_Coop),
 # 
 # reverse_chain = c(Duet ~ FemaleSong_Agg01, 
 #                   FemaleSong_Agg01 ~ Territory3 + Social.bond + HighConfidence_Coop),
 # 
 # For comparison
 # direct_effects = c(FemaleSong_Agg01 ~ Territory3 + Social.bond + HighConfidence_Coop + Duet),
  .common = c()
)

# Models with binary terr grouping classifications
fs_models_binary <- define_model_set(
  #fs_terr12v3 = c(FemaleSong_Agg01 ~ Territory_12vs3),
  fs_terr1v23 = c(FemaleSong_Agg01 ~ Territory_1vs23),
  fs_coop = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  fs_social = c(FemaleSong_Agg01 ~ Social.bond),
  fs_duet = c(FemaleSong_Agg01 ~ Duet),
  fs_chorus = c(FemaleSong_Agg01 ~ Chorus),
  #fs_terr12v3_coop = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop),
  fs_terr1v23_coop = c(FemaleSong_Agg01 ~ Territory_1vs23 + HighConfidence_Coop),
  #fs_terr12v3_duet = c(FemaleSong_Agg01 ~ Territory_12vs3 + Duet),
  fs_terr1v23_duet = c(FemaleSong_Agg01 ~ Territory_1vs23 + Duet),
  #fs_terr12v3_chorus = c(FemaleSong_Agg01 ~ Territory_12vs3 + Chorus),
  fs_terr1v23_chorus = c(FemaleSong_Agg01 ~ Territory_1vs23 + Chorus),
  #fs_terr12v3_coop_duet = c(FemaleSong_Agg01 ~ Territory_12vs3 + HighConfidence_Coop + Duet),
  fs_terr1v23_coop_duet = c(FemaleSong_Agg01 ~ Territory_1vs23 + HighConfidence_Coop + Duet),
  fs_terr1v23_coop_chorus = c(FemaleSong_Agg01 ~ Territory_1vs23 + HighConfidence_Coop + Chorus),
  fs_terr1v23_duet_social = c(FemaleSong_Agg01 ~ Territory_1vs23 + Duet + Social.bond),
  fs_terr1v23_chorus_social = c(FemaleSong_Agg01 ~ Territory_1vs23 + Chorus + Social.bond),
  fs_terr1v23_coop_chorus_social = c(FemaleSong_Agg01 ~ Territory_1vs23 + HighConfidence_Coop + Chorus + Social.bond),
  fs_terr1v23_coop_duet_social = c(FemaleSong_Agg01 ~ Territory_1vs23 + HighConfidence_Coop + Duet + Social.bond),
  duet_fs = c(Duet ~ FemaleSong_Agg01),
  duet_fs_terr1v23 = c(Duet ~ FemaleSong_Agg01 + Territory_1vs23),
  duet_fs_terr1v23_social = c(Duet ~ FemaleSong_Agg01 + Territory_1vs23 + Social.bond),
  duet_fs_terr1v23_social_coop = c(Duet ~ FemaleSong_Agg01 + Territory_1vs23 + Social.bond + HighConfidence_Coop),
  .common = c()
)

fs_models_IntegerTerr <- define_model_set(
  fs_terrInteger = c(FemaleSong_Agg01 ~ Territory),
  fs_duet = c(FemaleSong_Agg01 ~ Duet),
  fs_chorus = c(FemaleSong_Agg01 ~ Chorus),
  fs_coop = c(FemaleSong_Agg01 ~ HighConfidence_Coop),
  fs_social = c(FemaleSong_Agg01 ~ Social.bond),
  fs_terrInteger_coop = c(FemaleSong_Agg01 ~ Territory + HighConfidence_Coop),
  fs_terrInteger_duet = c(FemaleSong_Agg01 ~ Territory + Duet),
  fs_terrInteger_chorus = c(FemaleSong_Agg01 ~ Territory + Chorus),
  fs_terrInteger_social = c(FemaleSong_Agg01 ~ Territory + Social.bond),
  fs_chorus_coop = c(FemaleSong_Agg01 ~ Chorus + HighConfidence_Coop),
  fs_duet_coop = c(FemaleSong_Agg01 ~ Duet + HighConfidence_Coop),
  fs_social_coop = c(FemaleSong_Agg01 ~ Social.bond + HighConfidence_Coop),
  fs_terrInteger_social_coop = c(FemaleSong_Agg01 ~ Territory + Social.bond + HighConfidence_Coop),
  fs_terrInteger_duet_coop = c(FemaleSong_Agg01 ~ Territory + Duet + HighConfidence_Coop),
  fs_terrInteger_chorus_coop = c(FemaleSong_Agg01 ~ Territory + Chorus + HighConfidence_Coop),
  fs_terrInteger_duet_coop_social = c(FemaleSong_Agg01 ~ Territory + Duet + HighConfidence_Coop + Social.bond),
  fs_terrInteger_chorus_coop_social = c(FemaleSong_Agg01 ~ Territory + Chorus + HighConfidence_Coop + Social.bond),
  duet_fs = c(Duet ~ FemaleSong_Agg01),
  duet_fs_terr = c(Duet ~ FemaleSong_Agg01 + Territory),
  duet_fs_terr_social = c(Duet ~ FemaleSong_Agg01 + Territory + Social.bond),
  duet_fs_terr_social_coop = c(Duet ~ FemaleSong_Agg01 + Territory + Social.bond + HighConfidence_Coop),
  # full_chain = c(FemaleSong_Agg01 ~ Duet, 
  #                Duet ~ Territory + Social.bond + HighConfidence_Coop),
  # 
  #partial_chain1 = c(FemaleSong_Agg01 ~ Duet, 
  #                   Duet ~ Territory + Social.bond),
  # 
  # partial_chain = c(FemaleSong_Agg01 ~ Duet, 
  #                    Duet ~ Social.bond + HighConfidence_Coop),
  # 
  # reverse_chain = c(Duet ~ FemaleSong_Agg01, 
  #                   FemaleSong_Agg01 ~ Territory + Social.bond + HighConfidence_Coop),
  # 
  # # For comparison
  # direct_effects = c(FemaleSong_Agg01 ~ Territory + Social.bond + HighConfidence_Coop + Duet),
  .common = c()
)

# # Run analysis with binary classifications
# fs_results_bin_BM <- phylo_path(fs_models_binary, core_data, core_tree, model = "BM")
# summary_bin_BM <- summary(fs_results_bin_BM)

## model = lambda
# Run analysis with individual territory levels
fs_results_ind_lambda <- phylo_path(fs_models_individual, core_data, core_tree, model = "lambda")
(summary_ind_lambda <- summary(fs_results_ind_lambda))
plot(best(fs_results_ind_lambda))

# Run analysis with integer territory levels
fs_results_integerTerr_lambda <- phylo_path(fs_models_IntegerTerr, core_data, core_tree, model = "lambda")
(summary_int_lambda <- summary(fs_results_integerTerr_lambda))

## model = BM
# Run analysis with individual territory levels
fs_results_ind_BM <- phylo_path(fs_models_individual, core_data, core_tree, model = "BM")
summary_ind_BM <- summary(fs_results_ind_BM)

# Run analysis with integer territory levels
fs_results_integerTerr_BM <- phylo_path(fs_models_IntegerTerr, core_data, core_tree, model = "BM")
summary_int_BM <- summary(fs_results_integerTerr_BM)

## model = OU
# Run analysis with individual territory levels
fs_results_ind_OUrr <- phylo_path(fs_models_individual, core_data, core_tree, model = "OUrandomRoot")
summary_ind_OUrr <- summary(fs_results_ind_OUrr)

# Run analysis with integer territory levels
fs_results_integerTerr_OUrr <- phylo_path(fs_models_IntegerTerr, core_data, core_tree, model = "OUrandomRoot")
summary_int_OUrr <- summary(fs_results_integerTerr_OUrr)

## model = EB
# Run analysis with individual territory levels
fs_results_ind_EB <- phylo_path(fs_models_individual, core_data, core_tree, model = "EB")
summary_ind_EB <- summary(fs_results_ind_EB)

# Run analysis with integer territory levels
fs_results_integerTerr_EB <- phylo_path(fs_models_IntegerTerr, core_data, core_tree, model = "EB")
summary_int_EB <- summary(fs_results_integerTerr_EB)

summary_ind_lambda
summary_int_lambda
summary_ind_BM
summary_int_BM
summary_ind_OUrr
summary_int_OUrr
summary_ind_EB
summary_int_EB

lambdaIndPlot<- plot(best(fs_results_ind_lambda), text_size = 3)
lambdaIndPlot + labs(title = "lambda individualTerr - best")
BMIndPlot <- plot(best(fs_results_ind_BM), text_size = 3)
BMIndPlot + labs(title = "BM individualTerr - best")
OUrrIndPlot <- plot(best(fs_results_ind_OUrr), text_size = 3)
OUrrIndPlot + labs(title = "OUrr individualTerr - best")
EBIndPlot <- plot(best(fs_results_ind_EB), text_size = 3)
EBIndPlot + labs(title = "EB individualTerr - best")

lambdaIntPlot<- plot(best(fs_results_integerTerr_lambda), text_size = 3)
lambdaIntPlot + labs(title = "lambda integerTerr - best")
BMIntPlot <- plot(best(fs_results_integerTerr_BM), text_size = 3)
BMIntPlot + labs(title = "BM integerTerr - best")
OUrrIntPlot <- plot(best(fs_results_integerTerr_OUrr), text_size = 3)
OUrrIntPlot + labs(title = "OUrr integerTerr - best")
EBIntPlot <- plot(best(fs_results_integerTerr_EB), text_size = 3)
EBIntPlot + labs(title = "EB integerTerr - best")

# Compare best models from each analysis by their CICc values
# (Lower CICc indicates better fit)
best_ind_CICc <- min(summary_ind$CICc)
best_bin_CICc <- min(summary_bin$CICc)

print(paste("Best individual territory model CICc:", best_ind_CICc))
print(paste("Best binary territory model CICc:", best_bin_CICc))

# Models exploring cooperative breeding evolution
coop_models <- define_model_set(
  # Territory models
  coop_terr1 = c(HighConfidence_Coop ~ Territory1),
  coop_terr3 = c(HighConfidence_Coop ~ Territory3),
  coop_terr_strong = c(HighConfidence_Coop ~ Territory_12vs3),
  
  # Female song and social models
  coop_fs = c(HighConfidence_Coop ~ FemaleSong_Agg01),
  coop_social = c(HighConfidence_Coop ~ Social.bond),
  coop_duet = c(HighConfidence_Coop ~ Duet),
  coop_chorus = c(HighConfidence_Coop ~ Chorus),
  
  # Combination models
  coop_terr_social = c(HighConfidence_Coop ~ Territory3 + Social.bond),
  coop_terr_fs = c(HighConfidence_Coop ~ Territory3 + FemaleSong_Agg01),
  coop_duet_chorus = c(HighConfidence_Coop ~ Duet + Chorus),
  
  .common = c()
)

#### jiaying's solo and duet

fs_models_IntegerTerr <- define_model_set(
  Fsolo_terrInteger = c(JiayingSolo ~ Territory),
  #Fsolo_duet = c(JiayingSolo ~ JiayingDuet),
  Fsolo_coop = c(JiayingSolo ~ HighConfidence_Coop),
  Fsolo_social = c(JiayingSolo ~ Social.bond),
  Fsolo_terrInteger_coop = c(JiayingSolo ~ Territory + HighConfidence_Coop),
  #Fsolo_terrInteger_duet = c(JiayingSolo ~ Territory + JiayingDuet),
  Fsolo_terrInteger_social = c(JiayingSolo ~ Territory + Social.bond),
  #Fsolo_duet_coop = c(JiayingSolo ~ JiayingDuet + HighConfidence_Coop),
  Fsolo_social_coop = c(JiayingSolo ~ Social.bond + HighConfidence_Coop),
  Fsolo_terrInteger_social_coop = c(JiayingSolo ~ Territory + Social.bond + HighConfidence_Coop),
  #Fsolo_terrInteger_duet_coop = c(JiayingSolo ~ Territory + JiayingDuet + HighConfidence_Coop),
  #Fsolo_terrInteger_duet_coop_social = c(JiayingSolo ~ Territory + JiayingDuet + HighConfidence_Coop + Social.bond),
  duet_coop = c(JiayingDuet ~ HighConfidence_Coop),
  duet_social = c(JiayingDuet ~ Social.bond),
  duet_coop_social = c(JiayingDuet ~ HighConfidence_Coop + Social.bond),
  duet_terr_coop = c(JiayingDuet ~ Territory + HighConfidence_Coop),
  duet_terr_coop_social = c(JiayingDuet ~ Territory + HighConfidence_Coop + Social.bond),
  duet_terr_social = c(JiayingDuet ~ Territory + Social.bond),
  duet_terr = c(JiayingDuet ~ Territory),
  #duet_Fsolo = c(JiayingDuet ~ JiayingSolo),
  #duet_Fsolo_terr = c(JiayingDuet ~ JiayingSolo + Territory),
  #duet_Fsolo_terr_social = c(JiayingDuet ~ JiayingSolo + Territory + Social.bond),
  #duet_Fsolo_terr_social_coop = c(JiayingDuet ~ JiayingSolo + Territory + Social.bond + HighConfidence_Coop),
  #duet_Fsolo_terr_coop = c(JiayingDuet ~ JiayingSolo + Territory + HighConfidence_Coop),
  .common = c()
)

fs_models_binTerr <- define_model_set(
  Fsolo_terrInteger = c(JiayingSolo ~ Territory3),
  #Fsolo_duet = c(JiayingSolo ~ JiayingDuet),
  Fsolo_coop = c(JiayingSolo ~ HighConfidence_Coop),
  Fsolo_social = c(JiayingSolo ~ Territory2),
  Fsolo_terrInteger_coop = c(JiayingSolo ~ Territory + HighConfidence_Coop),
  #Fsolo_terrInteger_duet = c(JiayingSolo ~ Territory + JiayingDuet),
  Fsolo_terrInteger_social = c(JiayingSolo ~ Territory + Territory2),
  #Fsolo_duet_coop = c(JiayingSolo ~ JiayingDuet + HighConfidence_Coop),
  Fsolo_social_coop = c(JiayingSolo ~ Territory2 + HighConfidence_Coop),
  Fsolo_terrInteger_social_coop = c(JiayingSolo ~ Territory + Territory2 + HighConfidence_Coop),
  #Fsolo_terrInteger_duet_coop = c(JiayingSolo ~ Territory + JiayingDuet + HighConfidence_Coop),
  #Fsolo_terrInteger_duet_coop_social = c(JiayingSolo ~ Territory + JiayingDuet + HighConfidence_Coop + Territory2),
  duet_coop = c(JiayingDuet ~ HighConfidence_Coop),
  duet_social = c(JiayingDuet ~ Territory2),
  duet_coop_social = c(JiayingDuet ~ HighConfidence_Coop + Territory2),
  duet_terr_coop = c(JiayingDuet ~ Territory + HighConfidence_Coop),
  duet_terr_coop_social = c(JiayingDuet ~ Territory + HighConfidence_Coop + Territory2),
  duet_terr_social = c(JiayingDuet ~ Territory + Territory2),
  duet_terr = c(JiayingDuet ~ Territory),
  #duet_Fsolo = c(JiayingDuet ~ JiayingSolo),
  #duet_Fsolo_terr = c(JiayingDuet ~ JiayingSolo + Territory),
  #duet_Fsolo_terr_social = c(JiayingDuet ~ JiayingSolo + Territory + Territory2),
  #duet_Fsolo_terr_social_coop = c(JiayingDuet ~ JiayingSolo + Territory + Territory2 + HighConfidence_Coop),
  #duet_Fsolo_terr_coop = c(JiayingDuet ~ JiayingSolo + Territory + HighConfidence_Coop),
  .common = c()
)

fs_results_integerTerr <- phylo_path(fs_models_IntegerTerr, core_data, core_tree, model = "lambda")
summary_int_EB <- summary(fs_results_integerTerr)
