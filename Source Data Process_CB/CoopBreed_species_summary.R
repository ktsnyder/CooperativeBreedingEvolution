## Generate summary df for species cooperative breeding data
## Coded by Aleyna Loughran-Pierce and Kate Snyder
## Created 10/9/2020
## Edited: 8/18/2021 by Kate Snyder
## 3/8/2022 - add nonkin/kin etc, added Odom Female Song to output, added BOW data, changed MeanCoop to be only > 0.5 (not >=)
## 3/9/2022 - output all Classification methods as different columns in the same table, remove coopClassMethod arg, remove extra rows
## 3/10/2022 - correct Jetz et al classification to reflect more liberal use of Cockburn et al data; add OC data
## 5/30/2023 - commented out last write csv since song data already added; removed Female Song from output since not included in input yet
## 6/14/2023 - added Mikula and Dale columns to output
## 6/20/2023 - included Dale in MeanCoop calculation
## 9/13/2023 - added Cornwallis et al data handling
## 10/25/2023 - just made subsetlabel NULL for non-subsetted data
## 2/23/2024 - edited to work with new output of merge_data_allcolumns.R
## 2/24/2024 - edited to do Griesser et al 2023 binarization
## 5/13/2024 - added AnyNoncoopEqualsNoncoop, require dplyr
## 6/6/2024 - fixed AnyNoncoopEqualsNoncoop (had assigned anything containing "0" as "1" *facepalm*); made new output column "JetzCoopInclCockburn" to enable single-source testing of Jetz data; integrated code from end of merge_data_allcolumns.R that added the HighConfidence_Coop column to the database and associated HighConfCoopFile = "2024-05-13_CoopClassesWSourceColumns_HighConfCoopColumn.csv" arg
## 6/7/2024 - AnyNoncoopEqualsNoncoop - made it be NA if meanclass is NA 

#setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")

#coopbreedfile = "2024-02-24_Aggregate_CBSource_Data_AllColumns.csv"
#coopbreedfile = "2024-05-13_Aggregate_CBSource_Data_AllColumns.csv"
#CoopBreed_species_summary(coopbreedfile = coopbreedfile)


CoopBreed_species_summary <- function(coopbreedfile, songfile = "SongData_R_Update.csv", HighConfCoopFile = "2024-05-13_CoopClassesWSourceColumns_HighConfCoopColumn.csv") {
  require(dplyr)
  ourdf <- read.csv(file = coopbreedfile)
  
  summarydf <- set.seed(10)  #this creates a blank dataframe for the for loop to fill in
  
  for (i in 1:length(ourdf$BirdtreeSpecies)) {
    tempRowOut <- set.seed(10) # this creates a blank new row for our summary dataframe
    tempRowIn <- ourdf[i,] 
    species <- tempRowIn$BirdtreeSpecies
    
    if (is.na(tempRowIn$MateSys5_Dunn2015)) { 
      Dunn <- NA
    } else if (tempRowIn$MateSys5_Dunn2015 == "Coop") {  
      Dunn <- 1
    } else {
      Dunn <- NA  
    }
    
    ### This one I did a little differently because there were so many different "values" in this column, and those that contain the character string "CB" will get a 1, and any that don't will get a 0. 
    require(stringr)  # checks whether the package is loaded, and loads it if not - i.e. the equivalent of using library(stringr); it does need to be installed first 
    if (is.na(tempRowIn$BS_Biagolini2017)) {
      Biagolini <- NA
    } else if (str_detect(tempRowIn$BS_Biagolini2017, "CB", negate = TRUE)) {
      Biagolini <- 0
    } else if (str_detect(tempRowIn$BS_Biagolini2017, "CB", negate = FALSE)) {
      Biagolini <- 1
    } else {
      Biagolini <- NA
    }
    
    
    if (is.na(tempRowIn$breeding.system_Downing2015)) { 
      DowningCoop <- NA
      DowningKinCoop <- NA
    } else if (tempRowIn$breeding.system_Downing2015 == "noncooperative") {
      DowningCoop <- 0
      DowningKinCoop <- NA
    } else if (tempRowIn$breeding.system_Downing2015 %in% c("kinCooperator")) {  
      DowningCoop <- 1
      DowningKinCoop <- "K"
    } else if (tempRowIn$breeding.system_Downing2015 %in% c("nonkinCooperator")) {  
      DowningCoop <- 1
      DowningKinCoop <- "NK"
    } else {
      DowningCoop <- NA  
      DowningKinCoop <- NA
    }
    
    
    if (is.na(tempRowIn$Social_System_RubensteinLovette2007)) {
      Rubenstein <- NA
    } else if (tempRowIn$Social_System_RubensteinLovette2007 == "Noncooperative") {
      Rubenstein <- 0
    } else if (tempRowIn$Social_System_RubensteinLovette2007 %in% c("Cooperative", " Cooperative")) {
      Rubenstein <- 1
    } else {
      Rubenstein <- NA
    }
    
    
    if (is.na(tempRowIn$KnownParentalCare_Cockburn2006)) { 
      Cockburn<- NA
    } else if (tempRowIn$KnownParentalCare_Cockburn2006 %in% c("Pair", "Female_only", "Suspected", "Male_only")) {
      Cockburn <- 0
    } else if (tempRowIn$KnownParentalCare_Cockburn2006 == "Cooperation") {  
      Cockburn <- 1
    } else {
      Cockburn<- NA  
    }
    
    if (is.na(tempRowIn$System_Jetz2011)) {
      Jetz <- JetzInclCockburn <- NA
    } else if (tempRowIn$System_Jetz2011 == "Non-cooperative") {
      Jetz <- JetzInclCockburn <- 0
    } else if (tempRowIn$System_Jetz2011 == "Cooperative") {
      Jetz <- JetzInclCockburn <- 1
    } else {
      Jetz <- JetzInclCockburn <- NA
    }
    
    JetzSource <- as.character(tempRowIn$Source_Jetz2011)
    JetzSources <- as.character(-17:-2)
    
    if (!is.na(Jetz) & !is.na(Cockburn)) {  
      if  (JetzSource == "-1") {    # if both Jetz and Cockburn have a CoopBreed classification, and Jetz's source is Cockburn, Jetz data treated as redundant
        Jetz <- NA
      } else if (JetzSource %in% JetzSources) {
        Jetz <- Jetz
      }
    }
    
    if (is.na(tempRowIn$Social_breeding_system_when_cooperative_Riehl2013)) {
      ReihlCoop <- NA
    } else {
      ReihlCoop <- 1  # because all entries in this column are either "Obligate" or facultatively cooperative
    } 
    
    # Birds of the World
    BOWCoop <- as.numeric(tempRowIn$Cooperative.BOW)
    
    # Add Mikula and Dale data 06/2023
    #Mikula <- tempRowIn$Coop_breed
    
    if (is.na(tempRowIn$Cooperative_breeding_ppca_Dale2015)) {
      Dale = NA
    } else if (tempRowIn$Cooperative_breeding_ppca_Dale2015 < 0) {
      Dale = 0
    } else if (tempRowIn$Cooperative_breeding_ppca_Dale2015 > 0) {
      Dale = 1
    }
    
    # add Cornwallis et al data 9/2023
    if (is.na(tempRowIn$Breeding.system_Cornwallis2017)) {
      Cornwallis = NA
    } else if (tempRowIn$Breeding.system_Cornwallis2017 == "Noncooperative") {
      Cornwallis = 0
    } else if (tempRowIn$Breeding.system_Cornwallis2017 == "Cooperative") {
      Cornwallis = 1
    }
    
    
    if (is.na(tempRowIn$Kin_Riehl2013)) {
      RiehlKin <- NA
    } else if (tempRowIn$Kin_Riehl2013 %in% c("K", "NK", "M")) {
      RiehlKin <- tempRowIn$Kin_Riehl2013
    } else if (tempRowIn$Kin_Riehl2013 %in% c("K; M", "K; NK")) {
      RiehlKin <- "M"
    } else {
      RiehlKin <- NA
    }
 
    # Griesser 2017 - family living vs non family, and coop vs noncoop
    if (is.na(tempRowIn$social_system_incl_nk_coop_Griesser2017)) {
      Griesser2017_Coop <- NA
      Griesser2017_KinCoop <- NA
    } else if (tempRowIn$social_system_incl_nk_coop_Griesser2017 %in% c("no_fam", "family")) {
      Griesser2017_Coop <- 0
      Griesser2017_KinCoop <- NA
    } else if (tempRowIn$social_system_incl_nk_coop_Griesser2017 == "nk-coop") {
      Griesser2017_Coop <- 1
      Griesser2017_KinCoop <- "NK"
    } else if (tempRowIn$social_system_incl_nk_coop_Griesser2017 == "coop_families") {
      Griesser2017_Coop <- 1
      Griesser2017_KinCoop <- "K"
    } else {
      Griesser2017_Coop <- NA
      Griesser2017_KinCoop <- NA
    }
    
    if (is.na(tempRowIn$social_system_Griesser2017)) {
      Griesser2017_Familial <- NA
    } else if (tempRowIn$social_system_Griesser2017 %in% c("coop_fam", "fam")) {
      Griesser2017_Familial <- 1
    } else if (tempRowIn$social_system_Griesser2017 == "no_fam") {
      Griesser2017_Familial <- 0 
    }
    
    FemaleSong <- tempRowIn$FemaleSong_Agg01
    
    # Dale et al added 6/20/2023
    allsourcesvec <- c(Dunn, Biagolini, DowningCoop, Jetz, Rubenstein, Cockburn, ReihlCoop, Griesser2017_Coop, BOWCoop, Dale, Cornwallis) # Mikula not included because they only added one new species compared to Dale et al (which has many more species than Mikula et al) and otherwise totally agreed with Dale. The one species in Mikula but not Dale is already present in the dataset 
    allsourcesNoNA <- na.omit(allsourcesvec)
    meanclass <- sum(as.numeric(allsourcesvec), na.rm = T)/length(allsourcesNoNA) 
    
    
    numSourcesNonCoop <- sum(allsourcesvec == 0, na.rm = TRUE)
    numSourcesCoop <- sum(allsourcesvec == 1, na.rm = TRUE)
    
    if (numSourcesNonCoop != 0 & numSourcesCoop != 0) {
      SourceDiscrepancy <- 1
    } else {
      SourceDiscrepancy <- 0
    }
    
    if (!is.na(meanclass)) {
      if (meanclass > 0.5) { # 3/8/2022 changed from >=
        MeanCoopOmitTies <- 1
      } else if (meanclass < 0.5) {
        MeanCoopOmitTies <- 0
      } else {MeanCoopOmitTies <- NA}
      
      if (1 %in% allsourcesvec) {
        AnyCoopEqualsCoop = 1
      } else if (0 %in% allsourcesvec) {
        AnyCoopEqualsCoop = 0
      } else {AnyCoopEqualsCoop = NA}
      
      if (0 %in% allsourcesvec) {
        AnyNoncoopEqualsNoncoop = 0
      } else if (1 %in% allsourcesvec) {
        AnyNoncoopEqualsNoncoop = 1
      } else {AnyNoncoopEqualsNoncoop = NA}
      
      if (meanclass > 0.5) { # 3/8/2022 changed from >=
        MeanCoopTie2Noncoop <- 1
      } else if (meanclass <= 0.5) {
        MeanCoopTie2Noncoop <- 0
      } else {MeanCoopTie2Noncoop <- NA}
      
      if (meanclass >= 0.5) { # 3/8/2022 changed from >=
        MeanCoopTie2Coop <- 1
      } else if (meanclass < 0.5) {
        MeanCoopTie2Coop <- 0
      } else {MeanCoopTie2Coop <- NA}
    } else {  # end if !is.na(meanclass)
      MeanCoopOmitTies <- NA
      AnyCoopEqualsCoop <- NA
      MeanCoopTie2Noncoop <- NA
      MeanCoopTie2Coop <- NA
      AnyNoncoopEqualsNoncoop = NA
    }  # end else
        
    allkinsources <- c(Griesser2017_KinCoop, RiehlKin, DowningKinCoop)
    numKin <- sum(allkinsources == "K", na.rm = TRUE)
    numNonKin <- sum(allkinsources == "NK", na.rm = TRUE)
    numMixedKinNonKin <- sum(allkinsources == "M", na.rm = TRUE)
    if (numKin + numNonKin + numMixedKinNonKin == 0) {
      Kin_NK <- NA
    } else if (numKin > numNonKin) {
      Kin_NK <- "Kin"
    } else if (numKin < numNonKin) {
      Kin_NK <- "NonKin"
    } else if (numMixedKinNonKin > numKin + numNonKin) {
      Kin_NK <- "Mixed"
    }
    
    
    # this comes after all the if...else statements for each column - putting the whole thing together!
    tempRowOut <- c(species, MeanCoopOmitTies, MeanCoopTie2Noncoop, MeanCoopTie2Coop, AnyCoopEqualsCoop, AnyNoncoopEqualsNoncoop, Kin_NK, numSourcesNonCoop, numSourcesCoop, SourceDiscrepancy, Dunn, Biagolini, DowningCoop, Jetz, JetzInclCockburn, Rubenstein, Cockburn, ReihlCoop, Griesser2017_Coop, BOWCoop, Dale, Cornwallis, RiehlKin, DowningKinCoop, Griesser2017_KinCoop, Griesser2017_Familial, numKin, numNonKin, numMixedKinNonKin)  # removed FemaleSong, added Mikula and Dale; on 2/24/2024, removed Mikula
    summarydf <- rbind(summarydf, tempRowOut)
    colnames(summarydf) <- c("species","MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Kin_NK", "numSourcesNonCoop", "numSourcesCoop", "SourceDiscrepancy", "DunnCoop","BiagoliniCoop", "DowningCoop", "JetzCoop", "JetzCoopInclCockburn", "RubensteinCoop", "CockburnCoop", "ReihlCoop", "Griesser2017Coop", "BOWCoop", "DaleCoop", "CornwallisCoop", "RiehlKin", "DowningKinNKCoop", "Griesser2017KinCoop", "Griesser2017FamilialLiving", "numKin", "numNonKin", "numMixed") # removed FemaleSong, added Mikula and Dale; 2/24/24 removed Mikula
    
  }  # end for (i in 1:length(ourdf$species_in_birdtree))
  
  summarydf <- as.data.frame(summarydf)
  
  summaryNoDups <- summarydf[which(!duplicated(summarydf$species)),] 
  
  summaryNoDups <- summaryNoDups[which(!summaryNoDups$species %in% c(NA, "(non-passerine)", "(not_checked)")),]

  write.csv(summaryNoDups, file = paste0(Sys.Date(), "_working_coop_breed.csv"), row.names = FALSE)
  
  CoopSourcesdf = merge(summaryNoDups, ourdf[,c("BirdtreeSpecies", "FemaleSong_Agg01", "HighConfidence_FemaleSong", "O.C", "System_Jetz2011", "Source_Jetz2011", "breeding.system_Downing2015", "justification_Downing2015", "cooperation.references_Downing2015", "KnownParentalCare_Cockburn2006", "InferredParentalCare_Cockburn2006", "Source_Cockburn2006", "MateSys5_Dunn2015", "BS_Biagolini2017", "O.F_Riehl2013", "Kin_Riehl2013", "Frequency_Riehl2013", "Social_breeding_system_when_cooperative_Riehl2013", "Social_System_RubensteinLovette2007", "engage.in.misdirected.parental.care_Griesser2016", "parental.care.mode_Griesser2016", "number.of.Zoological.Record.entries_Griesser2016","Griesser2017_speciesnames_Griesser2017", "social_system_Griesser2017", "social_system_incl_nk_coop_Griesser2017", "social_system_assessment_Griesser2017", "zoological_record_hits_Griesser2017", "Cooperative.BOW", "BOW.Quote", "Cooperative_breeding_ppca_Dale2015", "Cooperation_Remes2015", "Breeding.system_Cornwallis2017", "Reference.for.Breeding.System.1_Cornwallis2017", "Reference.for.Breeding.System.2_Cornwallis2017", "Reference.for.Breeding.System.3_Cornwallis2017", "Family3_BirdtreeMatchSpecies2_AVONET", "Order_AVONET", "EPP.10..threshold","EPP.source.s.",  "EPP.Source.data", "promiscuity...._Downing2015", "promiscuity.reference_Downing2015", "EPP1_Biagolini2017", "EPB1_Biagolini2017", "EPP2_Biagolini2017", "EPB2_Biagolini2017", "Ref._Biagolini2017", "EPP_Remes2015", "EPY_Remes2015", "Refs_Paternity_Remes2015", "Percentage.of.extra.group.paternity_Cornwallis2017", "Reference.for.Parentage_Cornwallis2017", colnames(ourdf)[which(str_detect(colnames(ourdf), "Griesser2023"))] )], by.x = "species", by.y = "BirdtreeSpecies", all = T)
  print(CoopSourcesdf$species[which(duplicated(CoopSourcesdf$species))])
  write.csv(CoopSourcesdf, file = paste0(Sys.Date(),"_CoopClassesWSourceColumns.csv"), row.names = FALSE)
  
  # subset to Passeriformes and binarize sociality data
  ordercounts = CoopSourcesdf %>% group_by(Order_AVONET) %>% count
  allPass = CoopSourcesdf[which(CoopSourcesdf$Order_AVONET == "Passeriformes"),]
  
  unique(allPass$colonial_Griesser2023)
  allPass$Griesser2023.Colonial01 = NA
  allPass$Griesser2023.Colonial01[which(allPass$colonial_Griesser2023 == "acolonial")] <- 0
  allPass$Griesser2023.Colonial01[which(allPass$colonial_Griesser2023 == "colonial")] <- 1
  
  unique(allPass$caretakers_Griesser2023)
  allPass$Griesser2023.MoreThanTwoCaretakers = NA
  allPass$Griesser2023.MoreThanTwoCaretakers[which(allPass$caretakers_Griesser2023 <= 2 )] <- 0
  allPass$Griesser2023.MoreThanTwoCaretakers[which(allPass$caretakers_Griesser2023 > 2 )] <- 1
  allPass$Griesser2023.TwoOrMoreCaretakers = NA
  allPass$Griesser2023.TwoOrMoreCaretakers[which(allPass$caretakers_Griesser2023 < 2 )] <- 0
  allPass$Griesser2023.TwoOrMoreCaretakers[which(allPass$caretakers_Griesser2023 >= 2 )] <- 1
  
  unique(allPass$social_bonds_Griesser2023)
  allPass %>% group_by(social_bonds_Griesser2023) %>% count
  allPass$Griesser2023.LongSocialBonds = NA
  allPass$Griesser2023.LongSocialBonds[which(allPass$social_bonds_Griesser2023 %in% c("a-short", "b-season"))] <- 0
  allPass$Griesser2023.LongSocialBonds[which(allPass$social_bonds_Griesser2023 == "c-long")] <- 1
  allPass$Griesser2023.SeasonOrLongerSocialBonds = NA
  allPass$Griesser2023.SeasonOrLongerSocialBonds[which(allPass$social_bonds_Griesser2023 %in% c("a-short"))] <- 0
  allPass$Griesser2023.SeasonOrLongerSocialBonds[which(allPass$social_bonds_Griesser2023 %in% c("c-long", "b-season"))] <- 1
  
  unique(allPass$grouping_Griesser2023)
  allPass %>% group_by(grouping_Griesser2023) %>% count
  allPass$Griesser2023.GroupsLargerThanPair = NA
  allPass$Griesser2023.GroupsLargerThanPair[which(allPass$grouping_Griesser2023 %in% c("asocial", "pair"))] <- 0
  allPass$Griesser2023.GroupsLargerThanPair[which(allPass$grouping_Griesser2023 %in% c("small_groups", "large_groups"))] <- 1
  allPass$Griesser2023.Asocial0VsSocial1 = NA
  allPass$Griesser2023.Asocial0VsSocial1[which(allPass$grouping_Griesser2023 %in% c("asocial"))] <- 0
  allPass$Griesser2023.Asocial0VsSocial1[which(allPass$grouping_Griesser2023 %in% c("small_groups", "large_groups", "pair"))] <- 1
  allPass$Griesser2023.LargestGroupSizes = NA
  allPass$Griesser2023.LargestGroupSizes[which(allPass$grouping_Griesser2023 %in% c("small_groups", "asocial", "pair"))] <- 0
  allPass$Griesser2023.LargestGroupSizes[which(allPass$grouping_Griesser2023 %in% c( "large_groups"))] <- 1
  
  
  SnyderCreanzaData <- read.csv(file = songfile)
  SnyderCreanzaData$Family = NULL
  coopsongdf <- merge(allPass, SnyderCreanzaData, by.x = "species", by.y = "BirdtreeFormat", all = TRUE)
  
  write.csv(coopsongdf, paste0(Sys.Date(),"_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"), row.names = FALSE)
  
  # adapted from end of merge_data_allcolumns.R
  if (is.character(HighConfCoopFile)) {
    fulldf = read.csv(HighConfCoopFile)
    CoopSongSocPasser = merge(fulldf[,c("species", "HighConfidence_Coop")], coopsongdf, by = "species", all.y = T)
    write.csv(CoopSongSocPasser, paste0(Sys.Date(), "_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"), row.names = F)
  }

} # end function
