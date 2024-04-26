## Generate summary df for species cooperative breeding data
## Coded by Aleyna Loughran-Pierce and Kate Snyder
## Created 10/9/2020
## Last edited: 8/18/2021 by Kate Snyder
## 3/8/2022 - add nonkin/kin etc, added Odom Female Song to output, added BOW data, changed MeanCoop to be only > 0.5 (not >=)
## 3/9/2022 - output all Classification methods as different columns in the same table, remove coopClassMethod arg, remove extra rows
## 3/10/2022 - correct Jetz et al classification to reflect more liberal use of Cockburn et al data; add OC data
## 5/30/2023 - commented out last write csv since song data already added; removed Female Song from output since not included in input yet
## 6/14/2023 - added Mikula and Dale columns to output
## 6/20/2023 - included Dale in MeanCoop calculation
## 9/13/2023 - added Cornwallis et al data handling
## 10/25/2023 - just made subsetlabel NULL for non-subsetted data
## 2/23/2024 - archived this code copy

setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")

#CoopBreed_species_summary(coopbreedfile = "Aggregate_Source_Data_AllColumns.csv", allcoop = TRUE, OCdatafile = TRUE)
#CoopBreed_species_summary(coopbreedfile = "2023-05-30_Aggregate_Source_Data_AllColumns.csv", allcoop = TRUE, OCdatafile = TRUE)
#CoopBreed_species_summary(coopbreedfile = coopbreedfile, allcoop = TRUE, OCdatafile = TRUE)

#coopbreedfile = "2023-05-30_Aggregate_Source_Data_AllColumns.csv"
#coopbreedfile = "/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_Aggregate_Source_Data_AllColumns.csv"

CoopBreed_species_summary <- function(coopbreedfile, songfile = "SongData_R_Update.csv", allcoop = FALSE, OCdatafile = FALSE) {
  
  ourdf <- read.csv(file = coopbreedfile)
  
  summarydf <- set.seed(10)  #this creates a blank dataframe for the for loop to fill in
  
  for (i in 1:length(ourdf$BirdtreeSpecies)) {
    tempRowOut <- set.seed(10) # this creates a blank new row for our summary dataframe
    tempRowIn <- ourdf[i,] 
    species <- tempRowIn$BirdtreeSpecies
    
    if (is.na(tempRowIn$MateSys5)) { 
      Dunn <- NA
    } else if (tempRowIn$MateSys5 == "Coop") {  
      Dunn <- 1
    } else {
      Dunn <- NA  
    }
    
    ### This one I did a little differently because there were so many different "values" in this column, and those that contain the character string "CB" will get a 1, and any that don't will get a 0. The function str_detect() is in the package 'stringr'. Feel free to use this template if this kind of thing occurs in another column!
    require(stringr)  # checks whether the package is loaded, and loads it if not - i.e. the equivalent of using library(stringr); it does need to be installed first  -  in the bottom right quadrant of Rstudio, go to the "Packages" tab --> "Install" --> enter "stringr" and click "Install"
    if (is.na(tempRowIn$BS)) {
      Biagolini <- NA
    } else if (str_detect(tempRowIn$BS, "CB", negate = TRUE)) {
      Biagolini <- 0
    } else if (str_detect(tempRowIn$BS, "CB", negate = FALSE)) {
      Biagolini <- 1
    } else {
      Biagolini <- NA
    }
    
    
    if (is.na(tempRowIn$breeding.system)) { 
      DowningCoop <- NA
      DowningKinCoop <- NA
    } else if (tempRowIn$breeding.system == "noncooperative") {
      DowningCoop <- 0
      DowningKinCoop <- NA
    } else if (tempRowIn$breeding.system %in% c("kinCooperator")) {  
      DowningCoop <- 1
      DowningKinCoop <- "K"
    } else if (tempRowIn$breeding.system %in% c("nonkinCooperator")) {  
      DowningCoop <- 1
      DowningKinCoop <- "NK"
    } else {
      DowningCoop <- NA  
      DowningKinCoop <- NA
    }
    
    
    if (is.na(tempRowIn$Social_System)) {
      Rubenstien <- NA
    } else if (tempRowIn$Social_System == "Noncooperative") {
      Rubenstien <- 0
    } else if (tempRowIn$Social_System %in% c("Cooperative", " Cooperative")) {
      Rubenstien <- 1
    } else {
      Rubenstien <- NA
    }
    
    
    if (is.na(tempRowIn$KnownParentalCare)) { 
      Cockburn<- NA
    } else if (tempRowIn$KnownParentalCare %in% c("Pair", "Female_only", "Suspected", "Male_only")) {
      Cockburn <- 0
    } else if (tempRowIn$KnownParentalCare == "Cooperation") {  
      Cockburn <- 1
    } else {
      Cockburn<- NA  
    }
    
    if (is.na(tempRowIn$System)) {
      Jetz <- NA
    } else if (tempRowIn$System == "Non-cooperative") {
      Jetz <- 0
    } else if (tempRowIn$System == "Cooperative") {
      Jetz <- 1
    } else {
      Jetz <- NA
    }
    
    JetzSource <- as.character(tempRowIn$Source)
    JetzSources <- as.character(-17:-2)
    
    if (!is.na(Jetz) & !is.na(Cockburn)) {  
      if  (JetzSource == "-1") {    # if both Jetz and Cockburn have a CoopBreed classification, and Jetz's source is Cockburn
        Jetz <- NA
      } else if (JetzSource %in% JetzSources) {
        Jetz <- Jetz
      }
    }
    
    if (is.na(tempRowIn$Social_breeding_system_when_cooperative)) {
      ReihlCoop <- NA
    } else {
      ReihlCoop <- 1  # because all entries in this column are cooperative
    } 
    
    # Birds of the World
    BOWCoop <- tempRowIn$Cooperative.BOW
    
    # Add Mikula and Dale data 06/2023
    Mikula <- tempRowIn$Coop_breed
    
    if (is.na(tempRowIn$Cooperative_breeding_ppca)) {
      Dale = NA
    } else if (tempRowIn$Cooperative_breeding_ppca < 0) {
      Dale = 0
    } else if (tempRowIn$Cooperative_breeding_ppca > 0) {
      Dale = 1
    }
    
    # add Cornwallis et al data 9/2023
    if (is.na(tempRowIn$Breeding.system)) {
      Cornwallis = NA
    } else if (tempRowIn$Breeding.system == "Noncooperative") {
      Cornwallis = 0
    } else if (tempRowIn$Breeding.system == "Cooperative") {
      Cornwallis = 1
    }
    
    
    if (is.na(tempRowIn$Kin)) {
      RiehlKin <- NA
    } else if (tempRowIn$Kin %in% c("K", "NK", "M")) {
      RiehlKin <- tempRowIn$Kin
    } else if (tempRowIn$Kin %in% c("K; M", "K; NK")) {
      RiehlKin <- "M"
    } else {
      RiehlKin <- NA
    }
 
    # Griesser 2017 - family living vs non family, and coop vs noncoop
    if (is.na(tempRowIn$social_system_incl_nk_coop)) {
      Griesser2017_Coop <- NA
      Griesser2017_KinCoop <- NA
    } else if (tempRowIn$social_system_incl_nk_coop %in% c("no_fam", "family")) {
      Griesser2017_Coop <- 0
      Griesser2017_KinCoop <- NA
    } else if (tempRowIn$social_system_incl_nk_coop == "nk-coop") {
      Griesser2017_Coop <- 1
      Griesser2017_KinCoop <- "NK"
    } else if (tempRowIn$social_system_incl_nk_coop == "coop_families") {
      Griesser2017_Coop <- 1
      Griesser2017_KinCoop <- "K"
    } else {
      Griesser2017_Coop <- NA
      Griesser2017_KinCoop <- NA
    }
    
    if (is.na(tempRowIn$social_system)) {
      Griesser2017_Familial <- NA
    } else if (tempRowIn$social_system %in% c("coop_fam", "fam")) {
      Griesser2017_Familial <- 1
    } else if (tempRowIn$social_system == "no_fam") {
      Griesser2017_Familial <- 0 
    }
    
    FemaleSong <- tempRowIn$Female.Song.Score
    
    # Dale et al added 6/20/2023
    allsourcesvec <- c(Dunn, Biagolini, DowningCoop, Jetz, Rubenstien, Cockburn, ReihlCoop, Griesser2017_Coop, BOWCoop, Dale, Cornwallis) # Mikula not included because they only added one new species compared to Dale et al (which has many more species than Mikula et al) and otherwise totally agreed with Dale. The one species in Mikula but not Dale is already present in the dataset 
    allsourcesNoNA <- na.omit(allsourcesvec)
    meanclass <- sum(allsourcesNoNA)/length(allsourcesNoNA) 
    
    
    numSourcesNonCoop <- sum(allsourcesvec == 0, na.rm = TRUE)
    numSourcesCoop <- sum(allsourcesvec == 1, na.rm = TRUE)
    
    if (numSourcesNonCoop != 0 & numSourcesCoop != 0) {
      SourceDiscrepancy <- 1
    } else {
      SourceDiscrepancy <- 0
    }
    
    #if (coopClassMethod == "MeanCoopOmitTies") {
      if (!is.na(meanclass)) {
        if (meanclass > 0.5) { # 3/8/2022 changed from >=
          MeanCoopOmitTies <- 1
        } else if (meanclass < 0.5) {
          MeanCoopOmitTies <- 0
        } else {MeanCoopOmitTies <- NA}
    #  } else {MeanCoopOmitTies <- NA}
      
    #} else if (coopClassMethod == "AnyCoopEqualsCoop") {
      if (1 %in% allsourcesvec) {
        AnyCoopEqualsCoop = 1
      } else if (0 %in% allsourcesvec) {
        AnyCoopEqualsCoop = 0
      } else {AnyCoopEqualsCoop = NA}
    #}
    
    #if (coopClassMethod == "MeanCoopTie2NonCoop") {
     # if (!is.na(meanclass)) {
        if (meanclass > 0.5) { # 3/8/2022 changed from >=
          MeanCoopTie2Noncoop <- 1
        } else if (meanclass <= 0.5) {
          MeanCoopTie2Noncoop <- 0
        } else {MeanCoopTie2Noncoop <- NA}
      #} else {MeanCoopTie2Noncoop <- NA}
    
      #if (!is.na(meanclass)) {
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
  #  tempRowOut <- c(species, MeanCoopOmitTies, MeanCoopTie2Noncoop, MeanCoopTie2Coop, AnyCoopEqualsCoop, Kin_NK, numSourcesNonCoop, numSourcesCoop, SourceDiscrepancy, Dunn, Biagolini, DowningCoop, Jetz, Rubenstien, Cockburn, ReihlCoop, Griesser2017_Coop, BOWCoop, RiehlKin, DowningKinCoop, Griesser2017_KinCoop, Griesser2017_Familial, numKin, numNonKin, numMixedKinNonKin, FemaleSong, JetzSource)  
    tempRowOut <- c(species, MeanCoopOmitTies, MeanCoopTie2Noncoop, MeanCoopTie2Coop, AnyCoopEqualsCoop, Kin_NK, numSourcesNonCoop, numSourcesCoop, SourceDiscrepancy, Dunn, Biagolini, DowningCoop, Jetz, Rubenstien, Cockburn, ReihlCoop, Griesser2017_Coop, BOWCoop, Mikula, Dale, Cornwallis, RiehlKin, DowningKinCoop, Griesser2017_KinCoop, Griesser2017_Familial, numKin, numNonKin, numMixedKinNonKin, JetzSource)  # removed FemaleSong, added Mikula and Dale
    summarydf <- rbind(summarydf, tempRowOut)
   # colnames(summarydf) <- c("species","MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop", "Kin_NK", "numSourcesNonCoop", "numSourcesCoop", "SourceDiscrepancy", "Dunn","Biagolini", "DowningCoop", "Jetz", "Rubenstein", "Cockburn", "ReihlCoop", "Griesser2017Coop", "BOWCoop", "RiehlKin", "DowningKinNKCoop", "Griesser2017KinCoop", "Griesser2017FamilialLiving", "numKin", "numNonKin", "numMixed", "FemaleSong", "JetzSource") # must have the same length as tempRowOut
    colnames(summarydf) <- c("species","MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop", "Kin_NK", "numSourcesNonCoop", "numSourcesCoop", "SourceDiscrepancy", "Dunn","Biagolini", "DowningCoop", "Jetz", "Rubenstein", "Cockburn", "ReihlCoop", "Griesser2017Coop", "BOWCoop", "MikulaCoop", "DaleCoop", "CornwallisCoop", "RiehlKin", "DowningKinNKCoop", "Griesser2017KinCoop", "Griesser2017FamilialLiving", "numKin", "numNonKin", "numMixed", "JetzSource") # removed FemaleSong, added Mikula and Dale
    
  }  # end for (i in 1:length(ourdf$species_in_birdtree))
  
  summarydf <- as.data.frame(summarydf)
  #summarydfexpanded <- cbind(ourdf[,"species_in_birdtree"], #,"OurDatabaseOrder","OurDatabaseFamily","OurDatabaseEPP")], 
  #                           summarydf[,2:length(colnames(summarydf))])
  
  summaryNoDups <- summarydf[which(!duplicated(summarydf$species)),] 
  
  summaryNoDups <- summaryNoDups[which(!summaryNoDups$species %in% c(NA, "(non-passerine)", "(not_checked)")),]
  

  write.csv(summaryNoDups, file = paste0(Sys.Date(), "_working_coop_breed.csv"), row.names = FALSE)
  
  
  SnyderCreanzaData <- read.csv(file = songfile)
  coopsongdf <- merge(summaryNoDups, SnyderCreanzaData, by.x = "species", by.y = "BirdtreeFormat", all.x = allcoop, all.y = TRUE)
  
  if (OCdatafile != FALSE) {
  OCdata <- read.csv(file = "OCPaperData.csv")
  coopsongdf <- merge(coopsongdf, OCdata, by.x = "species", by.y = "BirdtreeFormat", all.x = allcoop, all.y = TRUE)
  }
  
  if (allcoop == FALSE) {
    subsetlabel <- "_NatCommsSubset"
  } else if (allcoop == TRUE) {
    subsetlabel <- ""   # "_All" # changed from _All 10/26/2023
  }
  
  
  write.csv(coopsongdf, file = paste0(Sys.Date(),"CoopSong", subsetlabel, "_R.csv"), row.names = FALSE)
  
} # end function
