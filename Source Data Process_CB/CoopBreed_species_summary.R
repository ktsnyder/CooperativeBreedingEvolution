## Generate summary df for species cooperative breeding data
## Coded by Aleyna Loughran-Pierce and Kate Snyder
## Created 10/9/2020
## Last edited: 8/18/2021 by Kate Snyder

setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")

CoopBreed_species_summary(coopbreedfile = "Aggregate_Source_Data.csv", coopClassMethod = "AnyCoopEqualsCoop", allcoop = FALSE)

CoopBreed_species_summary <- function(coopbreedfile, songfile = "SongData_R_Update.csv", coopClassMethod = c("MeanCoop", "AnyCoopEqualsCoop"), allcoop = FALSE) {
  
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
      Downing <- NA
    } else if (tempRowIn$breeding.system == "noncooperative") {
      Downing <- 0
    } else if (tempRowIn$breeding.system %in% c("kinCooperator","nonkinCooperator")) {  
      Downing <- 1
    } else {
      Downing <- NA  
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
    
    if (is.na(tempRowIn$Social_breeding_system_when_cooperative)) {
      Reihl <- NA
    } else {
      Reihl <- 1  # because all entries in this column are cooperative
    } 
    
    # ## Extract mating system from Biagolini
    # if (is.na(tempRowIn$OurDatabaseMateSys.Biagolini)) {
    #   BiagoliniPolygyny <- NA
    # } else if (tempRowIn$OurDatabaseMateSys.Biagolini %in% c("Social Monogamy", "Social Monogamy, Cooperative Breeder", "Social Monogamy,Cooperative Breeder")) {
    #   BiagoliniPolygyny <- 0
    # } else if (tempRowIn$OurDatabaseMateSys.Biagolini == "Polygyny") {
    #   BiagoliniPolygyny <- 1
    # } else {
    #   BiagoliniPolygyny <- NA
    # }
    
    allsourcesvec <- c(Dunn, Biagolini, Downing, Jetz, Rubenstien, Cockburn, Reihl)
    allsourcesNoNA <- na.omit(allsourcesvec)
    meanclass <- sum(allsourcesNoNA)/length(allsourcesNoNA) 
    
    numSourcesNonCoop <- sum(allsourcesvec == 0, na.rm = TRUE)
    numSourcesCoop <- sum(allsourcesvec == 1, na.rm = TRUE)
    
    if (numSourcesNonCoop != 0 & numSourcesCoop != 0) {
      SourceDiscrepancy <- 1
    } else {
      SourceDiscrepancy <- 0
    }
    
    if (coopClassMethod == "MeanCoop") {
      if (!is.na(meanclass)) {
        if (meanclass >= 0.5) {
          CoopBreed <- 1
        } else if (meanclass < 0.5) {
          CoopBreed <- 0
        } else {CoopBreed <- NA}
      } else {CoopBreed <- NA}
      
    } else if (coopClassMethod == "AnyCoopEqualsCoop") {
      if (1 %in% allsourcesvec) {
        CoopBreed = 1
      } else if (0 %in% allsourcesvec) {
        CoopBreed = 0
      } else {CoopBreed = NA}
    }
    
    # this comes after all the if...else statements for each column - putting the whole thing together!
    tempRowOut <- c(species, CoopBreed, numSourcesNonCoop, numSourcesCoop, SourceDiscrepancy, Dunn, Biagolini, Downing, Jetz, Rubenstien, Cockburn, Reihl)  
    summarydf <- rbind(summarydf, tempRowOut)
    colnames(summarydf) <- c("species","CoopBreed", "numSourcesNonCoop", "numSourcesCoop", "SourceDiscrepancy", "Dunn","Biagolini", "Downing", "Jetz", "Rubenstein", "Cockburn", "Reihl") # must have the same length as tempRowOut
    
  }  # end for (i in 1:length(ourdf$species_in_birdtree))
  
  summarydf <- as.data.frame(summarydf)
  #summarydfexpanded <- cbind(ourdf[,"species_in_birdtree"], #,"OurDatabaseOrder","OurDatabaseFamily","OurDatabaseEPP")], 
  #                           summarydf[,2:length(colnames(summarydf))])
  
  summaryNoDups <- summarydf[which(!duplicated(summarydf$species)),] #not necessary anymore
  
  write.csv(summaryNoDups, file = paste0(Sys.Date(), "_working_coop_breed.csv"), row.names = FALSE)
  
  
  SnyderCreanzaData <- read.csv(file = songfile)
  coopsongdf <- merge(summaryNoDups, SnyderCreanzaData, by.x = "species", by.y = "BirdtreeFormat", all.x = allcoop, all.y = TRUE)
  
  if (allcoop == FALSE) {
    subsetlabel <- "_NatCommsSubset"
  } else if (allcoop == TRUE) {
    subsetlabel <- "_All"
  }
  
  write.csv(coopsongdf, file = paste0(Sys.Date(),"CoopSong_", coopClassMethod, subsetlabel, ".csv"), row.names = FALSE)
  
} # end function
