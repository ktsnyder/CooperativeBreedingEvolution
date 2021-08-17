## Merge individual source datasets
## Created by Aleyna Loughran-Pierce
## Created 10/10/2020
## Last modified: 8/17/2021 by Kate Snyder

ourdatabase <- read.csv("20200602_Cooperative Breeding Bird Database.csv", stringsAsFactors = FALSE)
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
GreisserData <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE)
RiehlData <- read.csv("Riehl_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)

DowningData[c(1,2,3,4,5,10),]
removerows <- which(DowningData$species == "")
#removedrows <- DowningData[removerows,]
keeprows <- which(DowningData$species != "")
DowningSubset <- DowningData[keeprows,]

keeprows <- which(ourdatabase$Cooperative == "1")
keeprows <- which(ourdatabase$Cooperative %in% c("0", "1"))
OurdatabaseSubset <- ourdatabase[keeprows,]



newdf <- merge(OurdatabaseSubset, DowningData, by.y = "species", by.x = "SpeciesScientific", all = TRUE)
newdf <- merge(ourdatabase, DowningSubset, by.y = "species", by.x = "SpeciesScientific", all = TRUE)
newdf <- newdf[, 1:34]

SpeciesName <- set.seed(10) 
for (i in 1:length(JetzData$Index)) {
  SpeciesName[i] <- paste0(JetzData$Genus[i], "_", JetzData$Species[i])
}
SpeciesJetzData <- cbind(SpeciesName, JetzData)
SpeciesJetzData<- SpeciesJetzData[, 1:9]

newdf2 <- merge(newdf, SpeciesJetzData, by.x = "SpeciesScientific", by.y = "SpeciesName", all = TRUE)

SpeciesName <- set.seed(10)
for (i in 1:length(RubensteinLovetteData$Genus)) {
  SpeciesName[i] <- paste(RubensteinLovetteData$Genus[i], RubensteinLovetteData$Species[i], sep = '_')
}
SpeciesNameRubenstein <- str_remove(SpeciesName, " ")
SpeciesRubensteinLovetteData <- cbind(SpeciesNameRubenstein, RubensteinLovetteData)
SpeciesRubensteinLovetteData <- SpeciesRubensteinLovetteData[, 1:8]

newdf3 <- merge(newdf2, SpeciesRubensteinLovetteData, by.x = "SpeciesScientific", by.y = "SpeciesNameRubenstein", all = TRUE)


GreisserData$scientific.name[!(GreisserData$scientific.name %in% newdf3$SpeciesScientific)]

newdf4 <- merge(newdf3, GreisserData, by.y = "scientific.name", by.x  = "SpeciesScientific", all = TRUE)

CockburnDataAdditionalNames <- read.csv("cockburndatanewdf4missingnamesonbirdtree.csv", stringsAsFactors = FALSE)
newdf5 <- merge(CockburnDataAdditionalNames, RiehlData, by.x = "SpeciesScientific", by.y = "SpeciesName", all = TRUE)

newdf6 <- merge(newdf4, newdf5, by.x = "SpeciesScientific", by.y = "SpeciesScientific", all = TRUE)

write.csv(newdf6, file = "FinaldfCooperativeBreeding2020correct.csv")

#read in the two .csv files

#Replace misspelled birds in database
#
#

#setwd("~/Documents/Creanza Lab/Comparative Evolution/Cooperative Breeding")  #this only applies to my (Kate's) computer

misspelledbirds <- read.csv("misspelled_birds_to_check.csv")  # The version of this I sent you may not be the most updated version. if you have the most recent/correct version of this file, you can change this to whatever the name of that file is - a couple notes: I noticed that Phoenicurus_erythronota was mistakenly changed (I think by the code that I used to generate this file) to Phoenicurus_erythrogastrus, when it should be Phoenicurus_erythronotus. I did edit that in this file. Also, Apteryx_rowi, Ara_ambigua, Ara_manilata, and Ara_nobilis were I think also incorrectly corrected by my code, and I changed all of these to NA. You may very well have caught these already, and I def didn't go through the whole list, just wanted to make a note of what I actually manually changed.
ourdb <- read.csv("2020-09-11_GoogleDrive_FinaldfCooperativeBreeding2020correct.csv")

correctednamesdf <- merge(misspelledbirds,ourdb, by.x = "in_database", by.y = "SpeciesScientific", all = TRUE)
correctednamesdf2 <- correctednamesdf[,which(colnames(correctednamesdf) != "X.x" & colnames(correctednamesdf) != "X.y")] 
NAvector <- which(is.na(correctednamesdf2$species_in_birdtree))  # the rows where birds were spelled correctly (and now are NA in the corrected names column)
correctednamesdf2$species_in_birdtree[NAvector] <- correctednamesdf2$in_database[NAvector] # putting the correctly spelled names from in_database into the "NA" spots in species_in_birdtree

correctednamesdf3 <- correctednamesdf2[, 2:length(correctednamesdf2)]  # the next step doesn't work if the first column (w/ misspelled names) is there, this line gets rid of that column

# combine any duplicate species rows
require(dplyr)
correctednamesdf4 <- correctednamesdf3 %>% group_by(species_in_birdtree) %>% summarize_all(.funs=function(x) {if(length(na.omit(x))==0){NA} else {na.omit(x)}})

write.csv(correctednamesdf4, file = "v2020-09-11_GDrive_CoopBreeding_correctednames.csv")  # if you run this line, just make sure this .csv isn't open

SnyderCreanza_song <- read.csv("SnyderCreanza_NatComms2019_SupplementalData_R.csv", stringsAsFactors = FALSE)
CoopBreed_song <- merge(summaryNoDups, SnyderCreanza_song, by.x = "species_in_birdtree", by.y = "BirdtreeFormat", all = TRUE)
write.csv(CoopBreed_song, file = "Coop_Breed_song1.csv")
