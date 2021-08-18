## Merge individual source datasets
## Created by Aleyna Loughran-Pierce
## Created 10/10/2020
## Last modified: 8/17/2021 by Kate Snyder

ourdatabase <- read.csv("SongData_R_Update.csv", stringsAsFactors = FALSE)
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
DunnData <- read.csv("Dunn supp data 977sp.csv", stringsAsFactors = FALSE)
GreisserData <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE)
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
RiehlData <- read.csv("Riehl 2013 supp data columns lines_kts edited.csv", stringsAsFactors = FALSE)
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)

# merge Jetz + Downing
DowningSubset <- DowningData[,c(1,2,4,11,12,13,14)]
JetzSubset <- JetzData[, c(1,2,3,4,7,8)]
newdf <- merge(x = JetzSubset, y = DowningSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_Jetz", "_Downing"))
colnames(newdf)[2] <- "Species_Jetz"
colnames(newdf)[7] <- "Promiscuity_Downing"

# merge newdf and Cockburn
CockburnSubset <- CockburnData[, c(1,2,5,6,7)]
newdf2 <- merge(newdf, CockburnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Cockburn"))
colnames(newdf2)[which(colnames(newdf2) == "SpeciesScientific")] <- "SpeciesScientific_Cockburn"

# merge newdf2 and Dunn
DunnSubset <- DunnData[, c(1,2,10,11,12,13,14,15,16,17,18)]
newdf3 <- merge(newdf2, DunnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Dunn"))
colnames(newdf3)[which(colnames(newdf3) == "CommonNameBL")] <- "CommonName_Dunn"

# merge newdf3 and Biagolini
BiagoliniSubset <- BiagoliniData[, c(1,2,4,5,6,7,9)]
newdf4 <- merge(newdf3, BiagoliniSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Biagolini"))
colnames(newdf4)[which(colnames(newdf4) == "Species")] <- "Species_Biagolini"


# merge newdf4 and Riehl
RiehlSubset <- RiehlData[, c(1,2,3,5,6,7,8)]
newdf5 <- merge(newdf4, RiehlSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Riehl"))
colnames(newdf5)[which(colnames(newdf5) == "Dispersal")] <- "Dispersal_Riehl"


# merge newdf5 and Rubenstein
RubensteinSubset <- RubensteinLovetteData[, c(1,2,3,4)]
newdf6 <- merge(newdf5, RubensteinSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Rubenstein"))
colnames(newdf6)[which(colnames(newdf6) == "Species")] <- "Species_Rubenstein"


# merge newdf6 and Griesser
GriesserSubset <- GreisserData[, c(1,2,3,4,5,6)]
newdf7 <- merge(newdf6, GriesserSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Griesser"))
colnames(newdf7)[which(colnames(newdf7) == "common.name")] <- "common.name_Griesser"

# write file
write.csv(newdf7, file = "Aggregate_Source_Data.csv", row.names = FALSE)
