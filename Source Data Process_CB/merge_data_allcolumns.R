## Merge individual source datasets
## Created by Aleyna Loughran-Pierce
## Created 10/10/2020
## Last modified: 12/09/2021 by Kate Snyder - include all columns from each source file
## 3/7/2022 - add in Griesser 2017 data, female song (Odom 2014) data

ourdatabase <- read.csv("SongData_R_Update.csv", stringsAsFactors = FALSE)
ourdatabaserefs <- read.csv("SupplementDataRefs_Update.csv")
# colnames(ourdatabaserefs)[1] <- "BirdtreeSpecies"  # changed in csv
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
DunnData <- read.csv("Dunn supp data 977sp.csv", stringsAsFactors = FALSE)
GriesserData <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE)
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
RiehlData <- read.csv("Riehl 2013 supp data columns lines_kts edited.csv", stringsAsFactors = FALSE)
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
OdomFSData <- read.csv()
Griesser2017Data <- read.csv()

# merge Jetz + Downing
DowningSubset <- DowningData #[,c(1,2,4,11,12,13,14)]
JetzSubset <- JetzData #[, c(1,2,3,4,7,8)]
newdf <- merge(x = JetzSubset, y = DowningSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_Jetz", "_Downing"))
colnames(newdf)[which(colnames(newdf) == "Species")] <- "Species_Jetz"
colnames(newdf)[which(colnames(newdf) == "promiscuity....")] <- "Promiscuity_Downing"

# merge newdf + ourdatabase
newdf_wOurs <- merge(x = ourdatabaserefs, y = newdf, by = "BirdtreeSpecies", all = TRUE)


# merge newdf_wOurs and Cockburn
CockburnSubset <- CockburnData #[, c(1,2,5,6,7)]
colnames(CockburnSubset)[which(colnames(CockburnSubset) == "Species")] <- "Sp"
newdf2 <- merge(newdf_wOurs, CockburnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Cockburn"))
colnames(newdf2)[which(colnames(newdf2) == "SpeciesScientific")] <- "SpeciesScientific_Cockburn"

# merge newdf2 and Dunn
DunnSubset <- DunnData #[, c(1,2,10,11,12,13,14,15,16,17,18)]
newdf3 <- merge(newdf2, DunnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Dunn"))
colnames(newdf3)[which(colnames(newdf3) == "CommonNameBL")] <- "CommonName_Dunn"

# merge newdf3 and Biagolini
BiagoliniSubset <- BiagoliniData #[, c(1,2,4,5,6,7,9)]
newdf4 <- merge(newdf3, BiagoliniSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Biagolini"))
colnames(newdf4)[which(colnames(newdf4) == "Species")] <- "Species_Biagolini"


# merge newdf4 and Riehl
RiehlSubset <- RiehlData #[, c(1,2,3,5,6,7,8)]
newdf5 <- merge(newdf4, RiehlSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Riehl"))
colnames(newdf5)[which(colnames(newdf5) == "Dispersal")] <- "Dispersal_Riehl"


# merge newdf5 and Rubenstein
RubensteinSubset <- RubensteinLovetteData #[, c(1,2,3,4)]
newdf6 <- merge(newdf5, RubensteinSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Rubenstein"))
colnames(newdf6)[which(colnames(newdf6) == "Species")] <- "Species_Rubenstein"


# merge newdf6 and Griesser
GriesserSubset <- GriesserData #[, c(1,2,3,4,5,6)]
newdf7 <- merge(newdf6, GriesserSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Griesser2016"))
colnames(newdf7)[which(colnames(newdf7) == "common.name")] <- "common.name_Griesser2016"

# write file
write.csv(newdf7, file = "Aggregate_Source_Data_AllColumns.csv", row.names = FALSE)


# Additional (non-Cooperative Breeding) Sources
LislevandData <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Lislevand et al 2007 supp data.csv")
SpottiswoodeData <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Spottiswoode and Moller 2004_data_noWeirdCharacters.csv")
DatabaseUpdate2019 <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Song Database Update 2019_LastEdited_20200408_539PM_LislevandSpottiswoodeClementsOnly.csv")
Clements2021 <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/eBird-Clements-v2021-integrated-checklist-August-2021_kts_speciesonly.csv")
Clements2021subset <- Clements2021[,c(6,5,7,8,9,10,11)]
colnames(DatabaseUpdate2019)[which(colnames(DatabaseUpdate2019) == "BirdtreeFormat")] <- "BirdtreeSpecies"
DatabaseUpdate2019_Clements <- DatabaseUpdate2019[,c(1:7,9,10,12,13,14,15,16,17)]
DatabaseUpdate2019_LisSpo <- DatabaseUpdate2019[,c(1,18:64)]

#merge DatabaseUpdate2019_Clements and newdf7
newdf8 <- merge(DatabaseUpdate2019_Clements, newdf7, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_ClementsGuide", ""))
newdf9 <- merge(newdf8, DatabaseUpdate2019_LisSpo, by = "BirdtreeSpecies", all = TRUE)
newdf10 <- merge(Clements2021subset, newdf9, by.x = "Scientific", by.y = "BirdtreeSpecies", all.x = FALSE, all.y = TRUE, suffixes = c("ClementsGuide2021",))
write.csv(newdf10, file = "Aggregate_Source_Data_AllColumns_Incl2019UpdateSources_new.csv", row.names = FALSE)
  #edited file in excel to remove some unnecessary columns; copying to Google Sheet for Teo 12/9/2021
