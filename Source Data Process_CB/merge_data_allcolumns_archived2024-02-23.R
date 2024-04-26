## Merge individual source datasets
## Created by Aleyna Loughran-Pierce
## Created 10/10/2020
## Last modified: 12/09/2021 by Kate Snyder - include all columns from each source file
## 3/7/2022 - add in Griesser 2017 data, female song (Odom 2014) data
## 3/8/2022 - add in BOW data
## 9/7/2022 - female song dataset
## 9/13/2023 - cornwallis data
## 10/26/2023 - Griesser et al 2023 PNAS data, also added OC data in initial source compilation, also AVONET data; sociality data binarized
## 2/23/2024 - saved current version as an archive

require(ape)
require(phytools)
library(readxl)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/")

birdtree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/birdzillatreeMaybeConsensus.nex")

ourdatabase <- read.csv("SongData_R_Update.csv", stringsAsFactors = FALSE)
ourdatabaserefs <- read.csv("SupplementDataRefs_Update.csv")
# colnames(ourdatabaserefs)[1] <- "BirdtreeSpecies"  # changed in csv
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
DunnData <- read.csv("Dunn supp data 977sp.csv", stringsAsFactors = FALSE)
GriesserData <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE) # Griesser and Suzuki 2016
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
#RiehlData <- read.csv("Riehl 2013 supp data columns lines_kts edited.csv", stringsAsFactors = FALSE)
RiehlData <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Riehl_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
#OdomFSData <- read.csv("FemaleSongData_OdomEtal2014_PresentAbsentSubset_BirdTreeNames.csv")
Griesser2017Data <- read.csv("Griesser2017 matched species names.csv", stringsAsFactors = FALSE)
Griesser2017Data$BirdtreeFormat[which(Griesser2017Data$Griesser2017_speciesnames == "finchyellowthighed")] <- "Pselliophorus_tibialis"
#BOWData <- read.csv("BOW Cooperative Breeding Data.csv")
BOWData <- read.csv("BOW CB data - merged w species with CB source discrepancies and SongRepvalue 2024-02-22_kts-added-bow-data.csv")
library(readxl)
Mikula2 = read_xlsx("/Users/kate/Downloads/Mikula_data 2.xlsx", 2) # Primary data used in analysis contrasting duetting species versus non-duetting species (species with non-singing and solo singing females combined). Duets: duetting present (1) or absent (0), Coop_breed: cooperative breeding, Territory1: year-round territoriality;  Territory2: seasonal territoriality; Territory3: weak territoriality; Social_bond1: long-term social bonds; Social_bond2: short-term social bonds; Social_bond3: solitary; NDVI: maximum normalized difference vegetation index.
colnames(Mikula2)[which(colnames(Mikula2) == "Duets")] <- "Mikula1_Duetting_vs_NoDuetting"
Mikula3 = read_xlsx("/Users/kate/Downloads/Mikula_data 2.xlsx", 3) # Primary data used in analysis contrasting duetting species versus species with non-singing females. Duets: duetting present (1) or female song absent (0)
colnames(Mikula3)[which(colnames(Mikula3) == "Duets")] <- "Mikula2_Duetting_vs_NoFS"
Mikula4 = read_xlsx("/Users/kate/Downloads/Mikula_data 2.xlsx", 4) # Primary data used in analysis contrasting species where females produce only solo songs versus duetting species. Solo_song: female solo song present (1) or duetting present (0)
colnames(Mikula4)[which(colnames(Mikula4) == "Solo_song")] <- "Mikula3_FemaleSoloSong_vs_Duetting"
Mikula5 = read_xlsx("/Users/kate/Downloads/Mikula_data 2.xlsx", 5) # Primary data used in analysis contrasting species where females produce only solo songs versus species with non-singing females. Solo_song: female solo song present (1) or no female song present (0)
colnames(Mikula5)[which(colnames(Mikula5) == "Solo_song")] <- "Mikula4_FemaleSoloSong_vs_NoFS"
DaleData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/plumage_scores_Dale et al 2016.csv") # Cooperative breeding was scored as absent (0), suspected (0.5), or present (1) primarily based on the data provided in ref. 63, including both the known and inferred data categories. For species in our data set not present in ref. 63 (Cockburn 2006), we used standard references, in particular ref. 18, to obtain the additional cooperative breeding scores.
#DaleCockburn = merge(DaleData, CockburnData, by.x = "TipLabel", by.y = "BirdtreeSpecies", all = T)
#DaleCockburnCoopSummary = DaleCockburn %>% group_by(Cooperative_breeding_ppca, KnownParentalCare, InferredParentalCare) %>% summarize(n=n())
#CockburnCoopSummary = CockburnData %>% group_by(KnownParentalCare, InferredParentalCare) %>% summarize(n=n())
#MikulaDale = merge(Mikula2, DaleData, by.x = "Sci_name", by.y = "TipLabel", all = T)
#MikulaDaleCoopSub = MikulaDale[which(!is.na(MikulaDale$Coop_breed) | !(is.na(MikulaDale$Cooperative_breeding_ppca))),]
#MikulaDaleCoopMikulaSub = MikulaDale[which(!is.na(MikulaDale$Coop_breed)),]


#### start merge ----
# merge Jetz + Downing
DowningSubset <- DowningData #[,c(1,2,4,11,12,13,14)]
JetzSubset <- JetzData #[, c(1,2,3,4,7,8)]
newdf <- merge(x = JetzSubset, y = DowningSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_Jetz", "_Downing"))
colnames(newdf)[which(colnames(newdf) == "Species")] <- "Species_Jetz"
colnames(newdf)[which(colnames(newdf) == "promiscuity....")] <- "Promiscuity_Downing"

# merge newdf + ourdatabase
newdf_wOurs1 <- merge(x = ourdatabaserefs, y = newdf, by = "BirdtreeSpecies", all = TRUE)

# merge newdf_wOurs
OCin = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/OCPaperData.csv")
newdf_wOurs <- merge(newdf_wOurs1, OCin, by.x = "BirdtreeSpecies", by.y = "BirdtreeFormat", all = TRUE, suffixes = c("", "_RobinsonSnyderCreanza"))


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

# merge new7 + Griesser2017
newdf8 <- merge(newdf7, Griesser2017Data, by.x =  "BirdtreeSpecies", by.y = "BirdtreeFormat", all = TRUE, suffixes = c("", "_Griesser2017"))
#colnames(newdf8)[which(colnames(newdf8) == "common.name")] <- "common.name_Griesser2017"

# merge newdf8 + OdomFS
#newdf9 <- merge(newdf8, OdomFSData, by.x =  "BirdtreeSpecies", by.y = "Latin_binomial", all = TRUE, suffixes = c("","_Odom"))
newdf9 = newdf8

# merge newdf9 and BOW data
newdf10 <- merge(newdf9, BOWData, by.x =  "BirdtreeSpecies", by.y = "SpeciesScientific.BOW", all = TRUE, suffixes = c("","_BOW"))

newdf11 <- newdf10[which(!newdf10$BirdtreeSpecies %in% c("(non-Oscine)","(non-Passerine)","(non-passerine)","(not_checked)","NA",NA)),]

newdf11$BirdtreeSpecies[which(!newdf11$BirdtreeSpecies %in% birdtree$tip.label)]

# merge Mikula data
newdf12 <- merge(newdf11, Mikula2[,c("Sci_name", "Mikula1_Duetting_vs_NoDuetting")], by.x = "BirdtreeSpecies", by.y = "Sci_name", all = T)
newdf13 <- merge(newdf12, Mikula3[,c("Sci_name", "Mikula2_Duetting_vs_NoFS")], by.x = "BirdtreeSpecies", by.y = "Sci_name", all = T)
newdf14 <- merge(newdf13, Mikula4[,c("Sci_name", "Mikula3_FemaleSoloSong_vs_Duetting")], by.x = "BirdtreeSpecies", by.y = "Sci_name", all = T)
newdf15 <- merge(newdf14, Mikula5, by.x = "BirdtreeSpecies", by.y = "Sci_name", all = T, suffixes = c("","_Mikula"))

# Merge Dale data
colnames(DaleData)[which(colnames(DaleData) == "Scientific_name")] <- "Scientific_name_Dale2015"
newdf16 <- merge(newdf15, DaleData, by.x = "BirdtreeSpecies", by.y = "TipLabel", all = T, suffixes = c("","_Dale"))

# Merge Remes data - does not contain female song or cooperative breeding per se, but has EPP and biparental cooperation
RemesData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data.csv")
RemesClimateData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data_climatic.csv")
sum(RemesData$Species_name %in% newdf14$BirdtreeSpecies)
RemesDataAll = merge(RemesData, RemesClimateData)
colnames(RemesDataAll)[which(colnames(RemesDataAll) == "Species_name")] <- "Species_Name_Remes"
newdf17 = merge(newdf16, RemesDataAll, by.x = "BirdtreeSpecies", by.y = "Species_Name_Remes", all = T, suffixes = c("","_Remes"))

# Merge Cornwallis et al 2017 Nature EcoEvo
cornwallisData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Unaltered from publication/Cornwallis 2017 supp table13.csv")
cornwallisDataNoDups = cornwallisData[which(!duplicated(cornwallisData$Species)),] # all species that have 2+ rows have the same breeding system in both/all rows, so duplicates removed 
newdf17b = merge(newdf17, cornwallisDataNoDups, by.x = "BirdtreeSpecies", by.y = "Species", all = T, suffixes = c("", "_Cornwallis"))

# Merge Griesser et al 2023 PNAS
Griesser2023 <- read.csv("Griesser et al 2023 pnas.csv")
Griesser2023$X = NULL
newdf17c = merge(newdf17b, Griesser2023, by.x = "BirdtreeSpecies", by.y = "tip_label", all = T, suffixes = c("", "_Griesser2023"))

# Merge Female Song aggregated data
FSdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-09-13_Female Song Data_GSheetDownload_withHighConfFS01.csv")
colnames(FSdata)[which(colnames(FSdata) == "Song_data_source_Webb..NOTE..appears.to.be.reversed...ones.marked..del.Hoyo..were.actually.from.Odom..ones.marked..del.Hoyo..actually.from.HBW.")] = "Song_data_source_Webb"
FSdata$Odom_FemaleSong[which(FSdata$Odom_FemaleSong == "")] <- NA
FSdata$Webb_FemaleSong[which(FSdata$Webb_FemaleSong == "")] <- NA
newdf17d = merge(newdf17c, FSdata, by.x = "BirdtreeSpecies", by.y = "BirdtreeSpecies", all = T, suffixes = c("", "_FSgsheet"))

# Merge AVONET
AVONET = read_excel("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/AVONET Supplementary dataset 1.xlsx", sheet = 4)
colnames(AVONET)[which(colnames(AVONET) == "Family3")] <- "Family3_BirdtreeMatchSpecies2"
colnames(AVONET)[which(colnames(AVONET) == "Order3")] <- "Order_AVONET"
AVONET$Species3 = str_replace(AVONET$Species3, " ", "_")
newdf17e = merge(newdf17d, AVONET, by.x = "BirdtreeSpecies", by.y = "Species3", all = T, suffixes = c("", "_AVONET"))

# remove non-birdtrees
sum(!newdf17e$BirdtreeSpecies %in% birdtree$tip.label)
newdf18 = newdf17e[which(newdf17e$BirdtreeSpecies %in% birdtree$tip.label),]

# remove duplicated rows
newdf19 = newdf18[which(!duplicated(newdf18$BirdtreeSpecies)),]

# write file - then use file in CoopBreed_species_summary.R
source("CoopBreed_species_summary.R")
#write.csv(newdf19, file = paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllCoopBreedColumns.csv"), row.names = FALSE)
#CoopBreed_species_summary(paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllCoopBreedColumns.csv"), songfile = "SongData_R_Update.csv", allcoop = TRUE, OCdatafile = TRUE)
write.csv(newdf19, file = paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllColumns.csv"), row.names = FALSE)
CoopBreed_species_summary(paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllColumns.csv"), songfile = "SongData_R_Update.csv", allcoop = TRUE, OCdatafile = FALSE)
cbSong = read.csv(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/", Sys.Date(),"CoopSong_R.csv")) # created with above line

# Finally, merge classified species with rest of relevant data
colnames(newdf19)
relevant_cols = c("BirdtreeSpecies", "di_OC", "tri_OC", "cont_OC", "O.C", "Region", "Cover", "MaleCareYg", "Cavity01", "Habitat", "Predictability..P.", "Constancy..C.", "Contingency..M.", "engage.in.misdirected.parental.care",	"family.time..in.days.", "social_system_assessment", "longevity",	"zoological_record_hits", "movements",	"food_specialisation",	"nest_type", "Reported.latitude",	"Reported.longitude",	"Corrected.latitude",	"Corrected.longitude", "Mean.precipitation",	"Variation.in.precipitation",	"Colwell.s.predictability.in.precipitation","Mean.temperature",	"Variation.in.temperature",	"Colwell.s.predictability.in.temperature", "time_fed",	"clutch_size",	"caretakers", "social_bonds", "sedentariness",	"insularity",	"colonial",	"grouping", "FemaleSong_Agg01",	"HighConfidence_FemaleSong", "Family3_BirdtreeMatchSpecies2", "Order_AVONET", "Habitat_AVONET",	"Habitat.Density",	"Migration",	"Trophic.Level",	"Trophic.Niche",	"Primary.Lifestyle",	"Min.Latitude",	"Max.Latitude",	"Centroid.Latitude",	"Centroid.Longitude",	"Range.Size")
relevant_cols %in% colnames(newdf19)
relevant_newdf19 = newdf19[,relevant_cols]
cbSongEtc = merge(cbSong, relevant_newdf19, by.x = "species", by.y = "BirdtreeSpecies", all = T)
write.csv(cbSongEtc, paste0(Sys.Date(),"_CoopBreed-FemaleSong-Song-Sociality_Data_R.csv"), row.names = FALSE)

# subset to Passeriformes and binarize sociality data
ordercounts = cbSongEtc %>% group_by(Order_AVONET) %>% count
allPass = cbSongEtc[which(cbSongEtc$Order_AVONET == "Passeriformes"),]

unique(allPass$colonial)
allPass$Griesser2023.Colonial01 = NA
allPass$Griesser2023.Colonial01[which(allPass$colonial == "acolonial")] <- 0
allPass$Griesser2023.Colonial01[which(allPass$colonial == "colonial")] <- 1

unique(allPass$caretakers)
allPass$Griesser2023.MoreThanTwoCaretakers = NA
allPass$Griesser2023.MoreThanTwoCaretakers[which(allPass$caretakers <= 2 )] <- 0
allPass$Griesser2023.MoreThanTwoCaretakers[which(allPass$caretakers > 2 )] <- 1

unique(allPass$social_bonds)
allPass %>% group_by(social_bonds) %>% count
allPass$Griesser2023.LongSocialBonds = NA
allPass$Griesser2023.LongSocialBonds[which(allPass$social_bonds %in% c("a-short", "b-season"))] <- 0
allPass$Griesser2023.LongSocialBonds[which(allPass$social_bonds == "c-long")] <- 1

unique(allPass$grouping)
allPass %>% group_by(grouping) %>% count
allPass$Griesser2023.GroupsLargerThanPair = NA
allPass$Griesser2023.GroupsLargerThanPair[which(allPass$grouping %in% c("asocial", "pair"))] <- 0
allPass$Griesser2023.GroupsLargerThanPair[which(allPass$grouping %in% c("small_groups", "large_groups"))] <- 1

write.csv(allPass, paste0(Sys.Date(),"_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"), row.names = FALSE)


allPass %>% group_by(Griesser2023.GroupsLargerThanPair) %>% count


#### scratch/misc ----
# moved up to above workflow 10/26/2023
FSdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-09-13_Female Song Data_GSheetDownload_withHighConfFS01.csv")
colnames(FSdata)[which(colnames(FSdata) == "Song_data_source_Webb..NOTE..appears.to.be.reversed...ones.marked..del.Hoyo..were.actually.from.Odom..ones.marked..del.Hoyo..actually.from.HBW.")] = "Song_data_source_Webb"

FSdata$BirdtreeSpecies[which(!FSdata$BirdtreeSpecies %in% birdtree$tip.label)]
FSdata = FSdata[which(FSdata$BirdtreeSpecies %in% birdtree$tip.label),]
FSdata$Odom_FemaleSong[which(FSdata$Odom_FemaleSong == "")] <- NA
FSdata$Webb_FemaleSong[which(FSdata$Webb_FemaleSong == "")] <- NA
cbSong$species[which(!cbSong$species %in% birdtree$tip.label)]
FSCBSong = merge(cbSong, FSdata, by.x = "species", by.y = "BirdtreeSpecies", all = T)
write.csv(FSCBSong, file = paste0(Sys.Date(),"_CoopBreed-FemaleSong-Song_Data_R.csv"), row.names = FALSE)



## 9/12/2023 - add in learning window data - actually ended up doing this above using the CoopBreed_species_summary.R function
## what to do with Philesturnus rufusater? Also have data that's not present in the all-sources song features file from 6/20/2023 - philesturnus rufusater, passerculus sandwichensis, parus palustris - not sure how to add yet
# read most recent file
dfIn = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-20_Aggregate_Source_Data_AllCoopBreedColumns.csv")
# read Fem Song data
dfFS = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-20_Female Song Data_GSheetDownload.csv")
dfFSR = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv")
FSbin = merge(dfFS, dfFSR[,c("species", "FemaleSong_Agg01", "HighConfidence_FemaleSong")], all.x = T, by.x = "BirdtreeSpecies", by.y = "species")
#write.csv(FSbin, "2023-06-20_Female Song Data_GSheetDownload_withHighConfFS01.csv")

# read OC data
OCin = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/OCPaperData.csv")
# Read song data
dfSong = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/SongData_R_Update.csv")

#colnames(OCin) = paste(colnames(OCin), )

dfOC = merge(dfIn, OCin, by.x = "BirdtreeSpecies", by.y = "BirdtreeFormat", all = T)
dfOCFS = merge(dfOC, dfFS, by = "BirdtreeSpecies", all = T)
dfOCFSsong = merge(dfOCFS, dfSong, by.x = "BirdtreeSpecies", by.y= "BirdtreeFormat", all = T, suffixes = c(".NotUpdated",""))
colnames(dfOCFSsong)
dfOCFSsong$BirdtreeSpecies[which(dfOCFSsong$Song.rep.final != dfOCFSsong$Song.rep.final.Update)]
dfOCFSsong$BirdtreeSpecies[which(dfOCFSsong$Interval.final != dfOCFSsong$Interval.final.Update)]
dfOCFSsong$BirdtreeSpecies[which(dfOCFSsong$Duration.final != dfOCFSsong$Duration.final.Update)]

write.csv()


#### Additional (non-Cooperative Breeding) Sources ----
LislevandData <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Lislevand et al 2007 supp data.csv")
SpottiswoodeData <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Spottiswoode and Moller 2004_data_noWeirdCharacters.csv")
DatabaseUpdate2019 <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Song Database Update 2019_LastEdited_20200408_539PM_LislevandSpottiswoodeClementsOnly.csv")
Clements2021 <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/eBird-Clements-v2021-integrated-checklist-August-2021_kts_speciesonly.csv")
Clements2021subset <- Clements2021[,c(6,5,7,8,9,10,11)]
colnames(DatabaseUpdate2019)[which(colnames(DatabaseUpdate2019) == "BirdtreeFormat")] <- "BirdtreeSpecies"
DatabaseUpdate2019_Clements <- DatabaseUpdate2019[,c(1:7,9,10,12,13,14,15,16,17)]
DatabaseUpdate2019_LisSpo <- DatabaseUpdate2019[,c(1,18:64)]

#merge DatabaseUpdate2019_Clements and newdf7
# newdf8 <- merge(DatabaseUpdate2019_Clements, newdf7, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_ClementsGuide", ""))
# newdf9 <- merge(newdf8, DatabaseUpdate2019_LisSpo, by = "BirdtreeSpecies", all = TRUE)
# newdf10 <- merge(Clements2021subset, newdf9, by.x = "Scientific", by.y = "BirdtreeSpecies", all.x = FALSE, all.y = TRUE, suffixes = c("ClementsGuide2021",))
#write.csv(newdf10, file = "Aggregate_Source_Data_AllColumns_Incl2019UpdateSources_new.csv", row.names = FALSE)
  #edited file in excel to remove some unnecessary columns; copying to Google Sheet for Teo 12/9/2021


#### Female song ----
# most recent aggregated FS data file
FSdf = read.csv("2023-09-13_Female Song Data_GSheetDownload_withHighConfFS01.csv")

OdomFSData <- read.csv("FemaleSongData_OdomEtal2014_PresentAbsentSubset_BirdTreeNames.csv")
OdomFSData$Latin_binomial <- str_replace(OdomFSData$Latin_binomial, " ", "_")
OdomFSDataNonBT=read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/FemaleSongData_OdomEtal2014_PresentAbsentSubset.csv")
Odomdf = cbind(OdomFSDataNonBT$Latin_binomial, OdomFSData)
colnames(Odomdf)[1] = "OdomOGSpeciesName"
colnames(Odomdf)[2] = "BirdtreeSpecies"
colnames(Odomdf)[6] = "Odom_FemaleSong"
Odomdf=Odomdf[,c(1,2,3,6)]

WebbFSData <- read.csv("/Users/kate/Documents/Creanza Lab/Paper PDFs/Webb et al 2016 Female Song Plumage Data.csv")
colnames(WebbFSData)[which(colnames(WebbFSData) == "Female_song_score")] = "Webb_FemaleSong"
WebbFSData = WebbFSData[1:1314,]
FSData=merge(Odomdf, WebbFSData, by.x = "BirdtreeSpecies", by.y="TipLabel", all=T)
#write.csv(FSData, "FemaleSong_OdomWebb.csv")

GaraFSData <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Garamszegi et al 2006 - Appendix 1 - Female Song Data .csv")
GaraFSData <- GaraFSData[which(GaraFSData$Female.song %in% c(0,2)),c("BirdtreeSpecies","Female.song")]
GaraFSData$Female.song[which(GaraFSData$Female.song==2)] = "Evidence for FS"
GaraFSData$Female.song[which(GaraFSData$Female.song==0)] = "Evidence against FS"
colnames(GaraFSData)[2] = "Garamszegi_FemaleSong"
FSData = merge(FSData, GaraFSData, all = T)
#write.csv(FSData, "FemaleSong_OdomWebbGaramszegi.csv")
#


# 5/4/2023 - Find the right file and add female song stuff to it
library(phytools)
birdtree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/birdzillatreeMaybeConsensus.nex")
passertree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/2020-10-11ConsensusPasserineTreeHack100.nex")
# df = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Aggregate_Source_Data_AllColumns.csv")
# df2 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Aggregate_Source_Data_unifiedCommonNames.csv")
# df3 = read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Aggregate_Source_Data_AllColumns_Incl2019UpdateSources.csv")
# df2Common = df2[,c("BirdtreeSpecies","CommonName_aggregate")]
# cbdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2022-03-10_working_coop_breed.csv")
# colnames(df3)
# is.na(df3[,1])
# FSdataOdom = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/FemaleSong_OdomWebb.csv")

Aggdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_Aggregate_Source_Data_AllColumns.csv")
#colnames(Aggdf)[which(colnames(Aggdf) == "BirdtreeSpecies")] <- "species"
#cbdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-05-30CoopSong_All.csv")
cbdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_working_coop_breed.csv")

FSdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_Female Song Data_GSheetDownload.csv")
#colnames(FSdata) %in% colnames(Aggdf)
#colnames(Aggdf) %in% colnames(FSdata)
#sum(colnames(Aggdf) %in% colnames(cbdf))
#AggdfNewFS = Aggdf[which(!is.na(Aggdf$Mikula1_Duetting_vs_NoDuetting)),c("BirdtreeSpecies", "Mikula1_Duetting_vs_NoDuetting", "Mikula2_Duetting_vs_NoFS", "Mikula3_FemaleSoloSong_vs_Duetting", "Mikula4_FemaleSoloSong_vs_NoFS")]
#!AggdfNewFS$BirdtreeSpecies %in% FSdata$BirdtreeSpecies

FSdata$BirdtreeSpecies[which(!FSdata$BirdtreeSpecies %in% birdtree$tip.label)]
#FSdata = FSdata[,c("BirdtreeSpecies","Odom_FemaleSong","Webb_FemaleSong","FemaleSong_Aggregated")]
FSdata$Odom_FemaleSong[which(FSdata$Odom_FemaleSong == "")] <- NA
FSdata$Webb_FemaleSong[which(FSdata$Webb_FemaleSong == "")] <- NA
cbdf$species[which(!cbdf$species %in% birdtree$tip.label)]

orderedNewFS = merge(FSdata, AggdfNewFS, all = T)
#write.csv(orderedNewFS, file = "2023-06-16_Female Song Data_GSheet2023-05-30_addedMikula.csv", row.names = F)  #pasted this new file back into the Female Song Data Google Sheet https://docs.google.com/spreadsheets/d/1kEijryIyLgJ-Co40pVKkxqAPUDEoWGkCd6YY4tVTMVM/edit#gid=0; then manually changed FemaleSong_Aggregate column to include Mikula classifications on 6/20/2023 KTS
#sum(!orderedNewFS$BirdtreeSpecies %in% birdtree$tip.label)
#orderedNewFS %>% group_by(FemaleSong_Aggregated, Mikula1_Duetting_vs_NoDuetting, Mikula2_Duetting_vs_NoFS, Mikula3_FemaleSoloSong_vs_Duetting, Mikula4_FemaleSoloSong_vs_NoFS) %>% summarize(n=n())
newFS = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-20_Female Song Data_GSheetDownload.csv")
newFS %>% group_by(Mikula_SummaryClass, FemaleSong_Aggregated) %>% summarize(n=n())
#ReviewFS = newFS[which(newFS$Mikula_SummaryClass == "Absent Female Song" & newFS$FemaleSong_Aggregated == "Present" | (newFS$Mikula_SummaryClass == "Duetting" & newFS$FemaleSong_Aggregated == "Absent") | newFS$Mikula_SummaryClass == "Female Solo Song" & newFS$FemaleSong_Aggregated == "Absent"),]
#AddFS = newFS[which(newFS$Mikula_SummaryClass == "Absent Female Song" & is.na(newFS$FemaleSong_Aggregated) | newFS$Mikula_SummaryClass == "Female Solo Song" & is.na(newFS$FemaleSong_Aggregated)),]

cbFS = merge(cbdf, FSdata, by.x = "species", by.y = "BirdtreeSpecies", all=T)
cbFSAgg = merge(cbFS, Aggdf, by.x = "species", by.y = "BirdtreeSpecies", all=T)
#alldf = merge(cbdf, FSdata, by.x = "species", by.y = "BirdtreeSpecies", all = T)
#alldf = merge(df3, FSdata, by.x = "Scientific", by.y = "BirdtreeSpecies", all = T)


write.csv(cbFSAgg, paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/", Sys.Date(), "_CoopSongFS_AllSourceColumns.csv"), row.names = F)

colnames(cbFSAgg)

# did this early morning 6/14/23
alldfR = cbFSAgg[,c(1:30, 33, 34, 37, 40, 43, 46, 49, 52, 58, 60:71, 80:85, 87:104, 108:116, 118, 124:126, 128:139, 140, 145:153, 165, 168, 172:181, 183:198, 203, 205:211)]. 
write.csv(alldfR, paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/", Sys.Date(), "_CoopSongFS_RColumns.csv"), row.names = F)



# below here just exploring data
colnames(alldf)
alldf %>% group_by(MeanCoopTie2Noncoop, Webb_FemaleSong) %>% summarize(n=n())
passerdf = alldf[which(alldf$species %in% passertree$tip.label),]
passerdf %>% group_by(MeanCoopTie2Noncoop, FemaleSong_Aggregated) %>% summarize(n=n())
dfsub = passerdf[which(!is.na(passerdf$MeanCoopTie2Noncoop) & !is.na(passerdf$FemaleSong_Aggregated)),]
dfsub %>% group_by(MeanCoopTie2Noncoop, FemaleSong_Aggregated) %>% summarize(n=n())
chisq.test(dfsub$MeanCoopTie2Noncoop, dfsub$FemaleSong_Aggregated)

alldf = alldf[which(alldf$Scientific %in% birdtree$tip.label),]
cbdf$species[which(!cbdf$species %in% alldf$Scientific)]
alldfCoop = merge(alldf, cbdf, by.x = "Scientific", by.y = "species", all = T)
alldf = alldfCoop

TestColumns = c("Syllable.rep.final", "Song.rep.final", "MonogamyOrNot", "EPP10threshold", "MeanCoopTie2Coop")
outdf = set.seed(10)
for (i in 1:length(alldfCoop$Scientific)) {
  temprow = alldfCoop[i,]
  numNotNA= sum(!is.na(temprow[,TestColumns]))
  species = temprow$Scientific
  commonname = temprow$CommonName
  FemaleSong_Aggregated= temprow$FemaleSong_Aggregated
  tempnewrow = c(species, commonname, numNotNA, FemaleSong_Aggregated)
  outdf = rbind(outdf, tempnewrow)
}
colnames(outdf) = c("species", "CommonName", "numNotNA", "FemaleSong_Aggregated")
outdf = as.data.frame(outdf)
sum(outdf$numNotNA > 4)
speciesToGet = outdf[which(outdf$numNotNA > 3 & is.na(outdf$FemaleSong_Aggregated)),c("species","CommonName")]

dfSubsetSpeciesToGet = alldf[which(alldf$Scientific %in% speciesToGet$species),]

alldf$Scientific[which(!alldf$Scientific %in% birdtree$tip.label)]


alldf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_CoopSongFS_RColumns.csv")
alldf %>% group_by(FemaleSong_Aggregated, Mikula_FemaleSoloSong.vs.NoFS, Mikula_Duetting) %>% summarize(n=n())

cbdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14_working_coop_breed.csv")
cbdf %>% group_by(MeanCoopTie2Noncoop, MikulaCoop, DaleCoop) %>% summarize(n=n())

cbdf[which(cbdf$MeanCoopTie2Coop == 0 & cbdf$DaleCoop == 1),]
cbdf[which(cbdf$MeanCoopTie2Coop == 1 & cbdf$DaleCoop == 0),]
cbdf[which(is.na(cbdf$MeanCoopTie2Coop) & cbdf$DaleCoop == 0),]

coopsongdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-14CoopSong_All.csv")


