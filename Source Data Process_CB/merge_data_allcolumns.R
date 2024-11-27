## Merge individual source datasets
## Created by Kate Snyder
## Created 10/10/2020 - split from archived code 2/24/2024
## Last modified: 12/09/2021 by Kate Snyder - include all columns from each source file
## 3/7/2022 - add in Griesser 2017 data, female song (Odom 2014) data
## 3/8/2022 - add in BOW data
## 9/7/2022 - female song dataset
## 9/13/2023 - cornwallis data
## 10/26/2023 - Griesser et al 2023 PNAS data, also added OC data in initial source compilation, also AVONET data; sociality data binarized
## 2/23/2024 - removed unnecessary code/lines
## 2/24/2024 - removed merging Mikula data because the duetting info is in the FS Gsheet
## 5/13/2024 - changed ourdatabaserefs file from "SupplementDataRefs_Update.csv" to "SupplementDataRefs_Update_2024-05-13.csv"; made output files use current date

require(ape)
require(phytools)
library(readxl)

# setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/")

birdtree = read.nexus("birdzillatreeForTipNames.nex")

ourdatabase <- read.csv("SongData_R_Update.csv", stringsAsFactors = FALSE)
ourdatabaserefs <- read.csv("SupplementDataRefs_Update_2024-05-13.csv")
colnames(ourdatabaserefs)[which(colnames(ourdatabaserefs) == "Family")] <- "Family_SnyderCreanza2019"
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
DunnData <- read.csv("Dunn supp data 977sp.csv", stringsAsFactors = FALSE)
Griesser2016Data <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE)
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
RiehlData <- read.csv("Riehl_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
RiehlData$X = NULL
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
Griesser2017Data <- read.csv("Griesser2017 matched species names.csv", stringsAsFactors = FALSE)
Griesser2017Data$BirdtreeFormat[which(Griesser2017Data$Griesser2017_speciesnames == "finchyellowthighed")] <- "Pselliophorus_tibialis"
DaleData = read.csv("plumage_scores_Dale et al 2016.csv") # Cooperative breeding was scored as absent (0), suspected (0.5), or present (1) primarily based on the data provided in ref. 63, including both the known and inferred data categories. For species in our data set not present in ref. 63 (Cockburn 2006), we used standard references, in particular ref. 18, to obtain the additional cooperative breeding scores.
BOWData <- read.csv("BOW CB data - merged w species with CB source discrepancies and SongRepvalue 2024-02-22_kts-added-bow-data.csv")
BOWData = BOWData[,c("species","Cooperative.BOW","BOW.Quote")]


#Append author of source paper to each column
colnames(BiagoliniData) <- ifelse(colnames(BiagoliniData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(BiagoliniData), "Biagolini2017", sep = "_"))
colnames(CockburnData) <- ifelse(colnames(CockburnData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(CockburnData), "Cockburn2006", sep = "_"))
colnames(DowningData) <- ifelse(colnames(DowningData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(DowningData), "Downing2015", sep = "_"))
colnames(DunnData) <- ifelse(colnames(DunnData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(DunnData), "Dunn2015", sep = "_"))
colnames(Griesser2016Data) <- ifelse(colnames(Griesser2016Data) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(Griesser2016Data), "Griesser2016", sep = "_"))
colnames(JetzData) <- ifelse(colnames(JetzData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(JetzData), "Jetz2011", sep = "_"))
colnames(RiehlData) <- ifelse(colnames(RiehlData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(RiehlData), "Riehl2013", sep = "_"))
colnames(RubensteinLovetteData) <- ifelse(colnames(RubensteinLovetteData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(RubensteinLovetteData), "RubensteinLovette2007", sep = "_"))
colnames(Griesser2017Data) <- ifelse(colnames(Griesser2017Data) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(Griesser2017Data), "Griesser2017", sep = "_"))
colnames(DaleData) <- ifelse(colnames(DaleData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(DaleData), "Dale2015", sep = "_"))


#### start merge ----
# merge Snyder and Creanza 2019 data with Robinson, Snyder, Creanza 2019 data
OCin = read.csv("OCPaperData.csv")
OCin$Family = NULL
newdf <- merge(ourdatabaserefs, OCin, by.x = "BirdtreeSpecies", by.y = "BirdtreeFormat", all = TRUE, suffixes = c("", "_RobinsonSnyderCreanza"))

# merge Jetz + Downing
DowningSubset <- DowningData 
JetzSubset <- JetzData 
JetzDowning <- merge(x = JetzSubset, y = DowningSubset, by.x = "BirdtreeSpecies", all = TRUE, suffixes = c("_Jetz", "_Downing"))

# merge newdf + JetzDowning
newdf_wOurs1 <- merge(x = newdf, y = JetzDowning, by = "BirdtreeSpecies", all = TRUE)

# merge newdf_wOurs and Cockburn
CockburnSubset <- CockburnData
newdf2 <- merge(newdf_wOurs1, CockburnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Cockburn"))

# merge newdf2 and Dunn
DunnSubset <- DunnData 
newdf3 <- merge(newdf2, DunnSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Dunn"))

# merge newdf3 and Biagolini
BiagoliniSubset <- BiagoliniData 
newdf4 <- merge(newdf3, BiagoliniSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Biagolini"))

# merge newdf4 and Riehl
RiehlSubset <- RiehlData
newdf5 <- merge(newdf4, RiehlSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Riehl"))

# merge newdf5 and Rubenstein
RubensteinSubset <- RubensteinLovetteData 
newdf6 <- merge(newdf5, RubensteinSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Rubenstein"))

# merge newdf6 and Griesser2016
GriesserSubset <- Griesser2016Data 
newdf7 <- merge(newdf6, GriesserSubset, by = "BirdtreeSpecies", all = TRUE, suffixes = c("", "_Griesser2016"))

# merge new7 + Griesser2017
newdf8 <- merge(newdf7, Griesser2017Data, by.x =  "BirdtreeSpecies", by.y = "BirdtreeFormat_Griesser2017", all = TRUE, suffixes = c("", "_Griesser2017"))
newdf9 = newdf8 # vestigial code

# merge newdf9 and BOW data
newdf10 <- merge(newdf9, BOWData, by.x =  "BirdtreeSpecies", by.y = "species", all = TRUE, suffixes = c("","_BOW"))

# remove some species from Griesser et al 2017 that we didn't endeavor to find the BirdTree match to
newdf11 <- newdf10[which(!newdf10$BirdtreeSpecies %in% c("(non-Oscine)","(non-Passerine)","(non-passerine)","(not_checked)","NA",NA)),]

# check for species that do not have a match in BirdTree
newdf11$BirdtreeSpecies[which(!newdf11$BirdtreeSpecies %in% birdtree$tip.label)]
newdf15 <- newdf11 # vestigial

# Merge Dale data
newdf16 <- merge(newdf15, DaleData, by.x = "BirdtreeSpecies", by.y = "TipLabel_Dale2015", all = T, suffixes = c("","_Dale"))

# Merge Remes data - does not contain female song or cooperative breeding per se, but has EPP and biparental cooperation
#RemesData = read.csv("Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data.csv")
#RemesClimateData = read.csv("Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data_climatic.csv")
#RemesDataAll = merge(RemesData, RemesClimateData)
#colnames(RemesDataAll) <- paste(colnames(RemesDataAll), "Remes2015", sep = "_")
#newdf17 = merge(newdf16, RemesDataAll, by.x = "BirdtreeSpecies", by.y = "Species_name_Remes2015", all = T, suffixes = c("","_Remes"))

newdf17 = newdf16

# Merge Cornwallis et al 2017 Nature EcoEvo
cornwallisData = read.csv("Unaltered from publication/Cornwallis 2017 supp table13.csv")
cornwallisDataNoDups = cornwallisData[which(!duplicated(cornwallisData$Species)),] # all species that have 2+ rows have the same breeding system in both/all rows, so duplicates removed 
colnames(cornwallisDataNoDups) <- paste(colnames(cornwallisDataNoDups), "Cornwallis2017", sep = "_")
newdf17b = merge(newdf17, cornwallisDataNoDups, by.x = "BirdtreeSpecies", by.y = "Species_Cornwallis2017", all = T, suffixes = c("", "_Cornwallis"))

# Merge Griesser et al 2023 PNAS
Griesser2023 <- read.csv("Griesser et al 2023 pnas.csv")
Griesser2023$X = NULL
colnames(Griesser2023) <- paste(colnames(Griesser2023), "Griesser2023", sep = "_")
newdf17c = merge(newdf17b, Griesser2023, by.x = "BirdtreeSpecies", by.y = "tip_label_Griesser2023", all = T, suffixes = c("", "_Griesser2023"))

# Merge Female Song aggregated data
FSdata = read.csv("2023-09-13_Female Song Data_GSheetDownload_withHighConfFS01.csv")
colnames(FSdata)[which(colnames(FSdata) == "Song_data_source_Webb..NOTE..appears.to.be.reversed...ones.marked..del.Hoyo..were.actually.from.Odom..ones.marked..del.Hoyo..actually.from.HBW.")] = "Song_data_source_Webb"
FSdata$Odom_FemaleSong[which(FSdata$Odom_FemaleSong == "")] <- NA
FSdata$Webb_FemaleSong[which(FSdata$Webb_FemaleSong == "")] <- NA
FSdata[,c( "Flagged.for.review" , "NC.notes", "NC.assessment")] = NULL
newdf17d = merge(newdf17c, FSdata, by.x = "BirdtreeSpecies", by.y = "BirdtreeSpecies", all = T, suffixes = c("", "_FSgsheet"))

# Merge AVONET
AVONET = read_excel("AVONET Supplementary dataset 1.xlsx", sheet = 4)
colnames(AVONET)[which(colnames(AVONET) == "Family3")] <- "Family3_BirdtreeMatchSpecies2"
colnames(AVONET)[which(colnames(AVONET) == "Order3")] <- "Order"
AVONET$Species3 = str_replace(AVONET$Species3, " ", "_")
colnames(AVONET) <- paste(colnames(AVONET), "AVONET", sep = "_")
newdf17e = merge(newdf17d, AVONET, by.x = "BirdtreeSpecies", by.y = "Species3_AVONET", all = T, suffixes = c("", "_AVONET"))

# remove non-birdtree-matched species
sum(!newdf17e$BirdtreeSpecies %in% birdtree$tip.label)
newdf17e$BirdtreeSpecies[which(!newdf17e$BirdtreeSpecies %in% birdtree$tip.label)]
newdf18 = newdf17e[which(newdf17e$BirdtreeSpecies %in% birdtree$tip.label),]

# remove duplicated species - due to different subspecies in Griesser et al 2017 not differentiated in BirdTree; subspecies do not differ in categorical variables used in this study
newdf19 = newdf18[which(!duplicated(newdf18$BirdtreeSpecies)),]

# write file - then use file in CoopBreed_species_summary.R
write.csv(newdf19, file = paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllColumns.csv"), row.names = FALSE)

#### run CoopBreed source aggregation summary ----
source("CoopBreed_species_summary.R")
CoopBreed_species_summary(paste0(Sys.Date(),"_Aggregate_CBSource_Data_AllColumns.csv"), songfile = "SongData_R_Update.csv")
# pull out all species with a source discrepancy, write to csv
CoopSubset = CoopSourcesdf[which(CoopSourcesdf$SourceDiscrepancy == 1 & CoopSourcesdf$Order_AVONET == "Passeriformes"),colnames(CoopSourcesdf)[which(!str_detect(colnames(CoopSourcesdf), "Griesser2023"))]]
CoopSubset = CoopSubset[,1:63]
#write.csv(CoopSubset, "2024-02-24_CoopClassesWSourceColumns_subsetPasserineSourceDiscrepancy.csv", row.names = F)

# opened ^ in excel and manually added new column "HighConfidence_Coop"
CoopSubset = read.csv("2024-02-24_CoopClassesWSourceColumns_subsetPasserineSourceDiscrepancy_HighConfColumn.csv")
# open big file with all refs and merge in HighConfidence_Coop column
ourdf = read.csv(paste0(Sys.Date(),"_CoopClassesWSourceColumns.csv"))
fulldf = merge(ourdf, CoopSubset[,c("species", "HighConfidence_Coop")], by = "species", all =T)
# in line with disputed species, if Cockburn "KnownParentalCare" is "Suspected" and it is the only cooperative breeding source (which would have meant it was classified as "Noncoop"), will give value "NA" for cooperative breeding. 
fulldf$HighConfidence_Coop[which(fulldf$SourceDiscrepancy == 0)] <- fulldf$MeanCoopTie2Noncoop[which(fulldf$SourceDiscrepancy == 0)] # give everything not disputed its classification
# put NAs in where the only source is Cockburn "Suspected"
fulldf$HighConfidence_Coop[which(fulldf$SourceDiscrepancy == 0 & fulldf$numSourcesNonCoop == 1 & fulldf$KnownParentalCare_Cockburn2006 == "Suspected")] <- NA

write.csv(fulldf, paste0(Sys.Date(), "_CoopClassesWSourceColumns_HighConfCoopColumn.csv"), row.names = F)

fulldf %>% group_by(MeanCoopTie2Noncoop, FemaleSong_Agg01) %>% count
fulldf %>% group_by(HighConfidence_Coop, FemaleSong_Agg01) %>% count

# merge fulldf with rest of columns in Song + binarized sociality file from CoopBreed_species_summary()
CoopSongSoc = read.csv(paste0(Sys.Date(), "_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"))
SongSocCols = colnames(CoopSongSoc)[which(!colnames(CoopSongSoc) %in% colnames(fulldf))]
CoopSongSocPasser = merge(fulldf, CoopSongSoc[,c("species", SongSocCols)], by = "species", all.y = T)
write.csv(CoopSongSocPasser, paste0(Sys.Date(), "_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"), row.names = F)
