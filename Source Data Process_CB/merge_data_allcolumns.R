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
## 2/28/2025 - at end of script, added AVONET data addition directly to Data_R_Passerine_withTobias.csv
## 6/6/2025 - added Tobias Territoriality data with code originally written/performed in "scratch MCMCglmm 2.R"; AVONET data was originally the main thing added in this version compared to the version of "merge_data_allcolumns.R" in the Git repo - changed it here to reflect that the species data files made were the Oscine subset, not Passerine; 


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


#### add Tobias et al 2016 data to 1st submission dataset ----

newdata = "Data_R.csv"
newdf = read.csv(newdata)
rownames(newdf) <- newdf$species
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
tree = read.nexus(treefile)
sum(!newdf$species %in% tree$tip.label)
newdfsub = newdf[which(newdf$species %in% tree$tip.label),]

require(readxl)
require(stringr)
require(dplyr)
tobias = read_excel("Unaltered from publication/Tobias et al 2016 Territoriality Communal Signalling - Data Sheet 3.xlsx")
tobias$BirdTree = str_replace(tobias$Species, " ", "_")
#write.csv(tobias, "Tobias et al 2016 Territoriality Communal Signalling DataSheet3_underscoredSpecies.csv", row.names = F)
tobias = read.csv("Tobias et al 2016 Territoriality Communal Signalling DataSheet3_underscoredSpecies.csv")
tobias$TobiasSpeciesUnderscored <- tobias$BirdTree
inconsistentNames = read.csv("inconsistent_species_names_InclTobias.csv") ## see below for process of getting this list
inconsistentNames <- inconsistentNames[which(!is.na(inconsistentNames$species_in_birdtree)),]
for (i in 1:length(tobias$BirdTree)) {
  tobiasSpecies = tobias$BirdTree[i]
  if (!tobiasSpecies %in% newdf$species) {
    if (tobiasSpecies %in% inconsistentNames$in_database) {
      # if (length(inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database == tobiasSpecies)]) > 1) {
      #   print(paste(tobiasSpecies, "getting replaced with", inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database == tobiasSpecies)])) 
      #   print(paste("Actually getting replaced with", na.omit(inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database == tobiasSpecies)])))
      # }
      tobias$BirdTree[i] <- na.omit(inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database == tobiasSpecies)])
      print(paste(na.omit(inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database == tobiasSpecies)])))
    }
  }
}
#write.csv(tobias, "Tobias et al 2016 Territoriality Communal Signalling DataSheet3_FSxCB species matchBirdtree_withdups.csv", row.names = F)
# manually changed the BirdTree name of the duplicated species to NA

tobias <- read.csv("Tobias et al 2016 Territoriality Communal Signalling DataSheet3_FSxCB species matchBirdtree_dupsNA.csv")
newdf_tobias = merge(newdfsub, tobias, by.x = "species", by.y = "BirdTree", all.x = T)
newdf_tobias$Territory <- as.factor(newdf_tobias$Territory)
newdf_tobias$Territory_12vs3 <- NA
newdf_tobias$Territory_12vs3[which(newdf_tobias$Territory %in% c(1,2))] <- 0
newdf_tobias$Territory_12vs3[which(newdf_tobias$Territory %in% c(3))] <- 1
newdf_tobias$Territory_1vs23 <- NA
newdf_tobias$Territory_1vs23[which(newdf_tobias$Territory %in% c(1))] <- 0
newdf_tobias$Territory_1vs23[which(newdf_tobias$Territory %in% c(2,3))] <- 1
table(newdf_tobias$Territory_12vs3, newdf_tobias$Territory_1vs23)
#write.csv(newdf_tobias,"Data_R_Oscine_withTobias.csv", row.names = F)


## From above - This was what I did to get the Tobias species that weren't in the female song + coop breed set
# sum(!newdf$species %in% tobias$BirdTree)
# newdf[which(!newdf$species %in% tobias$BirdTree),] %>% group_by(HighConfidence_Coop, FemaleSong_Agg01) %>% count # gonna need to figure these out
# newdf[which(!newdf$species %in% tobias$BirdTree & !is.na(newdf$HighConfidence_Coop) & !is.na(newdf$FemaleSong_Agg01)),"species"]
# tobiasUnmatchedPasserines = tobias[which(!tobias$BirdTree %in% newdf$species & tobias$Order == "Passeriformes"),]
# 
# inconsistentNames = read.csv('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data Processing/inconsistent_species_names_BirdTree.csv')
# missingFSxCBspecies <- newdf_clean[which(!newdf_clean$species %in% tobias$BirdTree), c("species")]
# tobiasPasserines = tobias[which(tobias$Order == "Passeriformes"),]
# tobiasUnmatchedPasserines = tobias[which(!tobias$BirdTree %in% newdf$species & tobias$Order == "Passeriformes"),]
# sum(tobiasUnmatchedPasserines$BirdTree %in% inconsistentNames$in_database)
# missingFSxCBspecies[missingFSxCBspecies %in% inconsistentNames$species_in_birdtree[which(inconsistentNames$in_database %in% tobiasUnmatchedPasserines$BirdTree)]]
# 
# MissingFromTobias = as.data.frame(cbind("", missingFSxCBspecies, "", "Tobias"))
# colnames(MissingFromTobias) <- colnames(inconsistentNames)
# inconsistentNamesWTobiasMissing <- rbind(inconsistentNames, MissingFromTobias)
# duplicated(inconsistentNamesWTobiasMissing$species_in_birdtree)
# #write.csv(inconsistentNamesWTobiasMissing, "inconsistent_species_names_InclTobias.   .csv", row.names = F)

#### Add AVONET data ----
datapath = "Data_R_Oscine_withTobias.csv" 
data = read.csv(datapath)
data = data[,which(!colnames(data) %in% c("Family3_BirdtreeMatchSpecies2_AVONET", "Order_AVONET"))] # remove the two AVONET columns previously used

require(readxl)
require(stringr)
AVONET = read_excel("Unaltered from publication/AVONET Supplementary dataset 1.xlsx", sheet = 4)
colnames(AVONET)[which(colnames(AVONET) == "Family3")] <- "Family3_BirdtreeMatchSpecies2"
colnames(AVONET)[which(colnames(AVONET) == "Order3")] <- "Order"
AVONET$Species3 = str_replace(AVONET$Species3, " ", "_")
colnames(AVONET) <- paste(colnames(AVONET), "AVONET", sep = "_")
newdfAVONET = merge(data, AVONET, by.x = "species", by.y = "Species3_AVONET", all.x = T)
colnames(newdfAVONET)
write.csv(newdfAVONET, "Data_R_Oscine_withTobias_AVONET.csv", row.names = F)

#### Add special territoriality classifications performed by KTS and NC from BOW and literature searches ----
 # Initially performed in "integrate bow categorization data.R"
data = read.csv("Data_R_Oscine_withTobias_AVONET.csv")
ktsnotes = read.csv('Territory2 species with KTS notes - Sheet1_20250509.csv') # already has BirdTree matches, can skip to "Basic analyses"
ktsnotes = ktsnotes[which(!ktsnotes$SPECIES %in% c("Buru Friarbird", "Chihuahuan Meadowlark", "Lynes's Cisticola", "Chestnut-capped Warbler")),] # drop duplicate species (species that are split in BOW, but not recognized in BirdTree)
ktsnotes1 = rename(ktsnotes, PermissiveExclusive = Exclusive.vs.Not.exclusive.territoriality..should.be.Strong...Exclusive.and.Weak.Grouped.except.for.colonial.and.group.defense.species..which.would.become.permissive., WeakStrong = Together.assessment..Weak.Strong., Confidence = Confidence.together)
mergeddata = merge(data, ktsnotes1[,c("BirdtreeMatch", "PermissiveExclusive", "WeakStrong", "Confidence")], by.x = "species", by.y = "BirdtreeMatch", all.x = T)

# Make special classifications binary 
mergeddata$TerritorialityWeakVsStrong = NA
mergeddata$TerritorialityWeakVsStrong[which(mergeddata$Territory == "1")] <- "0"
mergeddata$TerritorialityWeakVsStrong[which(mergeddata$Territory == "3")] <- "1"
mergeddata$TerritorialityWeakVsStrong[which(mergeddata$WeakStrong == "Weak")] <- "0"
mergeddata$TerritorialityWeakVsStrong[which(mergeddata$WeakStrong == "Strong")] <- "1"

mergeddata$TerritorialityWeakVsStrongHighConf = NA
mergeddata$TerritorialityWeakVsStrongHighConf[which(mergeddata$Territory == "1")] <- "0"
mergeddata$TerritorialityWeakVsStrongHighConf[which(mergeddata$Territory == "3")] <- "1"
mergeddata$TerritorialityWeakVsStrongHighConf[which(mergeddata$WeakStrong == "Weak" & mergeddata$Confidence %in% c("Medium", "High"))] <- "0"
mergeddata$TerritorialityWeakVsStrongHighConf[which(mergeddata$WeakStrong == "Strong" & mergeddata$Confidence %in% c("Medium", "High"))] <- "1"

mergeddata$TerritorialityPermissiveExclusive = NA
mergeddata$TerritorialityPermissiveExclusive[which(mergeddata$Territory == "1")] <- "0"
mergeddata$TerritorialityPermissiveExclusive[which(mergeddata$Territory == "3")] <- "1"
mergeddata$TerritorialityPermissiveExclusive[which(mergeddata$PermissiveExclusive == "Permissive")] <- "0"
mergeddata$TerritorialityPermissiveExclusive[which(mergeddata$PermissiveExclusive == "Exclusive")] <- "1"

mergeddata$TerritorialityPermissiveExclusiveHighConf = NA
mergeddata$TerritorialityPermissiveExclusiveHighConf[which(mergeddata$Territory == "1")] <- "0"
mergeddata$TerritorialityPermissiveExclusiveHighConf[which(mergeddata$Territory == "3")] <- "1"
mergeddata$TerritorialityPermissiveExclusiveHighConf[which(mergeddata$PermissiveExclusive == "Permissive" & mergeddata$Confidence %in% c("Medium", "High"))] <- "0"
mergeddata$TerritorialityPermissiveExclusiveHighConf[which(mergeddata$PermissiveExclusive == "Exclusive" & mergeddata$Confidence %in% c("Medium", "High"))] <- "1"

mergeddata$TerritorialityPermissiveColonialCoopVsExclusive = NA
mergeddata$TerritorialityPermissiveColonialCoopVsExclusive[which(mergeddata$Territory == "1")] <- "0"
mergeddata$TerritorialityPermissiveColonialCoopVsExclusive[which(mergeddata$Territory == "3")] <- "1"
mergeddata$TerritorialityPermissiveColonialCoopVsExclusive[which(mergeddata$colonial_Griesser2023 == "colonial" | mergeddata$HighConfidence_Coop == 1)] <- "0"
mergeddata$TerritorialityPermissiveColonialCoopVsExclusive[which(mergeddata$PermissiveExclusive == "Permissive")] <- "0"
mergeddata$TerritorialityPermissiveColonialCoopVsExclusive[which(mergeddata$PermissiveExclusive == "Exclusive")] <- "1"

# write.csv(mergeddata, paste0("Data_R_Oscine_withTobias_AVONET_WeakStrong", Sys.Date(), ".csv"), row.names = F)


#### Add data for bias tests ----

dfOs = read.csv("Data_R_Oscine_withTobias_AVONET_WeakStrong2025-06-06.csv")
dfOs$HaveFSData = !is.na(dfOs$FemaleSong_Agg01)
dfOs$HaveCBData = !is.na(dfOs$HighConfidence_Coop)

## Holarctic vs Tropical
# Realm (from Jetz and Rubenstein 2011), Region (from Cockburn 2006) - some regions have been studied more intensively. If certain regions are understudied, species would be less likely to have a female song and/or cooperative breeding classification there
# Realm: AT - Afrotropics, PA - Palearctic, NA - Nearctic, NT - Neotropics, IM - Indomalay, AA - Australasia, OC - Oceania. 
# Region: Africa   Antarctic   Australia   Holarctic Indomalayan    Nearctic Neotropical  Palearctic     Unknown  Widespread. 
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
#CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
colnames(JetzData) <- ifelse(colnames(JetzData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(JetzData), "Jetz2011", sep = "_"))
#colnames(CockburnData) <- ifelse(colnames(CockburnData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(CockburnData), "Cockburn2006", sep = "_"))

dfOs2 = merge(dfOs, JetzData[,c("BirdtreeSpecies", "Realm_Jetz2011")], by.x = "species", by.y = "BirdtreeSpecies", all.x = T)
#dfOs2 = merge(dfOs2, CockburnData[,c("BirdtreeSpecies", "Region_Cockburn2006")], by.x = "species", by.y = "BirdtreeSpecies", all.x = T)

#dfOs2 %>% group_by(Realm_Jetz2011, Region_Cockburn2006) %>% count %>% print(n=50)
#dfOs2 %>% filter(is.na(Realm_Jetz2011)) %>% group_by(HighConfidence_Coop, FemaleSong_Agg01,Realm_Jetz2011, Region_Cockburn2006) %>% count %>% print(n=50)

dfOs2$GeographicRegion_Jetz = NA
dfOs2$GeographicRegion_Jetz[which(dfOs2$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
dfOs2$GeographicRegion_Jetz[which(dfOs2$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
#dfOs2$GeographicRegion_Cockburn = NA 
#dfOs2$GeographicRegion_Cockburn[which(dfOs2$Region_Cockburn2006 %in% c("Nearctic", "Holarctic", "Palearctic"))] <- "Holarctic"
#dfOs2$GeographicRegion_Cockburn[which(dfOs2$Region_Cockburn2006 %in% c("Indomalayan", "Australia",   "Neotropical", "Africa"))] <- "Tropical"

table(dfOs2$GeographicRegion_Jetz)

## Dale dichromatism
DaleData = read.csv("plumage_scores_Dale et al 2016.csv")
colnames(DaleData) <- ifelse(colnames(DaleData) == "BirdtreeSpecies", "BirdtreeSpecies", paste(colnames(DaleData), "Dale2015", sep = "_"))
dfOs2 = merge(dfOs2, DaleData[,c("TipLabel_Dale2015", "Female_plumage_score_Dale2015", "Male_plumage_score_Dale2015")], by.x = "species", by.y = "TipLabel_Dale2015", all.x = T)
dfOs2$MaleFemalePlumageDiff = dfOs2$Male_plumage_score_Dale2015-dfOs2$Female_plumage_score_Dale2015
dfOs2$MaleFemalePlumageDiffAbs = abs(dfOs2$MaleFemalePlumageDiff) # absolute value of difference Male-Female
dfOs2$logMaleFemalePlumageDiffAbs = log(abs(dfOs2$MaleFemalePlumageDiff)) # log absolute value of difference Male-Female


## Calculate AVONET dimorphism
avonet_data <- read.csv("AVONET Supplementary dataset 1_Raw Data.csv", 
                        stringsAsFactors = FALSE) # data per individual measured
avonet_data$SpeciesUnderscored = gsub(pattern = " ", "_", avonet_data$Species3_BirdTree)

# Define the measurement columns
measurement_cols <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", 
                      "Beak.Depth", "Tarsus.Length", "Wing.Length", 
                      "Kipps.Distance", "Secondary1", "Hand.wing.Index", 
                      "Tail.Length")

# Filter for records with valid sex data (M or F)
avonet_filtered <- avonet_data %>%
  filter(Sex %in% c("M", "F"))
# Convert measurement columns to numeric
avonet_filtered[measurement_cols] <- lapply(avonet_filtered[measurement_cols], 
                                            function(x) as.numeric(as.character(x)))
# log-transform
avonet_filtered[measurement_cols] <- lapply(avonet_filtered[measurement_cols], 
                                            function(x) log(x))

# 1. Calculate mean and median values by SpeciesUnderscored and Sex
species_sex_summary <- avonet_filtered %>%
  group_by(SpeciesUnderscored, Sex) %>%
  summarise(
    n_individuals = n(),
    across(all_of(measurement_cols), 
           list(mean = ~mean(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  ) %>%
  ungroup()

sum(dfOs2$species[complete.cases(dfOs2[,c("HighConfidence_Coop", "FemaleSong_Agg01")])] %in% species_sex_summary$SpeciesUnderscored)

# 2. Reshape data to have male and female measurements side by side
dimorphism_data <- species_sex_summary %>%
  filter(!is.na(SpeciesUnderscored)) %>%
  pivot_wider(
    id_cols = SpeciesUnderscored,
    names_from = Sex,
    values_from = c(n_individuals, contains("_mean"))
  )

dimorphism_data$PercentLogWingDimorphism_AVONET = (dimorphism_data$Wing.Length_mean_M - dimorphism_data$Wing.Length_mean_F) / dimorphism_data$Wing.Length_mean_F * 100

hist(dimorphism_data$PercentLogWingDimorphism_AVONET)

dfOs2 = merge(dfOs2, dimorphism_data, by.x = "species", by.y = "SpeciesUnderscored", all.x = T)

write.csv(dfOs2, paste0("Data_R_", Sys.Date(), ".csv"), row.names = F)

