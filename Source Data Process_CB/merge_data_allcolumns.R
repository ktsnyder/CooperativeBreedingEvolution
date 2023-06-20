## Merge individual source datasets
## Created by Aleyna Loughran-Pierce
## Created 10/10/2020
## Last modified: 12/09/2021 by Kate Snyder - include all columns from each source file
## 3/7/2022 - add in Griesser 2017 data, female song (Odom 2014) data
## 3/8/2022 - add in BOW data
## 9/7/2022 - female song dataset
## Last edited 5/30/2023


setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/")

ourdatabase <- read.csv("SongData_R_Update.csv", stringsAsFactors = FALSE)
ourdatabaserefs <- read.csv("SupplementDataRefs_Update.csv")
# colnames(ourdatabaserefs)[1] <- "BirdtreeSpecies"  # changed in csv
BiagoliniData <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
CockburnData <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv", stringsAsFactors = FALSE)
DowningData <- read.csv("Downing supp table2_kts edited.csv", stringsAsFactors = FALSE)
DunnData <- read.csv("Dunn supp data 977sp.csv", stringsAsFactors = FALSE)
GriesserData <- read.csv("Griesser supp table1.csv", stringsAsFactors = FALSE)
JetzData <- read.csv("Jetz_data_BirdTreeNames_nodups.csv", check.names = FALSE)
#RiehlData <- read.csv("Riehl 2013 supp data columns lines_kts edited.csv", stringsAsFactors = FALSE)
RiehlData <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Riehl_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
RubensteinLovetteData <- read.csv("RubensteinLovette_data_BirdTreeNames.csv", stringsAsFactors = FALSE)
#OdomFSData <- read.csv("FemaleSongData_OdomEtal2014_PresentAbsentSubset_BirdTreeNames.csv")
Griesser2017Data <- read.csv("Griesser2017 matched species names.csv")
Griesser2017Data$BirdtreeFormat[which(Griesser2017Data$Griesser2017_speciesnames == "finchyellowthighed")] <- "Pselliophorus_tibialis"
BOWData <- read.csv("BOW Cooperative Breeding Data.csv")
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

# Merge Remes data
RemesData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data.csv")
RemesClimateData = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Remes et al 2015 PNAS doi_10.5061_dryad.02jk0__v1/PNAS_data_climatic.csv")
sum(RemesData$Species_name %in% newdf14$BirdtreeSpecies)
RemesDataAll = merge(RemesData, RemesClimateData)
colnames(RemesDataAll)[which(colnames(RemesDataAll) == "Species_name")] <- "Species_Name_Remes"
newdf17 = merge(newdf16, RemesDataAll, by.x = "BirdtreeSpecies", by.y = "Species_Name_Remes", all = T, suffixes = c("","_Remes"))

# remove non-birdtrees
sum(!newdf17$BirdtreeSpecies %in% birdtree$tip.label)
newdf18 = newdf17[which(newdf17$BirdtreeSpecies %in% birdtree$tip.label),]

# remove duplicated rows
newdf19 = newdf18[which(!duplicated(newdf18$BirdtreeSpecies)),]

# write file - then this file used in CoopBreed_species_summary.R
write.csv(newdf19, file = paste0(Sys.Date(),"_Aggregate_Source_Data_AllColumns.csv"), row.names = FALSE)


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
# newdf8 <- merge(DatabaseUpdate2019_Clements, newdf7, by = "BirdtreeSpecies", all = TRUE, suffixes = c("_ClementsGuide", ""))
# newdf9 <- merge(newdf8, DatabaseUpdate2019_LisSpo, by = "BirdtreeSpecies", all = TRUE)
# newdf10 <- merge(Clements2021subset, newdf9, by.x = "Scientific", by.y = "BirdtreeSpecies", all.x = FALSE, all.y = TRUE, suffixes = c("ClementsGuide2021",))
#write.csv(newdf10, file = "Aggregate_Source_Data_AllColumns_Incl2019UpdateSources_new.csv", row.names = FALSE)
  #edited file in excel to remove some unnecessary columns; copying to Google Sheet for Teo 12/9/2021


#### Female song ----
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

FSdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-05-30_Female Song Data_GSheetDownload.csv")
colnames(FSdata) %in% colnames(Aggdf)
colnames(Aggdf) %in% colnames(FSdata)
sum(colnames(Aggdf) %in% colnames(cbdf))
AggdfNewFS = Aggdf[which(!is.na(Aggdf$Mikula1_Duetting_vs_NoDuetting)),c("BirdtreeSpecies", "Mikula1_Duetting_vs_NoDuetting", "Mikula2_Duetting_vs_NoFS", "Mikula3_FemaleSoloSong_vs_Duetting", "Mikula4_FemaleSoloSong_vs_NoFS")]
!AggdfNewFS$BirdtreeSpecies %in% FSdata$BirdtreeSpecies

FSdata$BirdtreeSpecies[which(!FSdata$BirdtreeSpecies %in% birdtree$tip.label)]
#FSdata = FSdata[,c("BirdtreeSpecies","Odom_FemaleSong","Webb_FemaleSong","FemaleSong_Aggregated")]
FSdata$Odom_FemaleSong[which(FSdata$Odom_FemaleSong == "")] <- NA
FSdata$Webb_FemaleSong[which(FSdata$Webb_FemaleSong == "")] <- NA
cbdf$species[which(!cbdf$species %in% birdtree$tip.label)]

orderedNewFS = merge(FSdata, AggdfNewFS, all = T)
#write.csv(orderedNewFS, file = "2023-06-16_Female Song Data_GSheet2023-05-30_addedMikula.csv", row.names = F)  #pasted this new file back into the Female Song Data Google Sheet https://docs.google.com/spreadsheets/d/1kEijryIyLgJ-Co40pVKkxqAPUDEoWGkCd6YY4tVTMVM/edit#gid=0
sum(!orderedNewFS$BirdtreeSpecies %in% birdtree$tip.label)
orderedNewFS %>% group_by(FemaleSong_Aggregated, Mikula1_Duetting_vs_NoDuetting, Mikula2_Duetting_vs_NoFS, Mikula3_FemaleSoloSong_vs_Duetting, Mikula4_FemaleSoloSong_vs_NoFS) %>% summarize(n=n())
newFS = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-16_Female Song Data_GSheetDownload.csv")
newFS %>% group_by(Mikula_SummaryClass, FemaleSong_Aggregated) %>% summarize(n=n())
ReviewFS = newFS[which(newFS$Mikula_SummaryClass == "Absent Female Song" & newFS$FemaleSong_Aggregated == "Present" | (newFS$Mikula_SummaryClass == "Duetting" & newFS$FemaleSong_Aggregated == "Absent") | newFS$Mikula_SummaryClass == "Female Solo Song" & newFS$FemaleSong_Aggregated == "Absent"),]
AddFS = newFS[which(newFS$Mikula_SummaryClass == "Absent Female Song" & is.na(newFS$FemaleSong_Aggregated) | newFS$Mikula_SummaryClass == "Female Solo Song" & is.na(newFS$FemaleSong_Aggregated)),]

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
