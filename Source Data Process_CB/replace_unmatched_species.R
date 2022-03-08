# Replace misspelled species in each source / standardize species names across sources
# Coded by Kate T Snyder
# Started 8/11/2021
#  Updated 8/12/2021 - updated misspelled species reference doc; ran all segments to give output files
#   Next: Check individual source files for duplicated species (or write code to reduce to just one entry per species)
# Last Update: 8/17/2021 - added post-check to ensure all duplicates successfully removed for those sources that had 2+ entries for some species
# Updated: 3/7/2022 - add Odom et al 2014

setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")

library(phytools)
birdtree <- read.nexus("birdzillatreeMaybeConsensus.nex")
thousandtrees <- read.tree("~/Desktop/CooperativeBreedingEvolution/BirdzillaHackett3.tre")
birdtree <- thousandtrees[[1]]
AllBirdtreeSpecies <- birdtree$tip.label

misspelledbirds <- read.csv("inconsistent_species_names_BirdTree.csv") 
duplicatemisspellings <- misspelledbirds$in_database[duplicated(misspelledbirds$in_database)]

### Biagolini
BiagoliniDF <- read.csv("Biagolini_SuppTable1.csv")
BiagoliniDF$Species[which(!BiagoliniDF$Species %in% AllBirdtreeSpecies)]
BirdtreeSpecies <- set.seed(10)
df <- BiagoliniDF
for (i in 1:length(df[,1])) {
  tempspecies <- df$Species[i]
  if (tempspecies %in% AllBirdtreeSpecies) {
    BirdtreeSpecies[i] <- tempspecies
  } else if (tempspecies %in% misspelledbirds$in_database) {
    speciesreplacement <- misspelledbirds$species_in_birdtree[which(misspelledbirds$in_database %in% tempspecies)]
    BirdtreeSpecies[i] <- speciesreplacement
    print(paste("", tempspecies, speciesreplacement, sep = ", "))
  } else {
    BirdtreeSpecies[i] <- NA
    print(paste("", tempspecies, "NOT IN MISSING TABLE", sep = ", "))
  }
}
dfWithBTCol <- cbind(BirdtreeSpecies, df)
write.csv(dfWithBTCol, file = "Biagolini_data_BirdTreeNames.csv")

### Cockburn
CockburnDF <- read.csv("Cockburn2006_data.csv")
cockburnWeirdSpecies<- CockburnDF$SpeciesScientific[which(!CockburnDF$SpeciesScientific %in% AllBirdtreeSpecies)]
# cockburnmisspelled <- cockburnWeirdSpecies %in% misspelledbirds$in_database
# missingCockburnSpecies <- cockburnWeirdSpecies[!cockburnmisspelled]
# missingCockburnSpeciesDF <- CockburnDF[which(CockburnDF$SpeciesScientific %in% missingCockburnSpecies),]
# write.csv(missingCockburnSpeciesDF, file = "misspelled_missing_Cockburn_species.csv")
BirdtreeSpecies <- set.seed(10)
df <- CockburnDF
for (i in 1:length(df[,1])) {
  tempspecies <- df$SpeciesScientific[i]
  if (tempspecies %in% AllBirdtreeSpecies) {
    BirdtreeSpecies[i] <- tempspecies
  } else if (tempspecies %in% misspelledbirds$in_database) {
    speciesreplacement <- misspelledbirds$species_in_birdtree[which(misspelledbirds$in_database %in% tempspecies)]
    BirdtreeSpecies[i] <- speciesreplacement
  } else {
    BirdtreeSpecies[i] <- NA
  }
}
dfWithBTCol <- cbind(BirdtreeSpecies, df)
write.csv(dfWithBTCol, file = "Cockburn2006_data_BirdTreeNames.csv")


### Downing
DowningDF <- read.csv("Downing supp table2_kts edited.csv")
DowningDF$species[which(DowningDF$species %in% AllBirdtreeSpecies)]
#All good


### Dunn
DunnDF <- read.csv("Dunn supp data 977sp.csv")
DunnDF$SciName_BL[which(!DunnDF$SciName_BL %in% AllBirdtreeSpecies)]
#All good

### Jetz
JetzDF <- read.csv("Jetz2011_data_clean.csv")
JetzDF$Species[which(!JetzDF$Species %in% AllBirdtreeSpecies)]
BirdtreeSpecies <- set.seed(10)
df <- JetzDF
for (i in 1:length(df[,1])) {
  tempspecies <- df$Species[i]
  if (tempspecies %in% AllBirdtreeSpecies) {
    BirdtreeSpecies[i] <- tempspecies
  } else if (tempspecies %in% misspelledbirds$in_database) {
    speciesreplacement <- misspelledbirds$species_in_birdtree[which(misspelledbirds$in_database %in% tempspecies)]
    BirdtreeSpecies[i] <- speciesreplacement
   # print(paste("", tempspecies, speciesreplacement, sep = ", "))
  } else {
    BirdtreeSpecies[i] <- NA
    print(paste("", tempspecies, "NOT IN MISSING TABLE", sep = ", "))
  }
}
dfWithBTCol <- cbind(BirdtreeSpecies, df)
write.csv(dfWithBTCol, file = "Jetz_data_BirdTreeNames.csv")

### Riehl
RiehlDF <- read.csv("Riehl 2013 supp data columns lines_kts edited.csv")
RiehlDF$SpeciesName[which(!RiehlDF$SpeciesName %in% AllBirdtreeSpecies)]
BirdtreeSpecies <- set.seed(10)
df <- RiehlDF
for (i in 1:length(df[,1])) {
  tempspecies <- df$SpeciesName[i]
  if (tempspecies %in% AllBirdtreeSpecies) {
    BirdtreeSpecies[i] <- tempspecies
  } else if (tempspecies %in% misspelledbirds$in_database) {
    speciesreplacement <- misspelledbirds$species_in_birdtree[which(misspelledbirds$in_database %in% tempspecies)]
    BirdtreeSpecies[i] <- speciesreplacement
    # print(paste("", tempspecies, speciesreplacement, sep = ", "))
  } else {
    BirdtreeSpecies[i] <- NA
    print(paste("", tempspecies, "NOT IN MISSING TABLE", sep = ", "))
  }
}
dfWithBTCol <- cbind(BirdtreeSpecies, df)
write.csv(dfWithBTCol, file = "Riehl_data_BirdTreeNames.csv")

### Rubenstein Lovette
RubLovDF <- read.csv("RubensteinLovette2007_data_clean.csv")
RubLovDF$Species[which(!RubLovDF$Species %in% AllBirdtreeSpecies)]
BirdtreeSpecies <- set.seed(10)
df <- RubLovDF
for (i in 1:length(df[,1])) {
  tempspecies <- df$Species[i]
  if (tempspecies %in% AllBirdtreeSpecies) {
    BirdtreeSpecies[i] <- tempspecies
  } else if (tempspecies %in% misspelledbirds$in_database) {
    speciesreplacement <- misspelledbirds$species_in_birdtree[which(misspelledbirds$in_database %in% tempspecies)]
    BirdtreeSpecies[i] <- speciesreplacement
    # print(paste("", tempspecies, speciesreplacement, sep = ", "))
  } else {
    BirdtreeSpecies[i] <- NA
    print(paste("", tempspecies, "NOT IN MISSING TABLE", sep = ", "))
  }
}
dfWithBTCol <- cbind(BirdtreeSpecies, df)
write.csv(dfWithBTCol, file = "RubensteinLovette_data_BirdTreeNames.csv")

### Griesser
GriesserDF <- read.csv("Griesser supp table1.csv")
GriesserDF$scientific.name[which(!GriesserDF$scientific.name %in% AllBirdtreeSpecies)]
#All good
GriesserDF$scientific.name[duplicated(GriesserDF$scientific.name)]

# add new column to sources with species misspelled/missing from BirdTree species

#dfWithBTCol <- cbind(BirdtreeSpecies, df)
#write.csv(dfWithBTCol, file = "Cockburn2006_data_BirdTreeNames.csv")
#speciesNotInBirdtree <- df$SpeciesScientific[which(is.na(BirdtreeSpecies))]
#


df <- read.csv("RubensteinLovette_data_BirdTreeNames.csv")
duplicatemisspellings <- df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]

df <- read.csv("Riehl_data_BirdTreeNames.csv")
duplicatemisspellings <- df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]

df <- read.csv("Jetz_data_BirdTreeNames.csv")
duplicatemisspellings <- df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
unique(duplicatemisspellings)
# a bunch of dups
df <- read.csv("Jetz_data_BirdTreeNames_nodups.csv")
df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
sum(df$BirdtreeSpecies %in% AllBirdtreeSpecies)
length(df$BirdtreeSpecies)


df <- read.csv("Cockburn2006_data_BirdTreeNames.csv")
duplicatemisspellings <- df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
# several dups
df <- read.csv("Cockburn2006_data_BirdTreeNames_nodups.csv")
df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
sum(df$BirdtreeSpecies %in% AllBirdtreeSpecies)
length(df$BirdtreeSpecies)


df <- read.csv("Biagolini_data_BirdTreeNames.csv")
duplicatemisspellings <- df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
# just Delichon_urbicum and a couple NAs
df <- read.csv("Biagolini_data_BirdTreeNames_nodups.csv")
df$BirdtreeSpecies[duplicated(df$BirdtreeSpecies)]
sum(df$BirdtreeSpecies %in% AllBirdtreeSpecies)
length(df$BirdtreeSpecies)



# Odom et al 2014 - female song
df <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/FemaleSongData_OdomEtal2014_PresentAbsentSubset.csv")
df$Latin_binomial[which(df$Latin_binomial == "Chlorophoneus_sulfureopectus")] <- "Telophorus_sulfureopectus"
df$Latin_binomial[which(df$Latin_binomial == "Phylidonyris_nigra")] <- "Phylidonyris_niger"
df$Latin_binomial[which(df$Latin_binomial == "Phylidonyris_pyrrhoptera")] <- "Phylidonyris_pyrrhopterus"
df$Latin_binomial[which(df$Latin_binomial == "Sugomel_niger")] <- "Certhionyx_niger"
df$Latin_binomial[which(df$Latin_binomial == "Carterornis_chrysomela" )] <- "Monarcha_chrysomela"
df$Latin_binomial[which(df$Latin_binomial == "Diphyllodes_respublica")] <- "Cicinnurus_respublica"
df$Latin_binomial[which(df$Latin_binomial == "Drepanornis_albertisi")] <- "Epimachus albertisi"
df$Latin_binomial[which(df$Latin_binomial == "Drepanornis_bruijnii")] <- "Epimachus bruijnii"
df$Latin_binomial[which(df$Latin_binomial ==  "Manucodia_atra")] <- "Manucodia_ater"
df$Latin_binomial[which(df$Latin_binomial == "Phonygammus_keraudrenii")] <- "Manucodia_keraudrenii"
df$Latin_binomial[which(df$Latin_binomial == "Seleucidis_melanoleuca")] <- "Seleucidis_melanoleucus"
df$Latin_binomial[which(df$Latin_binomial == "Amblyornis_inornatus")] <- "Amblyornis_inornata"
df$Latin_binomial[which(df$Latin_binomial == "Vireo_atricapillus")] <- "Vireo_atricapilla"
df$Latin_binomial[which(df$Latin_binomial == "Calamanthus_pyrrhopygius")] <- "Hylacola_pyrrhopygia"
df$Latin_binomial[which(df$Latin_binomial == "Pyrrholaemus_sagittatus")] <- "Chthonicola_sagittatus"
df$Latin_binomial[which(df$Latin_binomial == "Callaeas_cinerea")] <- "Callaeas_cinereus"
df$Latin_binomial[which(df$Latin_binomial == "Climacteris_melanura")] <- "Climacteris_melanurus"
df$Latin_binomial[which(df$Latin_binomial == "Climacteris_rufa")] <- "Climacteris_rufus"
df$Latin_binomial[which(df$Latin_binomial == "Coloeus_monedula")] <- "Corvus_monedula"
df$Latin_binomial[which(df$Latin_binomial == "Finschia_novaeseelandiae")] <- "Mohoua_novaeseelandiae"
df$Latin_binomial[which(df$Latin_binomial == "Bocagia_minuta")] <- "Tchagra_minutus"
df$Latin_binomial[which(df$Latin_binomial == "Chlorophoneus_bocagei")] <- "Telophorus_bocagei"
df$Latin_binomial[which(df$Latin_binomial == "Chlorophoneus_nigrifrons")] <- "Telophorus_nigrifrons"
df$Latin_binomial[which(df$Latin_binomial == "Tchagra_senegala")] <- "Tchagra_senegalus"
df$Latin_binomial[which(df$Latin_binomial == "Cissomela_pectoralis")] <- "Certhionyx_pectoralis"
df$Latin_binomial[which(df$Latin_binomial == "Foulehaio_carunculata")] <- "Foulehaio_carunculatus"
#df$Latin_binomial[which(df$Latin_binomial == "Philemon_plumigenis")] <- "Philemon_moluccensis" - already have this species in this data
df$Latin_binomial[which(df$Latin_binomial == "Diphyllodes_magnificus")] <- "Cicinnurus_magnificus"
df$Latin_binomial[which(df$Latin_binomial == "Petroica_boodang")] <- "Petroica_multicolor"
#df$Latin_binomial[which(df$Latin_binomial == "Rhipidura_albiscapa")] <- "Rhipidura_fuliginosa" - already have this species in this data

df$Latin_binomial[duplicated(df$Latin_binomial)]
write.csv(df, file = "FemaleSongData_OdomEtal2014_PresentAbsentSubset_BirdTreeNames.csv", row.names = FALSE)
