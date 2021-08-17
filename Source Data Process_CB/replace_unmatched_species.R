# Replace misspelled species in each source / standardize species names across sources
# Coded by Kate T Snyder
# Started 8/11/2021
#  Updated 8/12/2021 - updated misspelled species reference doc; ran all segments to give output files
#   Next: Check individual source files for duplicated species (or write code to reduce to just one entry per species)
# Last Update: 8/17/2021 - added post-check to ensure all duplicates successfully removed for those sources that had 2+ entries for some species

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

