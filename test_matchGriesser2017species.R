# Test ability to alter common names in our database to match the species names in Griesser et al 2017 PLOS Bio
# Kate T Snyder
# 10/13/2021
# 10/27/2021 - add matching names from https://www.birdwatching.com/software/birdlists/n_c_amer98.html

library(readxl)
Griesser2017 <- read_excel("/Users/kate/Documents/Creanza Lab/Paper PDFs/Griesser et al 2017_Family living sets the stage for cooperative breeding and ecological resilience in birds_PLOSbio data.xlsx")
NCSAmer <- read.csv("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Cooperative Breeding/bird_species_list_NorthCentralSouthAmerica_1998.csv")
eBird <- read_excel("/Users/kate/Documents/Creanza Lab/Comparative Evolution/Cooperative Breeding/eBird-Clements-v2021-integrated-checklist-August-2021_kts.xlsx")

ourdata <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Aggregate_Source_Data_unifiedCommonNames.csv")
ourdata <- NCSAmer # 10/27/21
ourdata <- eBird

library(stringr)

setwd("~/Desktop/CooperativeBreedingEvolution")

ourCommonNames <- ourdata$CommonName_aggregate[which(!is.na(ourdata$CommonName_aggregate))]
ourCommonNames <- NCSAmer$Common   # 10/27/21 
ourCommonNames <- eBird$CommonName



reorderedNamesdf <- set.seed(10)
for ( i in 1:length(ourCommonNames)) {
  tempname <- ourCommonNames[i]    #$CommonNames[i]
  tempnameLower <- str_to_lower(tempname)
  tempvec <- str_split(tempnameLower, " ")[[1]]
 # tempvec <- ournamesSplit[[i]]
  tempvec_ExceptLast <- tempvec[1:length(tempvec)-1]
  tempvec_ExceptLast_collapsed <- str_c(tempvec_ExceptLast, collapse = "")
  tempvec_Last <- tempvec[length(tempvec)]
  ReformattedName <- str_c(tempvec_Last, tempvec_ExceptLast_collapsed)
  ReformattedName_nopunct <- str_remove_all(ReformattedName, "[^a-zA-Z0-9_]")  
  temprow <- c(tempname, ReformattedName_nopunct)
  reorderedNamesdf <- rbind(reorderedNamesdf, temprow)
}
reorderedNamesdf <- as.data.frame(reorderedNamesdf)

# sum(reorderedNamesdf$V1 %in% Griesser2017$species)
sum(Griesser2017$species %in% reorderedNamesdf$V2)
# reorderedNames[which(!reorderedNames %in% Griesser2017$species)]
# 
# str_subset(ournames, "s$")  # note: some species aren't matching because the common name we have for them is plural - mostly tits, also a swallow, warbler, shrike, jay...


reorderedNamesMatched <- reorderedNamesdf[which(reorderedNamesdf$V2 %in% Griesser2017$species),]

dfout <- set.seed(10)
for (k in 1:length(reorderedNamesMatched$V1)) {
  temprow <- reorderedNamesMatched[k,]

  tempspecies <- temprow$V1
  tempGriesser <- temprow$V2
  
  #SciName <- ourdata$BirdtreeSpecies[which(ourdata$CommonName_aggregate %in% tempspecies)]
  SciName <- ourdata$Scientific[which(ourdata$CommonName %in% tempspecies)] # used 10/27
  
  temprowout <- c(tempGriesser, tempspecies, SciName)
  
  if (length(temprowout) != 3) {
    #print(temprowout)
  }
  
  dfout <- rbind(dfout, temprowout)
}
dfout <- as.data.frame(dfout)
colnames(dfout)[1:3] <- c("V1", "CommonName_eBird","Scientific_eBird")

#totalout <- merge(y = Griesser2017, x = dfout, by.y = "species", by.x = "V1", all.y = TRUE, all.x = FALSE)
#write.csv(totalout, file = "Griesser2017_data_matchedCommonNames.csv")

prevDF <- read.csv("Griesser2017_data_matchedCommonNames.csv") # 10/27/21 adding to this file
totalout <- merge(y = prevDF, x = dfout, by.y = "V1", by.x = "V1", all.y = TRUE, all.x = FALSE)
totalout2 <- totalout[,which(colnames(totalout) != "X")]

#only keep ones that match BirdTree species
totalout2$Scientific_eBird <- str_replace(totalout2$Scientific_eBird, pattern = " ", replacement = "_")
ourdata <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Aggregate_Source_Data_unifiedCommonNames.csv")  #for checking match
sum(totalout2$Scientific_eBird %in% ourdata$BirdtreeSpecies)

#totalout2$Scientific_birdlist[which(!totalout2$Scientific_birdlist %in% ourdata$BirdtreeSpecies)] <- NA
totalout3 <- cbind(totalout2, eBirdMatchbirdtree)
write.csv(totalout3, file = "Griesser2017_data_matchedCommonNames_eBird.csv") 
