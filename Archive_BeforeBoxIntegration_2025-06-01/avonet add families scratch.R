### scratch get method to extract clades
### Kate Snyder
### 10/3/2023
### 

setwd("~/Desktop/CooperativeBreedingEvolution")

library(readxl)
AvonetBirdTree = read_excel("/Users/kate/Downloads/AVONET Supplementary dataset 1.xlsx", sheet = 4)
AvonetCrosswalk = read_excel("/Users/kate/Downloads/AVONET Supplementary dataset 1.xlsx", sheet = 10)
AvonetBirdLife = read_excel("/Users/kate/Downloads/AVONET Supplementary dataset 1.xlsx", sheet = 2)
AvonetEBird = read_excel("/Users/kate/Downloads/AVONET Supplementary dataset 1.xlsx", sheet = 3)

sum(!AvonetCrosswalk$Species1 %in% AvonetBirdTree$Species3)
sum(!AvonetCrosswalk$Species2 %in% AvonetBirdTree$Species3)
sum(!AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species1)
sum(!AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species2)
sum(!AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species2 & !AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species1)
sum(AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species2 | AvonetBirdTree$Species3 %in% AvonetCrosswalk$Species1)


newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data_R.csv"
df = read.csv(newdata)
AvonetAll = merge(AvonetCrosswalk, AvonetBirdTree[,c("Species3","Family3", "Order3")], by.x = "Species1", by.y = "Species3", all = T)
AvonetAll = merge(AvonetAll, AvonetBirdTree[,c("Species3","Family3", "Order3")], by.x = "Species2", by.y = "Species3", all = T, suffixes = c("_BirdtreeMatchSpecies1", "_BirdtreeMatchSpecies2"))
AvonetAll = merge(AvonetAll, AvonetBirdLife[,c("Species1","Family1", "Order1")], by.x = "Species1", by.y = "Species1", all = T)
AvonetAll = merge(AvonetAll, AvonetEBird[,c("Species2","Family2", "Order2")], by.x = "Species2", by.y = "Species2", all = T)
sum(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) | !is.na(AvonetAll$Family3_BirdtreeMatchSpecies1))
sum(AvonetAll$Species2==AvonetAll$Species1 & is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) & !is.na(AvonetAll$Family3_BirdtreeMatchSpecies1), na.rm=T)
sum(AvonetAll$Species2 != AvonetAll$Species1 & is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) & !is.na(AvonetAll$Family3_BirdtreeMatchSpecies1), na.rm = T)
AvonetAll$BirdtreeSpecies = NA
AvonetAll$BirdtreeSpecies[which(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies1))] <- AvonetAll$Species1[which(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies1))]
AvonetAll$BirdtreeSpecies[which(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies2))] <- AvonetAll$Species2[which(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies2))]
AvonetAll = merge(AvonetAll, )
sum(is.na(AvonetAll$BirdtreeSpecies))
sum(!AvonetBirdTree$Species3 %in% AvonetAll$BirdtreeSpecies)
sum(AvonetAll$Family3_BirdtreeMatchSpecies2 == AvonetAll$Family3_BirdtreeMatchSpecies1, na.rm=T)
sum(AvonetAll$Family3_BirdtreeMatchSpecies2 != AvonetAll$Family3_BirdtreeMatchSpecies1, na.rm=T)
sum(is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) & !is.na(AvonetAll$Family3_BirdtreeMatchSpecies1), na.rm=T)
sum(!is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) & is.na(AvonetAll$Family3_BirdtreeMatchSpecies1), na.rm=T)
sum(is.na(AvonetAll$Family3_BirdtreeMatchSpecies2) & is.na(AvonetAll$Family3_BirdtreeMatchSpecies1), na.rm=T)
AvonetAll[which(AvonetAll$Family3_BirdtreeMatchSpecies2 != AvonetAll$Family3_BirdtreeMatchSpecies1),]

library(stringr)
library(dplyr)
AvonetAll$Species2 = str_replace(AvonetAll$Species2, pattern = " ", "_")
AvonetAll$Species1 = str_replace(AvonetAll$Species1, pattern = " ", "_")
length(df$species)
sum(df$species %in% AvonetAll$Species2)
sum(df$species %in% AvonetAll$Species1)
sum(AvonetAll$Species2 %in% df$species)
sum(AvonetAll$Species1 %in% df$species)

df2 = merge(df, AvonetAll, by.x = "species", by.y = "Species2")
length(unique(df2$species))
which(duplicated(df2$species))
tib1 = df2 %>% group_by(species, Family1, Family2) %>% count
which(duplicated(tib1$species))
which(tib1$n > 1)


df3 = df2[which(!duplicated(df2$species)),]

write.csv(df3, file = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv")
sum(is.na(df3$Family3_BirdtreeMatchSpecies2))
df3[which(is.na(df3$Family3_BirdtreeMatchSpecies2)),]




source("subsettreedata.R")
Hacktree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")
treelabel <- "Hackett"
subsettreedata(columns=c("MeanCoopTie2Noncoop", "FemaleSong_Agg01"), newdata = AvonetAll, newtree = Hacktree)



AvonetAll$Family.w.eBird = NA
AvonetAll$Family.w.BirdLife = NA



library(phytools)
getDescendants(Hacktree, 7092)
getSisters(Hacktree, 5)

subsettree
subsetdf
strahler = extract.strahlerNumber(subsettree, 4, plot=TRUE)
getDescendants(subsettree, 1200)
paintedsubtree = paintSubTree(subsettree,1200, state = 0)
plotSimmap(paintedsubtree, fsize = 0.1)

estDiversity(subsettree, ) #biogeographic area-based

trait1 = subsetdf$MeanCoopTie2Noncoop
trait2 = subsetdf$FemaleSong_Agg01
names(trait1) = names(trait2) = subsetdf$species
fitPagel(tree = subsettree, x = trait1, y = trait2, model = "ARD")  # This function fit's Pagel's (1994) model for the correlated evolution of two binary characters.


get_ancestor_nodes(6, subsettree)
get_common_ancestor_nodes(c(12:100), subsettree)

get_ancestor_nodes <- function(tip, phylo) {
  # Check if the input phylo is of class "phylo"
  if (class(phylo) != "phylo") stop("Input phylo must be of class 'phylo'")
  
  # Determine if the input tip is a label or a number
  if (is.character(tip)) {
    # Check if the tip label exists in the phylogeny
    if (!tip %in% phylo$tip.label) stop("Tip label not found in the phylogeny")
    # Get the node number of the tip label
    tip_node <- which(phylo$tip.label == tip)
  } else if (is.numeric(tip)) {
    # Check if the tip number is valid
    if (tip <= 0 || tip > length(phylo$tip.label)) stop("Invalid tip number")
    tip_node <- tip
  } else {
    stop("Input tip must be a character (tip label) or numeric (tip number)")
  }
  
  # Initialize a vector to store ancestor nodes
  ancestor_nodes <- numeric()
  
  # Loop to find all ancestor nodes
  while (tip_node != 0) {
    parent_node <- which(phylo$edge[, 2] == tip_node)
    if (length(parent_node) == 0) break  # Exit loop if no parent node is found
    ancestor_nodes <- c(ancestor_nodes, phylo$edge[parent_node, 1])
    tip_node <- phylo$edge[parent_node, 1]  # Update tip_node to the found parent node
  }
  
  return(ancestor_nodes)
}

get_common_ancestor_nodes <- function(tips, phylo) {
  # Check if the input phylo is of class "phylo"
  if (class(phylo) != "phylo") stop("Input phylo must be of class 'phylo'")
  
  # Initialize a list to store ancestor nodes for each tip
  all_ancestor_nodes <- list()
  
  # Loop over each tip in the input vector
  for (tip in tips) {
    # Determine if the input tip is a label or a number
    if (is.character(tip)) {
      # Check if the tip label exists in the phylogeny
      if (!tip %in% phylo$tip.label) stop(paste("Tip label", tip, "not found in the phylogeny"))
      # Get the node number of the tip label
      tip_node <- which(phylo$tip.label == tip)
    } else if (is.numeric(tip)) {
      # Check if the tip number is valid
      if (tip <= 0 || tip > length(phylo$tip.label)) stop(paste("Invalid tip number:", tip))
      tip_node <- tip
    } else {
      stop("Input tips must be characters (tip labels) or numerics (tip numbers)")
    }
    
    # Initialize a vector to store ancestor nodes
    ancestor_nodes <- numeric()
    
    # Loop to find all ancestor nodes
    while (tip_node != 0) {
      parent_node <- which(phylo$edge[, 2] == tip_node)
      if (length(parent_node) == 0) break  # Exit loop if no parent node is found
      ancestor_nodes <- c(ancestor_nodes, phylo$edge[parent_node, 1])
      tip_node <- phylo$edge[parent_node, 1]  # Update tip_node to the found parent node
    }
    
    # Add the ancestor nodes of the current tip to the list
    all_ancestor_nodes[[as.character(tip)]] <- ancestor_nodes
  }
  
  # Find and return the common ancestor nodes among all tips
  common_ancestor_nodes <- Reduce(intersect, all_ancestor_nodes)
  return(common_ancestor_nodes)
}
