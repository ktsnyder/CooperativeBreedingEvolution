## Plot phylogeny with multiple trait tips labels

treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

dfIn = read.csv("Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04.csv")
dfIn = read.csv("Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-04-29.csv")
dfIn = read.csv("Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-08.csv")
passertree = read.nexus('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex')
passertree <- multi2di(passertree)
passertree$node.label = NULL
oscinetree = read.nexus(treefile)

source("subsettreedata.R")
library(phylosignal)

dfIn$JiayingDuet01 = NA
dfIn$JiayingDuet01[dfIn$JiayingDuet == "absent"] <- 0
dfIn$JiayingDuet01[dfIn$JiayingDuet == "present"] <- 1
dfIn$JiayingSolo01 = NA
dfIn$JiayingSolo01[dfIn$JiayingSolo == "absent"] <- 0
dfIn$JiayingSolo01[dfIn$JiayingSolo == "present"] <- 1

# Filter out NAs and prune tree
complete_cases <- !is.na(dfIn$JiayingSolo01)
dfIn_J = dfIn[complete_cases,]
pruned_tree <- drop.tip(passertree, setdiff(passertree$tip.label, dfIn_J$species))

trait_data = dfIn_J$NF_Solo_Duet
names(trait_data) <- dfIn_J$species
binary_Solo <- dfIn_J$JiayingSolo01
names(binary_Solo) <- dfIn_J$species
binary_Duet <- dfIn_J$JiayingDuet01
names(binary_Duet) <- dfIn_J$species




binarySolo_df <- data.frame(
  species = names(binary_Solo),
  binary_trait = binary_Solo
)
binaryDuet_df <- data.frame(
  species = names(binary_Duet),
  binary_trait = binary_Duet
)

phylo_d_result <- phylo.d(
  data = binaryDuet_df, 
  phy = pruned_tree,
  names.col = species,  # This needs to be the column name with species names
  binvar = binary_trait  # This needs to be the column name with the trait
)

library(phytools)

# Convert to numeric vector named by species 
trait_vec <- as.numeric(dfIn_J$NF_Solo_Duet)+1
names(trait_vec) <- dfIn_J$species

# Calculate Pagel's lambda for the ordinal trait
lambda_resultPlus1 <- phylosig(tree_J, trait_vec, method="lambda", test=TRUE)
print(lambda_result)

# Alternatively, calculate Blomberg's K
k_result <- phylosig(pruned_tree, pruned_trait, method="K", test=TRUE)
print(k_result)

library(caper)
phylo_d_result <- phylo.d(data=dfIn_J, phy=tree_J, names.col=species, binvar=binary_Solo)
print(phylo_d_result)

subsetdf = dfIn[which(!is.na(dfIn$FemaleSong_Agg01) | !is.na(dfIn$NF_Solo_Duet)),]
sum(!df$species %in% oscinetree$tip.label)
sum(!df$species %in% passertree$tip.label)


subsettree <- ape::drop.tip(passertree, setdiff(passertree$tip.label, df$species))



#### make for loop ----
subsetCBFS = subsettreedata(columns = c("HighConfidence_Coop", "FemaleSong_Agg01"), newdata = dfIn, newtree = treefile)
subsetdf = subsetCBFS$subsetdf
subsettree = subsetCBFS$subsettree

ToPlotColorSet1 = c("HighConfidence_Coop")
ToPlotColorSet2 = "FemaleSong_Agg01"
ToPlotColorSet3 = c("Territory")
ToPlotColorSet4 = c("TerritorialityWeakVsStrong", "TerritorialityWeakVsStrongHighConf", "TerritorialityPermissiveExclusive", "TerritorialityPermissiveExclusiveHighConf") 
ToPlotColorSet3 = c("NF_Solo_Duet")
ToPlotColorSet4 = c("Duet") # Tobias Duet
AllColsToPlot = c(ToPlotColorSet1, ToPlotColorSet2, ToPlotColorSet3, ToPlotColorSet4)

# duet stuff
py1 = c("blue","red", "white")
py2 = c("purple", "darkgreen", "white")
py3 = c("purple","orange", "green", "white")
py4 = c("purple","green", "white")

# territory stuff
py1 = c("blue","red", "white") # coop breed
py2 = c("purple", "darkgreen", "white") # female song
py3 = c("brown","orange", "hotpink", "white") # 3-state territory
py4 = c("brown","hotpink", "white") # 2-state territory

# Using the adjustcolor() function from base R
adjust_transparency <- function(color_vector, alpha = 0.5, white_alpha = 0.5) {
  sapply(color_vector, function(color) {
    if (color == "white") {
      adjustcolor(color, alpha.f = white_alpha)
    } else {
      adjustcolor(color, alpha.f = alpha)
    }
  })
}

# Apply to your color vectors
py1 <- adjust_transparency(py1)
py2 <- adjust_transparency(py2)
py3 <- adjust_transparency(py3)
py4 <- adjust_transparency(py4)

offsetDenom = 2
tipsize = 0.1

#pdf("/Users/kate/phylo HighConfidence_Coop FemaleSong_Agg01 Duet transparent tips_wSpecies.pdf", width = 12, height = 13)
pdf("phylo HighConfidence_Coop FemaleSong_Agg01 TerritoryWeakStrongPermissive transparent tips_wSpeciesAligned.pdf", width = 15, height = 16)
#par(mar=c(4,4,1,1))
plot.phylo(subsettree, type = "f", show.tip.label = T, align.tip.label = T, cex = 0.1)

checktips = subsettree$tip.label

for (i in 1:length(AllColsToPlot)) {
  tempCol = AllColsToPlot[i]
  
  print(unique(subsetdf[,tempCol]))
  print(class(subsetdf[,tempCol]))
  
  if (0 %in% subsetdf[,tempCol]) {
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 1)]
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 2)])] <- 2
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, tempCol]))])] = max(subsetdf[, tempCol], na.rm = T) + 1
    print("0 present in")
  } else {
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] != 1)]
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 2)])] <- 1
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 3)])] <- 2
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, tempCol]))])] = max(subsetdf[, tempCol], na.rm = T) + 1
    print("0 not present in")
  }
  
  print(tempCol)
  print(table(treetiplabels))
  
  if (tempCol %in% ToPlotColorSet1) {
    colorsToUse = py1
  } else if (tempCol %in% ToPlotColorSet2) {
    colorsToUse = py2
  } else if (tempCol %in% ToPlotColorSet3) {
    colorsToUse = py3
  } else if (tempCol %in% ToPlotColorSet4) {
    colorsToUse = py4
  }
  
  colorindex = as.numeric(treetiplabels)+1
  colorvector = colorsToUse[as.numeric(treetiplabels)+1]
  checktips = cbind(checktips, treetiplabels, colorindex, colorvector)
  
  print(length(treetiplabels))
  tiplabels(pch=21,bg=colorvector, col = colorvector, cex=tipsize, offset = i/offsetDenom)
}

checktips = as.data.frame(checktips)
colnames(checktips)[2:4] <- paste0(colnames(checktips)[2:4],".Coop")
colnames(checktips)[5:7] <- paste0(colnames(checktips)[5:7],".FS")
colnames(checktips)[8:10] <- paste0(colnames(checktips)[8:10],".Terr")
colnames(checktips)[11:13] <- paste0(colnames(checktips)[11:13],".TerrWS")
colnames(checktips)[14:16] <- paste0(colnames(checktips)[14:16],".TerrWSHC")

checktips %>% group_by(treetiplabels.Coop, colorindex.Coop, colorvector.Coop) %>% count
checktips %>% group_by(treetiplabels.Terr, colorindex.Terr, colorvector.Terr) %>% count

# Combine color definitions for a legend
allpy = c(py1[1:2], py2[1:2], py3[1:3], py4)# , py2nonsig[3])

#alllabs = c("Noncooperative", "Cooperative", "Female Song Absent", "Female Song Present", "Solo FS", "Duet")#

# Define legend labels for both sets of traits
alllabs = c("Noncooperative", "Cooperative", "Female Song Absent", "Female Song Present", "Terr1", "Terr2", "Terr3", "Weak/PermissiveTerr", "Strong/ExclusiveTerr")#, "Simmap Overlap Not Significant")

# Add a legend to the plot
legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")
# Add legend of all tip rows 
legend("topleft", legend = AllColsToPlot, cex = 0.9, bty="n", title = "Tip states - innermost to outermost")
dev.off()




# using dotTree
?dotTree
subsetdf$HighConfidence_Coop[which(subsetdf$HighConfidence_Coop == 0)] <- "Non-coop"
subsetdf$HighConfidence_Coop[which(subsetdf$HighConfidence_Coop == 1)] <- "Coop"
subsetdf$FemaleSong_Agg01[which(subsetdf$FemaleSong_Agg01 == 0)] <- "Female Song Absent"
subsetdf$FemaleSong_Agg01[which(subsetdf$FemaleSong_Agg01 == 1)] <- "Female Song Present"
subsetdf$TerritorialityWeakVsStrong[which(subsetdf$TerritorialityWeakVsStrong == 0)] <- "Weak"
subsetdf$TerritorialityWeakVsStrong[which(subsetdf$TerritorialityWeakVsStrong == 1)] <- "Strong"
subsetdf$TerritorialityWeakVsStrongHighConf[which(subsetdf$TerritorialityWeakVsStrongHighConf == 0)] <- "Weak"
subsetdf$TerritorialityWeakVsStrongHighConf[which(subsetdf$TerritorialityWeakVsStrongHighConf == 1)] <- "Strong"
  
subsetMatrix = as.matrix(cbind(subsetdf$HighConfidence_Coop, subsetdf$FemaleSong_Agg01)) #, subsetdf$Territory)) #, subsetdf$TerritorialityWeakVsStrong, subsetdf$TerritorialityWeakVsStrongHighConf))
rownames(subsetMatrix) <- subsetdf$species

colors<-setNames(c("blue","red"),c(0,1))

pdf("test dotplot.pdf", height = 145, width = 25)
dotTree(subsettree, x = subsetMatrix, fsize = 0.1, length = 80, data.type="discrete", cex = 2.7, colors = colors, pch = 19)
dev.off()
