# # Double-tip ACE trees
# 

source("subsettreedata.R")
#newdata = "2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
newdata = "Data_R_2025-06-09.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

specialspecies = c("Malurus_splendens", "Campylorhynchus_griseus",  "Prionops_plumatus", "Lamprotornis_superbus", "Turdoides_gymnogenys",  "Mohoua_ochrocephala", "Cisticola_robustus", "Agelaioides_badius",  "Myrmecocichla_formicivora", "Philetairus_socius", "Melanodryas_cucullata") #"Climacteris_affinis","Platysteira_peltata","Sericornis_frontalis", "Malurus_cyaneus", "Campylorhynchus_nuchalis",  "Sitta_pusilla", "Manorina_melanocephala",
maluridae Malurus_cyaneus Malurus_splendens
meliphagidae Manorina_melanocephala
timaliidae Turdoides_gymnogenys
troglodytes Campylorhynchus_griseus Campylorhynchus_nuchalis
malaconotidae Prionops_plumatus Prionops_retzii
acanthizidae Sericornis_frontalis Mohoua_ochrocephala
cisticolidae Cisticola_chiniana(maybe, slight discrepancy); Cisticola_anonymus Cisticola_robustus (all non-FS)
icteridae/grackle? - Agelaioides_badius
sturnidae Lamprotornis_superbus
climacteridae Climacteris_affinis Climacteris_melanurus Climacteris_picumnus Climacteris_rufus
muscicapidae Myrmecocichla_formicivora
passeridae Philetairus_socius (non-FS)
petroicidae Melanodryas_cucullata (non-FS)
platysteridae Platysteira_peltata
sittidae Sitta_pusilla







columns = c("FemaleSong_Agg01", "HighConfidence_Coop")
subsetout <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile)
subsetdf = subsetout$subsetdf
subsetdf$AnyNoncoopEqualsNoncoop = subsetdf$MeanCoopTie2Noncoop
subsetdf$AnyNoncoopEqualsNoncoop[which(subsetdf$SourceDiscrepancy == 1)] = 0
subsettree = subsetout$subsettree

# write.nexus(subsettree, file = "subsettree FemaleSong_Agg01 HighConfidence_Coop 2024-04-12di.nex")
subsettree = read.nexus("subsettree FemaleSong_Agg01 HighConfidence_Coop 2024-04-12di.nex")

tipsize = 0.1

# Construct the filename for saving the plot, based on specified columns
#filename = paste(columns[1], columns[2], "FS fan phylo multi tips.pdf")
filename = paste("fan phylo multi tips names SomeMoreNodesFSxCB5 FemaleSong_Agg01 HighConfidence_Coop 2024-04-12subtree.pdf")
#filename = "FemaleSong_Agg01 fan phylo HC-Coop tips_candidateLabelSpecies4.pdf"

# Generate filename base for both PDF and PNG
filename_base <- gsub("\\.pdf$", "", filename)

# Save the plot to a PDF file
pdf(filename, height = 8, width = 12)
par(mar = c(0,0,0,0))
offsetDenom = 1 # 0.5 #3

# Plot the phylogenetic tree without showing tip labels, in a fan layout
plot.phylo(subsettree, type = "f", show.tip.label = TRUE, align.tip.label = TRUE, cex = 0.08, show.node.label = TRUE)

# Define colors for the first set of traits
py = c("purple","orange")

# Identify which tips match the first trait condition and set their colors
treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[1]] == 1)]

# Add the first set of tip labels with custom colors based on the first trait
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=tipsize, offset = 1/offsetDenom)

# Define colors for the second set of traits
#py2 = c("blue","red", "white")
py2 = c("gray", "red")
#py2nonsig = c("blue","red", "gray")

# Identify which tips match the second trait condition and set their colors
#treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[2]]
treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"HighConfidence_Coop"] == 1)]
treetiplabels2[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[,"HighConfidence_Coop"]))])] = 2

# Add the second set of tip labels with custom colors based on the second trait, offset from the first set
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=tipsize, offset = 2/offsetDenom)

whichNodes = sort(c( 1042:1062, 1105, 1118:1120, 1156, 1217, 1238:1241, 1291:1294, 1317, 1330, 1355:1357, 1374:1380, 1447:1448, 1473, 1550, 1495:1497, 1559, 1561:1568, 1570, 1636:1638, 1667, 1725:1726, 1780, 1789:1791, 1792:1794, 1802:1807, 1821:1822, 1846:1847, 1894, 1900:1906, 1916, 1921, 1977, 1989, 1986, 1991:1995, 2034, 2051, 2053:2054, 2066:2068, 2071, 2080,   1353, 1358, 1506, 1476, 1592, 1581, 1917, 1923, 1960, 1981, 2043, 2062, 2072, 1157, 1130, 2076,2081, 1111, 1108, 1158:1160, 1173, 1479, 2036, 1990, 1953, 1954, 1878:1879, 1719:1721, 1610, 2024,2028,2004, 1978, 2023, 1973, 1974, 1924,  1948, 1895, 1848, 1880, 1882, 1866, 1849, 1825, 1842, 1795, 1812, 1783, 1730, 1748, 1756, 1714, 1715, 1707, 1708, 1668, 1660, 1653, 1644, 1729, 1591, 1575, 1554, 1522, 1557, 1498, 1475, 1477, 1470, 1393, 1359, 1408, 1332, 1334, 1346, 1318, 1300, 1296, 1284, 1278, 1286, 1218, 1212, 1211, 1201, 1144, 1125, 1112, 1090, 1516, 1440, 1423, 1383, 1397, 1333, 1277, 1285, 1188, 1175, 1121, 1151, 1081, 1976, 1459, 1460))
#nodelabels(node = whichNodes, pch = NULL, thermo = NULL, cex = 0.4, height = 1.2, width = 1.2) # text node numbers

discretetraitvec1 = subsetdf[,columns[1]]
names(discretetraitvec1) = subsetdf$species
circles=ace(x=discretetraitvec1,phy=subsettree,type="discrete",model="ARD")
#nodelabels(node = whichNodes, thermo=circles$lik.anc, piecol=py, frame = "circle", height = 1.2, width = 1.2)
circlesdf = as.data.frame(circles$lik.anc)
circlesdf = cbind(rownames(circlesdf), circlesdf)
colnames(circlesdf)[1] <- "nodeNum"
nodelabels(node = whichNodes, thermo=circlesdf[which(circlesdf$nodeNum %in% whichNodes),2:3], piecol=py, frame = "circle", height = 0.8, width = 1.8, adj = c(0,0.4), horiz = TRUE)


discretetraitvec2 = subsetdf[,columns[2]]
names(discretetraitvec2) = subsetdf$species
circles2=ace(x=discretetraitvec2,phy=subsettree,type="discrete",model="ARD")
circlesdf2 = as.data.frame(circles2$lik.anc)
circlesdf2 = cbind(rownames(circlesdf2), circlesdf2)
colnames(circlesdf2)[1] <- "nodeNum"
nodelabels(node = whichNodes, thermo=circlesdf2[which(circlesdf2$nodeNum %in% whichNodes),2:3], piecol=py2, frame = "circle", height = 0.8, width = 1.8, adj = c(0,-0.4), horiz = TRUE)
#nodelabels(thermo = circles2$lik.anc, piecol=py2, frame = "circle", height = 0.6, width = 1.8, horiz = TRUE)
#nodelabels(pch = NULL, thermo = NULL, cex = 0.1, frame = "none") # text node numbers


specialFamilies = unique(subsetdf$Family3_BirdtreeMatchSpecies2_AVONET[which(subsetdf$species %in% specialspecies)])
candidateSpecies= subsetdf[which(subsetdf$species %in% specialspecies), c("species", "Family3_BirdtreeMatchSpecies2_AVONET", "HighConfidence_Coop", "FemaleSong_Agg01")]
# add ring of species to label
treetiplabels3 = subsettree$tip.label %in% specialspecies  
# tiplabsdf= as.data.frame(cbind(subsettree$tip.label, treetiplabels3))
# tiplabsdf = merge(tiplabsdf, subsetdf[,c("species", "Family3_BirdtreeMatchSpecies2_AVONET")], by.x = "V1", by.y = "species")
# tiplabsdf$Family3_BirdtreeMatchSpecies2_AVONET[which(tiplabsdf$treetiplabels3 == "FALSE")] = "FALSE"
# treetiplabels3Factor = factor(tiplabsdf$Family3_BirdtreeMatchSpecies2_AVONET, levels = c("FALSE", specialFamilies))
# #names(treetiplabels3Factor) = tiplabsdf$V1
# treetiplabels3Numeric =as.numeric(treetiplabels3Factor)
# View(cbind(tiplabsdf, treetiplabels3Factor, treetiplabels3Numeric))

#pySpecial = c("white", "darkgreen")
#tiplabels(pch=21,bg=pySpecial[as.numeric(treetiplabels3)+1], col = pySpecial[as.numeric(treetiplabels3)+1], cex=tipsize*1.5, offset = 3/offsetDenom)
##pySpecial = c("white", "black", "gray", "brown", rainbow(10))
##tiplabels(pch=21,bg=pySpecial[as.numeric(treetiplabels3Factor)], col = pySpecial[as.numeric(treetiplabels3Factor)], cex=tipsize*2, offset = 3/offsetDenom)


#alllabs = c("FALSE", specialFamilies)

# Add a legend to the plot
#legend("bottomleft", legend = alllabs, cex = 0.9, fill=pySpecial, bty="n")

dev.off()

# Now create PNG output with same content
png(file.path("Outputs/Figures/PNG", paste0(filename_base, ".png")), 
    width = 12*150, height = 8*150, res = 150)
par(mar = c(0,0,0,0))

# Plot the phylogenetic tree without showing tip labels, in a fan layout
plot.phylo(subsettree, type = "f", show.tip.label = TRUE, align.tip.label = TRUE, cex = 0.08, show.node.label = TRUE)

# Define colors for the first set of traits
py = c("purple","orange")

# Identify which tips match the first trait condition and set their colors
treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[1]] == 1)]

# Add the first set of tip labels with custom colors based on the first trait
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=tipsize, offset = 1/offsetDenom)

# Define colors for the second set of traits
py2 = c("blue","red")

# Identify which tips match the second trait condition
treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[2]] == 1)]

# Add the second set of tip labels with custom colors based on the second trait
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=tipsize, offset = 2/offsetDenom)

dev.off()


# Add more cooperative breeding classification dots
treetiplabels3 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"MeanCoopTie2Coop"] == 1)]  
treetiplabels3[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "MeanCoopTie2Coop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels3)+1], col = py2[as.numeric(treetiplabels3)+1], cex=tipsize, offset = 3/offsetDenom)

treetiplabels4 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"CornwallisCoop"] == 1)]
treetiplabels4[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "CornwallisCoop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels4)+1], col = py2[as.numeric(treetiplabels4)+1], cex=tipsize, offset = 4/offsetDenom)

treetiplabels5 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"CockburnCoop"] == 1)]
treetiplabels5[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "CockburnCoop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels5)+1], col = py2[as.numeric(treetiplabels5)+1], cex=tipsize, offset = 5/offsetDenom)

treetiplabels6 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"BiagoliniCoop"] == 1)]
treetiplabels6[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "BiagoliniCoop"]))])] = 2
tiplabels(pch=21,bg=py2nonsig[as.numeric(treetiplabels6)+1], col = py2nonsig[as.numeric(treetiplabels6)+1], cex=tipsize, offset = 6/offsetDenom)

treetiplabels7 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"DowningCoop"] == 1)]
treetiplabels7[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "DowningCoop"]))])] = 2
tiplabels(pch=21,bg=py2nonsig[as.numeric(treetiplabels7)+1], col = py2nonsig[as.numeric(treetiplabels7)+1], cex=tipsize, offset = 7/offsetDenom)

treetiplabels8 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"Griesser2017Coop"] == 1)]
treetiplabels8[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "Griesser2017Coop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels8)+1], col = py2[as.numeric(treetiplabels8)+1], cex=tipsize, offset = 8/offsetDenom)

treetiplabels9 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"System_Jetz2011"] == "Cooperative")]
treetiplabels9[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "System_Jetz2011"]))])] = 2
tiplabels(pch=21,bg=py2nonsig[as.numeric(treetiplabels9)+1], col = py2nonsig[as.numeric(treetiplabels9)+1], cex=tipsize, offset = 9/offsetDenom)

treetiplabels10 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"AnyNoncoopEqualsNoncoop"] == 1)]
treetiplabels10[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "AnyNoncoopEqualsNoncoop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels10)+1], col = py2[as.numeric(treetiplabels10)+1], cex=tipsize, offset = 10/offsetDenom)

treetiplabels11 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,"AnyCoopEqualsCoop"] == 1)]
treetiplabels11[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, "AnyCoopEqualsCoop"]))])] = 2
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels11)+1], col = py2[as.numeric(treetiplabels11)+1], cex=tipsize, offset = 11/offsetDenom)


# Combine color definitions for a legend
allpy = c(py2[1:2], py, py2nonsig[3])

# Define legend labels for both sets of traits
alllabs = c("Noncooperative", "Cooperative", "Female Song Absent", "Female Song Present", "Simmap Overlap Not Significant")

# Add a legend to the plot
legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")

# Uncomment the next line to finish writing to the PDF file and close it
dev.off()



#### make for loop ----
ToPlotColorSet1 = "FemaleSong_Agg01"
ToPlotColorSet2 = c("HighConfidence_Coop", "MeanCoopTie2Coop", "MeanCoopTie2Noncoop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "DaleCoop", "AnyNoncoopEqualsNoncoop", "AnyCoopEqualsCoop")
ToPlotColorSetNonsig = c("JetzCoop", "DowningCoop", "BiagoliniCoop")

AllColsToPlot = c("FemaleSong_Agg01", "HighConfidence_Coop", "MeanCoopTie2Coop", "MeanCoopTie2Noncoop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "DaleCoop", "JetzCoop", "DowningCoop", "BiagoliniCoop", "AnyNoncoopEqualsNoncoop", "AnyCoopEqualsCoop")

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

AllColsToPlot = AllColsToPlot[1:2]

offsetDenom = 1
tipsize = 0.1


columnSubsetTree = c("HighConfidence_Coop","FemaleSong_Agg01", "logMass_AVONET")
subsetout <- subsettreedata(columns = columnSubsetTree, newdata = jackbest, newtree = treefile)
subsetdf = subsetout$subsetdf
#subsetdf$AnyNoncoopEqualsNoncoop = subsetdf$MeanCoopTie2Noncoop
#subsetdf$AnyNoncoopEqualsNoncoop[which(subsetdf$SourceDiscrepancy == 1)] = 0
subsettree = subsetout$subsettree
py = c("purple", "orange", "white")
py2 = c("blue","red", "white")
py3 = c("green","black", "white")
py4 = c()
py2nonsig = c("blue","red", "gray")


# Generate filename base for both PDF and PNG
filename_base2 <- "phylo_HighConfidence_Coop_FemaleSong_Agg01_tips_FSxCB_tree"

pdf(paste0(filename_base2, ".pdf"), width = 9, height = 8)
plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)

for (i in 1:length(AllColsToPlot)) {
  tempCol = AllColsToPlot[i]
  treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 1)]
  treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, tempCol]))])] = 2
  
  if (tempCol %in% ToPlotColorSet1) {
    colorsToUse = py
  } else if (tempCol %in% ToPlotColorSet2) {
    colorsToUse = py2
  } else if (tempCol %in% ToPlotColorSetNonsig) {
    colorsToUse = py2nonsig
  } else if (tempCol %in% ToPlotColorSet3) {
    colorsToUse = py3
  }
  
  tiplabels(pch=21,bg=colorsToUse[as.numeric(treetiplabels)+1], col = colorsToUse[as.numeric(treetiplabels)+1], cex=tipsize, offset = i/offsetDenom)
  
}
# Combine color definitions for a legend
allpy = c(py2[1:2], py[1:2])# , py2nonsig[3])

# Define legend labels for both sets of traits
alllabs = c("Noncooperative", "Cooperative", "Female Song Absent", "Female Song Present")#, "Simmap Overlap Not Significant")

# Add a legend to the plot
legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")
dev.off()

# Add legend of all tip rows 
legend("topleft", legend = AllColsToPlot, cex = 0.9, bty="n", title = "Tip states - innermost to outermost")


#### jackknife phylopath models ----

# Load required libraries
library(ape)
library(dplyr)

newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/Outputs/PhylopathJackknife/JackknifedSpeciesData_detailed_models_JackknifeSpecies_n875_2025-07-08.csv"
jackknife_data = read.csv(newdata)
jackknife_data <- jackknife_data %>%
  group_by(RemovedSpecies) %>%
  mutate(
    FullDatasetBestModelInTop = as.numeric(any(model == "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS" & delta_CICc < 2)),
    DownsampledDataBestModelInTop = as.numeric(any(model == "Y2_TERR→FS_FS→COOP_MASS→FS_MASS→COOP" & delta_CICc < 2))
  ) %>%
  ungroup()
jackbest = jackknife_data[which(jackknife_data$delta_CICc == 0),]
jackbest= jackbest[!duplicated(jackbest$RemovedSpecies),]
jackbest$species = jackbest$RemovedSpecies
class(jackbest)
jackbest = as.data.frame(jackbest)

columnSubsetTree = c("HighConfidence_Coop","FemaleSong_Agg01", "logMass_AVONET")
subsetout <- subsettreedata(columns = columnSubsetTree, newdata = jackbest, newtree = treefile)
subsetdf = subsetout$subsetdf
subsettree = subsetout$subsettree

Ntip <- length(subsettree$tip.label)

# Define columns to plot
ToPlotColorSet1 = "FemaleSong_Agg01"
ToPlotColorSet2 = "HighConfidence_Coop"
ToPlotColorSet3 = "TerritorialityWeakVsStrong"
ToPlotColorSet4 = "model_labels"  # Using simplified model labels
ToPlotColorSet5 = "FullDatasetBestModelInTop"
ToPlotColorSet6 = "DownsampledDataBestModelInTop"
AllColsToPlot = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong", "FullDatasetBestModelInTop", "DownsampledDataBestModelInTop", "model_labels")

subsetdf$model_labels <- gsub("→", "to", subsetdf$model)
subsetdf$model_labels <- gsub("_TERRtoFS", "", subsetdf$model_labels)
subsetdf$model_labels <- gsub("_MASStoFS", "", subsetdf$model_labels)
subsetdf$model_labels <- gsub("X1_", "", subsetdf$model_labels)
subsetdf$model_labels <- gsub("Y2_", "", subsetdf$model_labels)
subsetdf$model_labels <- gsub("COOP", "CB", subsetdf$model_labels)


# Define color schemes for each column
py = c("purple", "orange", "white")  # FemaleSong_Agg01
py2 = c("blue", "red", "white")      # HighConfidence_Coop
py3 = c("green", "black", "white")   # TerritorialityWeakVsStrong

# Define colors for the two main models (model column)
py4 = c("darkcyan", "sienna", "white")  # Model 1 CBtoFS_TERRtoCB, Model 2 FStoCB_MASStoCB, NA
py5 = c("gray", "darkcyan") # Model 1 in top models <2 CICc?
py6 = c("gray", "sienna") # Model 2 in top models <2 CICc?
# Model 1: CBtoFS_TERRtoCB (darkviolet) - corresponds to original X1 model
# Model 2: FStoCB_MASStoCB (goldenrod) - corresponds to original Y2 model

# Define the two main models that will be plotted as dots
main_models = c("CBtoFS_TERRtoCB", "FStoCB_MASStoCB")

# Define tip size and offset parameters
tipsize = 0.2
offsetDenom = 2  # Adjust this to control spacing between columns

# Create PDF
filename_base2 <- "phylo_HighConfidence_Coop_FemaleSong_Agg01_Territory_Model_TopModels_tips"
pdf(paste0(filename_base2, ".pdf"), width = 12, height = 11)
par(mar = c(6, 4, 4, 4), xpd = TRUE)  # Increase margins and allow plotting outside

# Plot phylogeny
plot.phylo(subsettree, type = "f", show.tip.label = TRUE, align.tip.label = TRUE, cex = 0.1)

# Get the last tip positions for calculating angles
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xx <- lastPP$xx
yy <- lastPP$yy

# Track which tips need text labels for rare models
tips_with_text_models = list()

# Plot each column
for (i in 1:length(AllColsToPlot)) {
  tempCol = AllColsToPlot[i]
  
  if (tempCol == "model_labels") {
    # Special handling for model column
    # Create numeric vector for model states (like binary traits)
    # 0 = CBtoFS_TERRtoCB, 1 = FStoCB_MASStoCB, 2 = rare models (no dot), 3 = NA
    model_states = rep(3, length(subsettree$tip.label))  # Default to NA
    
    for (j in 1:length(subsettree$tip.label)) {
      tip_species = subsettree$tip.label[j]
      if (tip_species %in% subsetdf$species) {
        idx = which(subsetdf$species == tip_species)
        if (length(idx) > 0) {
          model_val = subsetdf$model_labels[idx]
          if (!is.na(model_val)) {
            if (model_val == main_models[1]) {
              model_states[j] = 0  # CBtoFS_TERRtoCB model
            } else if (model_val == main_models[2]) {
              model_states[j] = 1  # FStoCB_MASStoCB model
            } else {
              model_states[j] = 2  # Rare models (will be text)
              # Store tip index and model name for text labeling
              tips_with_text_models[[length(tips_with_text_models) + 1]] = 
                list(index = j, model = model_val, offset = i/offsetDenom)
            }
          }
        }
      }
    }
    
    # Create color vector based on states
    # Use transparent color for rare models so no dot appears
    py4_extended = c(py4[1:2], "transparent", py4[3])  # darkviolet, goldenrod, transparent, white
    
    # Plot all dots at once using the same method as binary traits
    tiplabels(pch = 21, bg = py4_extended[model_states + 1], 
              col = py4_extended[model_states + 1], 
              cex = tipsize, offset = i/offsetDenom)
    
  } else {
    # Standard binary trait handling
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[, tempCol] == 1)]
    treetiplabels[which(subsettree$tip.label %in% subsetdf$species[which(is.na(subsetdf[, tempCol]))])] = 2
    
    if (tempCol %in% ToPlotColorSet1) {
      colorsToUse = py
    } else if (tempCol %in% ToPlotColorSet2) {
      colorsToUse = py2
    } else if (tempCol %in% ToPlotColorSet3) {
      colorsToUse = py3
    } else if (tempCol %in% ToPlotColorSet5) {
      colorsToUse = py5
    } else if (tempCol %in% ToPlotColorSet6) {
      colorsToUse = py6
    }
    
    tiplabels(pch = 21, bg = colorsToUse[as.numeric(treetiplabels) + 1], 
              col = colorsToUse[as.numeric(treetiplabels) + 1], 
              cex = tipsize, offset = i/offsetDenom)
  }
}

# Add text labels for rare models with proper angles
if (length(tips_with_text_models) > 0) {
  # Get the maximum x coordinate (where aligned tips are)
  max_x <- max(xx[1:Ntip])
  
  for (item in tips_with_text_models) {
    # Get x, y coordinate for this tip
    tip_x_original <- xx[item$index]
    tip_y <- yy[item$index]
    
    # Calculate angle from origin
    angle_rad <- atan2(tip_y, tip_x_original)  # Angle from horizontal
    angle_deg <- angle_rad * 180 / pi
    
    
    # Adjust text angle so it reads outward
    text_angle <- angle_deg
    if (angle_deg > 90 || angle_deg < -90) {
      text_angle <- angle_deg + 180
      adj_val <- 1  # Right-align text on left side
    } else {
      adj_val <- 0  # Left-align text on right side
    }
    
    # Calculate radial offset distance
    # item$offset is i/offsetDenom where i=4 for the 4th column
    # So offset_distance is how far beyond max_x to place the text
    offset_distance = item$offset
    extra_offset = 1  # Additional offset for text beyond dots
    
    # Calculate radius (distance from center)
    radius = max_x + offset_distance + extra_offset
    text_x <- radius * cos(angle_rad)
    text_y <- radius * sin(angle_rad)
    
    # Add text label at the calculated position
    text(text_x, text_y,
         labels = item$model,
         srt = text_angle,
         adj = adj_val,
         cex = 0.5,
         xpd = TRUE)
  }
}

# Create combined legend
# Binary traits
allpy_binary = c(py2[1:2], py[1:2], py3[1:2])
alllabs_binary = c("Noncooperative", "Cooperative", 
                   "Female Song Absent", "Female Song Present",
                   "Weak/No Territory", "Strong Territory")

# Model dots
allpy_models = py4[1:2]
alllabs_models = c("CBtoFS_TERRtoCB (best in full dataset, n=324)", 
                   "FStoCB_MASStoCB (best in downsamples, n=507)")

# Add legends
legend("bottomleft", legend = alllabs_binary, cex = 0.8, 
       fill = allpy_binary, bty = "n", title = "Binary Traits")

legend("bottomright", legend = alllabs_models, cex = 0.8, 
       fill = allpy_models, bty = "n", title = "Main Models (dots - also contain TERRtoFS and MASStoFS)")

# Add note about text labels
legend("topright", 
       legend = c("CBtoFS_TERRtoCB_MASStoCB (n=1)",
                  "CBtoFS_TERRtoCB_MASStoTERR (n=1)",
                  "FStoCB (n=29)",
                  "FStoCB_MASStoTERR (n=4)",
                  "FStoCB_TERRtoCB (n=9)"),
       cex = 0.7, bty = "n", 
       title = "Other models shown as text\n(all models contain TERRtoFS and MASStoFS)")

# Add column order legend
legend("topleft", legend = paste(1:6, c("HighConfidence_Coop", "FemaleSong_Agg01", 
                                        "TerritorialityWeakVsStrong", "BestFullDatasetModel_InTopModels", "BestModelFromDownsamples_InTopModels", "Model"), 
                                 sep = ". "), 
       cex = 0.8, bty = "n", title = "Column order (inner to outer)")

dev.off()

# Print summary of rare models for reference
cat("\nRare models that appear as text labels:\n")
rare_models = unique(subsetdf$model_labels[!subsetdf$model_labels %in% main_models])
for (rm in rare_models) {
  count = sum(subsetdf$model_labels == rm, na.rm = TRUE)
  cat(sprintf("%s: %d species\n", rm, count))
}

