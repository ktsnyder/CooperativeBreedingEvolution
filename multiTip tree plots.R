# # Double-tip ACE trees
# 

source("subsettreedata.R")
#newdata = "2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
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
offsetDenom = 1
tipsize = 0.1

AllColsToPlot = AllColsToPlot[1:2]

columnSubsetTree = c("HighConfidence_Coop","FemaleSong_Agg01")
subsetout <- subsettreedata(columns = columnSubsetTree, newdata = newdata, newtree = treefile)
subsetdf = subsetout$subsetdf
subsetdf$AnyNoncoopEqualsNoncoop = subsetdf$MeanCoopTie2Noncoop
subsetdf$AnyNoncoopEqualsNoncoop[which(subsetdf$SourceDiscrepancy == 1)] = 0
subsettree = subsetout$subsettree
py = c("purple", "orange", "white")
py2 = c("blue","red", "white")
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
