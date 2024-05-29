# Multitree analyses
# Kate Snyder
# 4/26/2024

require(phytools)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution")
source("subsettreedata.R")

multitree = read.tree("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/Birdsong - Life History Evolution/BirdzillaHackett4_Stage2_1000trees.tre")
OscineTree = read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
#newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
dfIn = read.csv(newdata)


#### Brownie ----
otherlabel = "_Hackett4Oscine_fulltreeQ_UpdatedSongData_"
otherlabel = "_Hackett4Oscine_consensustreeQ_"
nTreesToSample = 300
nSimsPerTree = 20
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

# get global Q rates from consensus tree
source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = OscineTree)
qrates= Qout$qrates

source("browniefunction.R")
trait1vec = c("HighConfidence_Coop", "HighConfidence_Coop")
trait2vec = c("Song.rep.final", "Syllable.rep.final")

for (f in 1:length(trait1vec)) {
  trait1 = trait1vec[f]
  trait2 = trait2vec[f]
#  subsetout = subsettreedata(columns = c(trait1, trait2), newdata = df, newtree = multitrees)
#  subsetmultitree = subsetout$subsettree
#  subsetdf = subsetout$subsetdf
  
  allbrownie = set.seed(10)
  for (i in 1:nTreesToSample) {
    print(paste("tree number", i))
#    temptree = subsetmultitree[[i]]
    temptree = multitrees[[i]]
    brownieout = browniefunction(columns = c(trait1, trait2), newdata = df, newtree = temptree, nsim = nSimsPerTree, islog = trait2, plotsimmaps = FALSE, otherlabel = otherlabel, phylanovaP = i, setQrates = qrates)
    
    allbrownie = rbind(allbrownie, brownieout)
  }
  colnames(allbrownie)[which(colnames(allbrownie) == "phylanovaP")] <- "TreeNum"
  write.csv(allbrownie, file = paste0(Sys.Date(), " brownie multitree ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree ", otherlabel, " ", trait1, " ", trait2, ".csv"), row.names = F)
  
} # end brownie cycle through traits

source("plotbrownie.R")
plotbrownie(data = allbrownie, columns = c("HighConfidence_Coop", "Syllable.rep.final"), discreteCategoryLabels = c("Non-cooperative", "Cooperative"), islog = T, newpdf = T, otherlabel = "_MultitreeHackett4Oscine_fulltreeQ_UpdatedSongData_")
allbrownie = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/2024-05-16 brownie multitree 300trees 20simsPerTree _Hackett4Oscine_fulltreeQ_UpdatedSongData_ HighConfidence_Coop Song.rep.final.csv")
plotbrownie(data = allbrownie, columns = c("HighConfidence_Coop", "Song.rep.final"), discreteCategoryLabels = c("Non-cooperative", "Cooperative"), islog = T, newpdf = T, otherlabel = "_MultitreeHackett4Oscine_fulltreeQ_UpdatedSongData_")

#### Simmap overlap ----
source("test_trait_overlap_simmaps.R")
otherlabel = "_Hackett4Oscine_multitree_setQtree"
nTreesToSample = 300
nSimsPerTree = 20
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

# Fake multitree with many consensus trees
#multitrees= as.list(rep(consTree, 50))
#otherlabel = "_multiConsensustree"
#nTreesToSample = 10
#nSimsPerTree = 3
# end Fake multitree with many consensus trees

alldfout = set.seed(10)
alldummy = set.seed(10)
for (i in 1:(0+nTreesToSample)) {
  temptree = multitrees[[i]]
  print(paste(Sys.time(), "tree #", i))
  dfout4 <- CharacterSimmaps(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = FALSE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL, setQratesTree = OscineTree)
  dfDummy4 <- CharacterSimmaps(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = TRUE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL, dummyMethod = "makeSimmap", setQratesTree = OscineTree)
  colnames(dfout4)[which(colnames(dfout4) == "treenum")] <- "simmapNum"
  dfout4$treenum = i
  colnames(dfDummy4)[which(colnames(dfDummy4) == "treenum")] <- "simmapNum"
  dfDummy4$treenum = i
  
  alldfout = rbind(alldfout, dfout4)
  alldummy = rbind(alldummy, dfDummy4)
  
  if (i %in% (0+c(2,25,50,75,100, 125, 150, 175, 200, 225, 250, 275))) {
    write.csv(alldfout, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"), row.names = FALSE)
    write.csv(alldummy, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All DUMMY.csv"), row.names = FALSE)
    print(paste("Wrote csvs at tree number", i, "--- Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"))

  }
}
write.csv(alldfout, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"), row.names = FALSE)
write.csv(alldummy, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All DUMMY.csv"), row.names = FALSE)
print(Sys.time())
Huelout = calcHuel(alldfout, alldummy)
Huelout$p2

alldfout = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/HighConfidence_Coop FemaleSong_Agg01 300trees 20simsPerTree_Hackett4Oscine_multitree_setQtree_All REAL.csv")
alldummy = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/HighConfidence_Coop FemaleSong_Agg01 300trees 20simsPerTree_Hackett4Oscine_multitree_setQtree_All DUMMY.csv")
sum(alldummy$ObsProp0Absent < median(alldfout$ObsProp0Absent))/length(alldummy$simmapNum)
sum(alldummy$ObsProp1Absent < median(alldfout$ObsProp1Absent))/length(alldummy$simmapNum)
sum(alldummy$ObsProp0Present < median(alldfout$ObsProp0Present))/length(alldummy$simmapNum)
sum(alldummy$ObsProp1Present < median(alldfout$ObsProp1Present))/length(alldummy$simmapNum)
Huelout = calcHuel(alldfout, alldummy, otherlabel = "_Hackett4Oscine_multitree_300trees-20simsPerTree")
pdf(file = paste("simmap overlap states", trait1, trait2, otherlabel, ".pdf"))
Huelout$p3
dev.off()

require(gridExtra)
plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", " ", otherlabel, " withTransCounts.pdf")
Huelout2 = Huelout[c("p1", "p2", "p3", "p4","p5","p6")]
nPlots = 6
m3 <- marrangeGrob(Huelout2, ncol = 1, nrow = nPlots)
ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")

#### Visualize different trees ----
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
df = read.csv(newdata)
subsetdf = df[which(!is.na(df$FemaleSong_Agg01)),]

tree1 = multitree[[1]]
tree2 = multitree[[2]]

subtree1 = drop.tip(tree1, tip = which(!tree1$tip.label %in% subsetdf$species))
subtree2 = drop.tip(tree2, tip = which(!tree2$tip.label %in% subsetdf$species))
comparePhylo(subtree1, subtree2, plot=T, cex = 0.1)

