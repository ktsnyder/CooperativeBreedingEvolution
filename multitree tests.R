# Multitree analyses
# Kate Snyder
# 4/26/2024
# Edited 9/13/2025 - streamlined, made cross-platform functionality

require(phytools)

source("subsettreedata.R")

multitree = read.tree("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/Birdsong - Life History Evolution/BirdzillaHackett4_Stage2_1000trees.tre")
OscineTree = read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
newdata = "Data_R.csv"
dfIn = read.csv(newdata)


#### Brownie ----
#otherlabel = "_Hackett4Oscine_fulltreeQ_UpdatedSongData_"
otherlabel = "_Hackett4Oscine_consensustreeQ_"
# nTreesToSample = 200
# nSimsPerTree = 20
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

# get global Q rates from consensus tree
source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = OscineTree)
qrates= Qout$qrates

source(file.path("Brownie_functions", "browniefunction.R"))
trait1vec = c("HighConfidence_Coop")
trait2vec = c("Song.rep.final")

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

source(file.path("Brownie_functions", "plotbrownie.R"))
plotbrownie(data = allbrownie, columns = c("HighConfidence_Coop", "Syllable.rep.final"), discreteCategoryLabels = c("Non-cooperative", "Cooperative"), islog = T, newpdf = T, otherlabel = "_MultitreeHackett4Oscine_fulltreeQ")

#### Simmap overlap ----
source(file.path("Simmap_Overlap_functions", "CharacterSimmaps_modified.R"))
otherlabel = "_multitree_setQtree"
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

alldfout = set.seed(10)
alldummy = set.seed(10)
for (i in 1:(nTreesToSample)) {
  temptree = multitrees[[i]]
  print(paste(Sys.time(), "tree #", i))
  dfout4 <- CharacterSimmaps_modified(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = FALSE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL, setQratesTree = OscineTree, save_simmaps = F, return_simmaps = F)
  dfDummy4 <- CharacterSimmaps_modified(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = TRUE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL, dummyMethod = "makeSimmap", setQratesTree = OscineTree, save_simmaps = F, return_simmaps = F)
  colnames(dfout4)[which(colnames(dfout4) == "treenum")] <- "simmapNum"
  dfout4$treenum = i
  colnames(dfDummy4)[which(colnames(dfDummy4) == "treenum")] <- "simmapNum"
  dfDummy4$treenum = i
  
  alldfout = rbind(alldfout, dfout4)
  alldummy = rbind(alldummy, dfDummy4)
  
  # if (i %in% (0+c(2,25,50,75,100, 125, 150, 175, 200, 225, 250, 275))) {
  #   write.csv(alldfout, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"), row.names = FALSE)
  #   write.csv(alldummy, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All DUMMY.csv"), row.names = FALSE)
  #   print(paste("Wrote csvs at tree number", i, "--- Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"))
  # }
}
write.csv(alldfout, file.path("Simmap_Overlap_Outputs", paste0(trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv")), row.names = FALSE)
write.csv(alldummy, file.path("Simmap_Overlap_Outputs", paste0(trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All DUMMY.csv")), row.names = FALSE)
print(Sys.time())
source(file.path("Simmap_Overlap_functions","calcHuel_corrected.R"))
Huelout = calcHuel_corrected(alldfout, alldummy, otherlabel = otherlabel)


pdf(file = file.path("Simmap_Overlap_Outputs", paste("Simmap Overlap Outputs/simmap overlap states", trait1, trait2, otherlabel, nTreesToSample, "trees", nSimsPerTree, "simsPerTree", ".pdf")))
Huelout$plots$p3
dev.off()



