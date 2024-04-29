# Multitree analyses
# Kate Snyder
# 4/26/2024

require(phytools)

source("subsettreedata.R")

multitree = read.tree("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/Birdsong - Life History Evolution/BirdzillaHackett4_Stage2_1000trees.tre")
OscineTree = read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"


#### Brownie ----
otherlabel = "_Hackett4Oscine_fulltreeQ_"
nTreesToSample = 400
nSimsPerTree = 20
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

source("browniefunction.R")
trait1vec = c("HighConfidence_Coop", "HighConfidence_Coop", "HighConfidence_Coop")
trait2vec = c("Song.rep.final", "Syllable.rep.final", "Syll.song.final")

for (f in 1:length(trait1vec)) {
  trait1 = trait1vec[f]
  trait2 = trait2vec[f]
#  subsetout = subsettreedata(columns = c(trait1, trait2), newdata = df, newtree = multitrees)
#  subsetmultitree = subsetout$subsettree
#  subsetdf = subsetout$subsetdf
  
  allbrownie = set.seed(10)
  for (i in 200:nTreesToSample) {
    print(paste("tree number", i))
#    temptree = subsetmultitree[[i]]
    temptree = multitrees[[i]]
    brownieout = browniefunction(columns = c(trait1, trait2), newdata = df, newtree = temptree, nsim = nSimsPerTree, islog = trait2, plotsimmaps = FALSE, otherlabel = otherlabel, phylanovaP = i)
    
    allbrownie = rbind(allbrownie, brownieout)
  }
  colnames(allbrownie)[which(colnames(allbrownie) == "phylanovaP")] <- "TreeNum"
  write.csv(allbrownie, file = paste0(Sys.Date(), " brownie multitree ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree ", otherlabel, " ", trait1, " ", trait2, ".csv"), row.names = F)
  
} # end brownie cycle through traits


#### Simmap overlap ----
source("test_trait_overlap_simmaps.R")
otherlabel = "_Hackett4Oscine_multitree"
nTreesToSample = 200
nSimsPerTree = 50
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
df = read.csv(newdata)
df = df[which(df$species %in% OscineTree$tip.label),]

multitrees = multitree[1:nTreesToSample]

alldfout = set.seed(10)
alldummy = set.seed(10)
for (i in 1:nTreesToSample) {
temptree = multitrees[[i]]

dfout4 <- CharacterSimmaps(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = FALSE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c(trait1,trait2), df = df, tree =  temptree, dummy = TRUE, nsims = nSimsPerTree, treelabel = otherlabel, datalabel = NULL, dummyMethod = "makeSimmap")
colnames(dfout4)[which(colnames(dfout4) == "treenum")] <- "simmapNum"
dfout4$treenum = i
colnames(dfDummy4)[which(colnames(dfDummy4) == "treenum")] <- "simmapNum"
dfDummy4$treenum = i

alldfout = rbind(alldfout, dfout4)
alldummy = rbind(alldummy, dfDummy4)
}
write.csv(alldfout, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All REAL.csv"), row.names = FALSE)
write.csv(alldummy, paste0("Simmap Overlap Outputs/", trait1, " ", trait2, " ", nTreesToSample, "trees ", nSimsPerTree, "simsPerTree", otherlabel, "_All DUMMY.csv"), row.names = FALSE)
print(Sys.time())
Huelout = calcHuel(alldfout, alldfdummy)
