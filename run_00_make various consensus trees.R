# Make various consensus trees
# Split and updated from script "make various consensus trees.R"

require(phytools)

if (!exists("passermultitree") & file.exists("PasserineMultiphy1000Hackett4_nondicho.nex")) {
  passermultitree = read.nexus("PasserineMultiphy1000Hackett4_nondicho.nex")  
  print("successfully read file PasserineMultiphy1000Hackett4_nondicho.nex")
}

if (!exists("nTrees")) {
  nTrees = 10
  print("performing consensus phylogeny computation for 10 trees for the sake of example. To use another number of trees (at least 2, up to 1000), define nTrees = #.")
}

## script used to subset trees to just Passerine species - commented out but provided for posterity
#thousandtrees = read.tree("BirdzillaHackett4_Stage2_1000trees.tre")
#onepassertree <- read.nexus("2021-08-31ConsensusPasserineTreeEricson10_1000.nex")
#datanames <- onepassertree$tip.label

# matezillamultitree <- list()
# for (i in 1:1000) {
#   dropformatezilla <- set.seed(10)
#   temptree <- thousandtrees[[i]]
#   havedatavec <- temptree$tip.label %in% as.character(datanames) 
#   dropformatezilla <- which(havedatavec == FALSE)
#   matezillamultitree[[i]] <- drop.tip(phy = temptree,tip=dropformatezilla)
# }
# print(paste("done trimming trees", Sys.time()))

multitree <- passermultitree[1:nTrees]

# consensus using mean.edge, use default absent edges = "zero" in edge means
startconsense <- Sys.time()
startconsense
matezillaconsensus <- consensus.edges(multitree, method = "mean.edge", if.absent = "zero")
endconsense <- Sys.time()
endconsense
endconsense-startconsense
write.nexus(matezillaconsensus,file=paste(Sys.Date(),"ConsensusPasserineTreeHackett4_", nTrees,"_mean-edge_absent-zero.nex",sep=""))

# consensus using mean.edge, ignore absent edges in edge means
startconsense <- Sys.time()
startconsense
matezillaconsensus <- consensus.edges(multitree, method = "mean.edge", if.absent = "ignore")
endconsense <- Sys.time()
endconsense
endconsense-startconsense
write.nexus(matezillaconsensus,file=paste(Sys.Date(),"ConsensusPasserineTreeHackett4_", nTrees,"_mean-edge_ignore-absent.nex",sep=""))

