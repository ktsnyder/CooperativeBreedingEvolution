# Make various consensus trees


require(phytools)
multitree = read.tree("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/Birdsong - Life History Evolution/BirdzillaHackett4_Stage2_1000trees.tre")
passermultitree = read.nexus("/Users/kate/Desktop/PhyloBiology/PasserineMultiphy1000Hackett4_nondicho.nex")
OscineTree = read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

onepassertree <- read.nexus("2021-08-31ConsensusPasserineTreeEricson10_1000.nex")

datanames <- onepassertree$tip.label

matezillamultitree <- list()
for (i in 1:1000) {
  dropformatezilla <- set.seed(10)
  temptree <- thousandtrees[[i]]
  havedatavec <- temptree$tip.label %in% as.character(datanames) 
  dropformatezilla <- which(havedatavec == FALSE)
  matezillamultitree[[i]] <- drop.tip(phy = temptree,tip=dropformatezilla)
}
print(paste("done trimming trees", Sys.time()))



# consensus using mean.edge, ignore absent edges in edge means
startconsense <- Sys.time()
startconsense
matezillaconsensus <- consensus.edges(passermultitree, method = "mean.edge", if.absent = "ignore")
endconsense <- Sys.time()
endconsense
endconsense-startconsense
#matezilladicho <- multi2di(matezillaconsensus)
#matezilladicho$edge.length[matezilladicho$edge.length == 0] <- 0.000000000000001
write.nexus(matezillaconsensus,file=paste(Sys.Date(),"ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex",sep=""))





# consensus using least squares
startconsenseLS <- Sys.time()
startconsenseLS
matezillaconsensusLS <- ls.consensus(passermultitree)
endconsenseLS <- Sys.time()
endconsenseLS-startconsenseLS
write.nexus(matezillaconsensusLS,file=paste(Sys.Date(),"ConsensusPasserineTreeHackett4_1000_ls-consensus.nex",sep=""))