# scratch compare multitree and regular simmap overlap

# Regular
dummyFile = "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine .csv"
realFile = "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine .csv"
dummyRegdf = read.csv(dummyFile)
realRegdf = read.csv(realFile)


# MultiTree
realMultitreeFile = "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/HighConfidence_Coop FemaleSong_Agg01 300trees 20simsPerTree_Hackett4Oscine_multitree_All REAL.csv"
realMultidf = read.csv(realMultitreeFile)

par(mfrow = c(5,1))
hist(realMultidf$coopQ01, xlim = c(0,0.07))
mean(realRegdf$coopQ01)
abline(v=mean(realRegdf$coopQ01), col = "red", lwd = 3)

hist(realMultidf$coopQ10, xlim = c(0,0.07))
mean(realRegdf$coopQ10)
abline(v=mean(realRegdf$coopQ10), col = "red", lwd = 3)

hist(realMultidf$FSQAbsPres, xlim = c(0,0.07))
mean(realRegdf$FSQAbsPres)
abline(v=mean(realRegdf$FSQAbsPres), col = "red", lwd = 3)

hist(realMultidf$FSQPresAbs, xlim = c(0,0.07))
mean(realRegdf$FSQPresAbs)
abline(v=mean(realRegdf$FSQPresAbs), col = "red", lwd = 3)



treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
consTree = read.nexus(treefile)
consTreedi = multi2di(consTree)
oneTree = multitrees[[301]]
oneTreesub = drop.tip(oneTree, which(!oneTree$tip.label %in% consTree$tip.label))
oneTreesubdi = multi2di(oneTreesub)

length(consTree$edge.length)
length(oneTreesub$edge.length)
length(oneTreesubdi$edge.length)
consTree$Nnode
oneTreesub$Nnode
oneTreesubdi$Nnode
sum(consTree$edge.length)
sum(oneTreesub$edge.length)
sum(oneTreesubdi$edge.length)

treeEdges = NA
for (i in 1:1000) {
  if (i %in% seq(0,1000, by=100)) {print(i)}
  oneTree = multitrees[[i]]
  oneTreesub = drop.tip(oneTree, which(!oneTree$tip.label %in% consTree$tip.label))
  treeEdges[i] <- sum(oneTreesub$edge.length)
}
hist(treeEdges, xlim = c(25000, 43000))
abline(v = sum(consTree$edge.length), col = "red", lwd = 3)

newdatafile = "2024-06-07_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
PasserineData = read.csv(newdatafile)
newdata = PasserineData[which(PasserineData$species %in% OscineTree$tip.label),]

for (i in 1:length(colnames(newdata))) {
  tempcolumn = colnames(newdata)[i]
  if (length(unique(newdata[,tempcolumn])) < 5) {
    print(tempcolumn)
    print(table(newdata[,tempcolumn]))
  }
}
3859+635 
4245+440 
3940+554 

newdata[which(!is.na(newdata$AnyNoncoopEqualsNoncoop) & is.na(newdata$AnyCoopEqualsCoop)),]
