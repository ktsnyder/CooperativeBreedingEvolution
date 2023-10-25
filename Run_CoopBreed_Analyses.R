## Execute analysis code
## Coded by Kate T Snyder
## Created 8/18/2021
## Last Edited 6/29/2023
## Copied to Creanza Lab Server 6/30/2023 to run bayestraits Discrete

setwd("~/Desktop/CooperativeBreedingEvolution")

songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate", "Continuity")
#newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
#newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
#newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-01_CoopSongFS_RColumns.csv"
olddata = read.csv("2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv") # pre-cornwallis
newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data_R.csv"
dataNoSongless = read.csv(newdata)
dataNoSongless$X = NULL
#treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
tree = read.nexus(treefile)
write.tree(tree, file = "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nwk")
currentlabel <- "Hackett-Tie2Noncoop"
CBcolumn = "MeanCoopTie2Noncoop"
#columns = c("HighConfidence_FemaleSong", CBcolumn)
columns = c("FemaleSong_Agg01", CBcolumn)


# test ER/ARD brownie
source("findQrates.R")
#findQrates(columns = "CoopBreed", plot=TRUE, newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = "PasserineTreeEricson-MeanCoop")
findQrates(columns = columns, plot=TRUE, newtree = treefile, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
findQrates(columns = CBcolumn, plot=TRUE, newtree = treefile, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
findQrates(columns = columns[1], plot=TRUE, newtree = treefile, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)


#### Brownie ----
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
#below run with "2020-10-11ConsensusPasserineTreeHack100.nex" 8/19/21
#run with "birdzillatreeMaybeConsensus.nex" 8/24/21
songfeatures <- c("Syllable.rep.final", "Syllable.rep.max", "Syllable.rep.min", "Syll.song.min", "Syll.song.max", "Syll.song.final", "Song.rep.final", "Song.rep.min", "Song.rep.max", "Duration.final", "Interval.final", "Duration.min", "Duration.max", "Interval.min", "Interval.max", "Song.rate","Continuity")
#CBcolumn = "HighConfidence_FemaleSong"
CBcolumn = "MeanCoopTie2Noncoop"

#discreteCatLabels = c("Female Song Absent", "Female Song Present")
discreteCatLabels = c("Noncooperative", "Cooperative")
currentlabel <- "_Hackett"
nsim = 100
for (k in 7) {
  print(Sys.time())
  newdata = newdata
  treefile = treefile
  currentlabel <- currentlabel
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c(CBcolumn, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  
  if (file.exists(paste0("OutputFiles/",Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"))) {
    print("file exists")
    plotbrownie(data = paste0(Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
  } else if (file.exists(paste0("OutputFiles/",Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"))) {
    plotbrownie(data = paste0(Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
    print("yesterday's file exists")
  } else {
    print("file does not exist")
    print(paste0("OutputFiles/",Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"))
  }
}


#### Scatterboxes, ACEtree ----
source("scatterboxes.R")
#scatterboxes(DiscreteTrait = CBcolumn, newdata = newdata, newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", otherlabel = "AnyCoop")
scatterboxes(DiscreteTrait = CBcolumn, newdata = newdata, newtree = treefile, otherlabel = currentlabel)


source("plotACEtree.R")
# Song Features
for (k in 1) {
  feature <- songfeatures[k]
plotACEtree(columns = c(CBcolumn, feature), cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = newdata, newtree = newtree, islog = feature, discretelabels = c("Non-cooperative","Cooperative"), discretemodel = "ARD", otherlabel = currentlabel)
} 

# Cooperative Breeding and Female Song
columns = c("MeanCoopTie2Coop", "HighConfidence_FemaleSong")
columns = c("AnyCoopEqualsCoop", "FemaleSong_Agg01")
#currentlabel <- c("Tie2Noncoop FSHighConf Hackett")
subsetout <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile)
subsetdf = subsetout$subsetdf
subsettree = subsetout$subsettree
#findQrates(columns = columns[1], plot=TRUE, newtree = subsettree, newdata = subsetdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
#findQrates(columns = columns[2], plot=TRUE, newtree = subsettree, newdata = subsetdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
#plotACEtree(columns = "MeanCoopTie2Noncoop", cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = subsetdf, newtree = subsettree, discretelabels = c("Non-cooperative","Cooperative"), discretemodel = "ARD", otherlabel = currentlabel)
#plotACEtree(columns = "HighConfidence_FemaleSong", cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = subsetdf, newtree = subsettree, discretelabels = c("Female Song Absent","Female Song Present"), discretemodel = "ARD", otherlabel = currentlabel)

# Double-tip ACE trees
whichnodes = "FS"  #"FS" # "Coop" # "no"
tipsize = 0.1
filename = paste(columns[1], columns[2], "fan phylo double tips", whichnodes, "nodes.pdf")
pdf(filename, height = 8, width = 9)
plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
py = c("black","orange")
treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf$MeanCoopTie2Coop == 1)]
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=tipsize, offset = 1)
py2 = c("blue","red")
treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf$HighConfidence_FemaleSong == 1)]
tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=tipsize, offset = 2)
if (whichnodes == "Coop") {
  discretetraitvec <- subsetdf[,columns[1]]
  names(discretetraitvec) <- subsetdf[,1]
  pynodes = py
  circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
  nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
} else if (whichnodes == "FS") {
  discretetraitvec <- subsetdf[,columns[2]]
  names(discretetraitvec) <- subsetdf[,1]
  pynodes = py2
  circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
  nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
} 
allpy = c(py, py2)
alllabs = c("Noncooperative", "Cooperative", "Female Song Absent", "Female Song Present")
legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")
dev.off()


#### Simple bayestraits discrete tests for Female Song and CoopBreed ----
source("btwDiscreteKTS.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4"
#.BayesTraitsPath = "~/Documents/BayesTraitsV3"
dataIn = read.csv(newdata)
dataIn$X = NULL
CBcolumn = "Kin_NK"
columns = c(CBcolumn, "FemaleSong_Agg01")
#currentlabel <- "Hackett-TieNoncoop-FSAgg"
subsetbtw <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = TRUE)
subsettree <- subsetbtw$subsettree
subsetdf <- subsetbtw$subsetdf
subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])


## Using my altered btw::bayestraits function
source("btwV2bayestraitsKTS.R")
seeds = 101:350
Version = "V4"
Method = "ML"
MLtries = 100
currentlabel = paste0(currentlabel, " mlt", MLtries, "_noRes") # no restrictions

outputdf = set.seed(10)
# Do Independent first
for (i in seeds) {
  print(paste("Independent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Model = "Independent"
  commandVector = c("2", "1", paste("mlt", MLtries), paste("Se", i)) # mlt number of tries
  
  outInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = Version, remove_files = T, BTdirpath = "~/Documents")
  resultsInd = outInd$Log$results
  temprow = cbind(Seed, Version, Model, Method, MLtries, resultsInd)
  tempdf = rbind(tempdf, temprow)
  
  outdf = tempdf[,c("Seed","Version","Model","Method","MLtries","Tree.No", "Lh")]
  outdf$q12 = tempdf$alpha2
  outdf$q13 = tempdf$alpha1
  outdf$q21 = tempdf$beta2
  outdf$q24 = tempdf$alpha1
  outdf$q31 = tempdf$beta1
  outdf$q34 = tempdf$alpha2
  outdf$q42 = tempdf$beta1
  outdf$q43 = tempdf$beta2
  outdf = cbind(outdf, tempdf[,c("Root...P.0.0.", "Root...P.0.1.", "Root...P.1.0.", "Root...P.1.1.")])
  
  outputdf = rbind(outputdf,outdf)
}
write.csv(outputdf, paste0("BayesTraitsDiscrete_", currentlabel, ".csv"))

# Then do dependent
for (i in seeds[1:length(seeds)]) {
  print(paste("Dependent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Model = "Dependent"
  
  commandVector = c("3", "1", paste("mlt", MLtries), paste("Se", i))
  
  outDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = Version, remove_files = T, BTdirpath = "~/Documents")
  resultsDep = outDep$Log$results
  temprow = cbind(Seed, Version, Model, Method, MLtries, resultsDep)
  tempdf = rbind(tempdf, temprow)
  
  outputdf = rbind(outputdf, tempdf)
  
  if (Seed %in% c(110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 240, 260, 280, 300, 350, 400, 450, 500)) {
    write.csv(outputdf, paste0("BayesTraitsDiscrete_", currentlabel, ".csv"))
  }
}
write.csv(outputdf, paste0("BayesTraitsDiscrete_", currentlabel, ".csv"))



## Using btw V2 - run 5/31/2023
require("btw")
setwd("/Users/kate/Documents")
tree = subsettree
df = subsetdf
currentlabel <- "V2-Hackett-TieNoncoop-FSAgg-q31q42Lim15"
simplebtwV2Out <- set.seed(10)
nsim = 50

for (n in 1:nsim) {
  print(paste("v2 #", n))
  commandIndML <- c("2","1", "mlt 10")#, "res q31 q42 15") # replace nocorrD
  IndMLout <- bayestraits(df,tree,commandIndML)
  
  commandDepML <- c("3","1", "mlt 10")#, "res q31 q42 15")  # replace corrD
  DepMLout <- bayestraits(df,tree,commandDepML) # , silent = FALSE
  
  lrtestresults <- lrtest(DepMLout, IndMLout)
  LRstat <- lrtestresults$LRstat
  LRpval <- lrtestresults$pval
  corrD = DepMLout$Log$results
  nocorrD = IndMLout$Log$results
  transandp <- cbind(corrD,LRstat,LRpval, nocorrD)
  simplebtwV2Out <- rbind(simplebtwV2Out, transandp)
}
means <- apply(X = simplebtwV2Out,MARGIN = 2,FUN = mean)
meansdf <- as.data.frame(rbind(means,means))
pvalMed <- median(simplebtwV2Out$LRpval)
pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
plotdiscrete(meansdf[1,3:10], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed)) # plotdiscrete from btwV1 now in btwDiscreteKTS.R
dev.off()


# Using simplebtwDiscrete.R - UNFINISHED
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = dfnew, newtree = treefile, treelabel = treelabel, nsim = nsim, savecsvs = TRUE, KeepBTInputFiles = TRUE)
simplebtwOutput1 = simplebtwOutput[]
plotDiscreteBayes(columns=columns, simplebtwOut = simplebtwOutput, nsim = 100, newpdf = T, treelabel = "Hackett", )


#### BayesTraits song features ----
#source("btwfunction.R")
source("~/Desktop/CooperativeBreedingEvolution/btw2function_DiscreteKTS.R")
source("BayesPlots_choosebin.R")
songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate", "Song.rep.min", "Song.rep.max")
#newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
#newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
#treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
currentlabel <- "PasserineTreeHackett-OmitTies"
CBcolumn = "MeanCoopOmitTies"
nsim = 100
for (k in c(7:11, 16)) { 
  feature <- songfeatures[k]
  btwfunction(columns = c(CBcolumn, feature), plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = dataNoSongless)
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"Bayes",CBcolumn,feature, nsim, "reps.csv") 
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  #colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = CBcolumn,SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3, arrowmod = 0.05)
} 

#bayestraitsKTS(columns = columns, plot = FALSE, )
CBcolumn = "MeanCoopTie2Noncoop"
feature = "Song.rep.final"
filename <- paste0(Sys.Date()-9,"Bayes",CBcolumn,feature, nsim, "reps.csv") 
BTdf <- read.csv(filename)
transitionBinplots(MateParam = CBcolumn, SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3, arrowmod = 0.05)



#### Brownie - Do sylls/song brownie with fake CoopBreed data (randomly assign 16 species as Coop) ----
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
feature <- "Syll.song.final"
subset <- subsettreedata(columns = "CoopBreed", newtree = "birdzillatreeMaybeConsensus.nex", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
pdf(paste0(Sys.Date(),"Brownie_Syll.song.final_fakeCoopData5.pdf"), width = 12, height = 12)
par(mfrow = c(3,3))
columns <- c("CoopBreed", feature)
nsims = 100
for (j in 13:20) {
  set.seed(j)
  df <- subset$subsetdf
  fakecoops <- sample(which(!is.na(df$Syll.song.final)), 16)
  df$CoopBreed[which(!is.na(df$Syll.song.final))] <- rep(0, times = length(which(!is.na(df$Syll.song.final))))
  df$CoopBreed[fakecoops] <- 1
 # print(df[,c("Syll.song.final", "CoopBreed")])
  browniefunction(columns = c("CoopBreed", feature), newdata = df, newtree = "birdzillatreeMaybeConsensus.nex", nsim = nsims, islog = feature, plotsimmaps = FALSE, otherlabel = paste0("fakeCoopData", j))
plotbrownie(data = paste0(Sys.Date(),columns[1], columns[2],"fakeCoopData",j, "_brownie",nsims,"sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = "fakeCoopData", newpdf = FALSE, nsim = nsims, islog = TRUE)

# plot tree
  fakedatasubset <- subsettreedata(columns = c("CoopBreed", feature), newtree = "birdzillatreeMaybeConsensus.nex", newdata = df)
  faketree <- fakedatasubset$subsettree
  fakedata <- fakedatasubset$subsetdf
  discretetraitvector <-fakedata[,columns[1]]
  names(discretetraitvector) <- fakedata[,1] 
treetiplabels <- faketree$tip.label %in% names(discretetraitvector[discretetraitvector==1]) 
plot(faketree, cex = 0.25)
py <- c("black", "red")
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.3)
}
dev.off()


### Re-plot fake data brownies with real data first and example simmaps
subset <- subsettreedata(columns = c("CoopBreed", feature), newtree = "birdzillatreeMaybeConsensus.nex", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
df <- subset$subsetdf
pdf(paste0(Sys.Date(),"Brownie_Syll.song.final_RealPlusFakeCoopData20reps.pdf"), width = 12, height = 12)
par(mfrow = c(3,3))
columns <- c("CoopBreed", feature)
py = c("black","red")
pynamed <- py
names(pynamed) <- c(0,1)

#Real Data
plotbrownie(data = "2021-08-27CoopBreedSyll.song.finalRealData_UpdatedQrates_brownie500sim.csv", columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = paste0("Real Coop Data"), newpdf = FALSE, nsim = 500, islog = TRUE)

simmapQset <-read.simmap("OutputFiles/2021-08-27CoopBreedSyll.song.finalRealData_UpdatedQrates500simmaps.txt",format="phylip",version=1)
plotSimmap(simmapQset,fsize=0.2,lwd=0.8, colors = pynamed)
#add tips
discretetraitvec <- df[,columns[1]]
names(discretetraitvec) <- df[,1]
treetiplabels <- simmapQset$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
Qout <- findQrates(columns = columns, plot = FALSE, newtree = "birdzillatreeMaybeConsensus.nex", newdata = df, otherlabel = "RealData")
numrates <- Qout$qrates
title(main=paste(" ","\nE.g. ARDmodel","Qrates (output for Brownie):",numrates[2],numrates[3]),cex.main = 0.5)

#Fake Data
nsims = 100
for (j in c(13:20,1:12)) {
plotbrownie(data = paste0(Sys.Date(),columns[1], columns[2],"fakeCoopData",j, "_brownie",nsims,"sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = paste0("Fake Coop Data"), newpdf = FALSE, nsim = nsims, islog = TRUE)
  
set.seed(j)
df <- subset$subsetdf
fakecoops <- sample(which(!is.na(df$Syll.song.final)), 16)
df$CoopBreed[which(!is.na(df$Syll.song.final))] <- rep(0, times = length(which(!is.na(df$Syll.song.final))))
df$CoopBreed[fakecoops] <- 1
fakesimmap <-read.simmap(paste0("OutputFiles/",Sys.Date(),columns[1], columns[2],"fakeCoopData",j,nsims,"simmaps.txt"),format="phylip",version=1)
plotSimmap(fakesimmap,fsize=0.2,lwd=0.8, colors = pynamed)
#add tips
discretetraitvec <- df[,columns[1]]
names(discretetraitvec) <- df[,1]
treetiplabels <- fakesimmap$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
Qout <- findQrates(columns = columns, plot = FALSE, newtree = "birdzillatreeMaybeConsensus.nex", newdata = df, otherlabel = "FAKEDATA")
numrates <- Qout$qrates
title(main=paste(" ","\nE.g. ARDmodel","Qrates (output for Brownie):",numrates[2],numrates[3]),cex.main = 0.5)
}



#### Misc snippets from along the way ----
# 8/31/2021
songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Syll.song.min", "Syll.song.max")
newdata = "~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2021-08-18CoopSong_MeanCoop_NatCommsSubset.csv"

source("subsettreedata.R")
source("browniefunction.R")
#below run with "2020-10-11ConsensusPasserineTreeHack100.nex" 8/19/21
#run with "birdzillatreeMaybeConsensus.nex" 8/24/21
# songfeatures <- c("Syll.song.min", "Syll.song.max")
for (k in 2) {
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c("CoopBreed", feature), newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", nsim = 500, islog = feature, plotsimmaps = TRUE,otherlabel = "Ericson10Consensus" )
} 

source("plotbrownie.R")
for (i in 2) {
  print(i)
  feature <- songfeatures[i]
  plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, "Ericson10Consensus_brownie500sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = "Ericson10Consensus", newpdf = TRUE, nsim = 500, islog = TRUE)
}

# Brownie continued - Change each "Cooperative" species to "Noncooperative" iteratively to test robustness of results
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
subset <- subsettreedata(columns = c("CoopBreed", "Syll.song.final"), newtree = "birdzillatreeMaybeConsensus.nex", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", islog = "Syll.song.final")
subsetdf <- subset$subsetdf
subsettree <- subset$subsettree
cooperativebirds <- subsetdf$species[which(subsetdf$CoopBreed == 1)]
feature <- "Syll.song.final"

pdf("Brownie_Syllsongfinal_eachcoopSwitched_500.pdf", width = 8, height = 12)
par(mfrow = c(4,2))
for (coopbird in cooperativebirds) {
  #thisdf <- subset$subsetdf
  #thisdf <- read.csv("2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
  thisdf[which(thisdf$species == coopbird), "CoopBreed"] <- 0
  browniefunction(columns = c("CoopBreed", feature), newdata = thisdf, newtree = treefile, nsim = 500, islog = feature, plotsimmaps = FALSE)
  plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, "brownie500sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = paste0(coopbird,"NonCoop"), newpdf = FALSE, nsim = 500, islog = TRUE)
}
dev.off()


## Using btw V1 - run 3/9/2022; replaced with function names from btwDiscreteKTS.R 6/1/23
simplebtwV1Out <- set.seed(10)
currentlabel <- "V1-Hackett-TieNoncoop-FSAgg"
nsim = 50
for (n in 1:nsim) {
  print(paste("v1 #", n))
  nocorrD <- DiscreteKTS(subsettree, subsetdf)
  corrD <- DiscreteKTS(subsettree, subsetdf, dependent=TRUE)
  lrtestresults <- lrtestV1(corrD, nocorrD)
  tempRow <- cbind(corrD, lrtestresults, nocorrD)
  simplebtwV1Out <- rbind(simplebtwV1Out, tempRow)
}
simplebtwV1Out <- as.data.frame(simplebtwV1Out)
write.csv(simplebtwV1Out, file = paste0(Sys.Date(), "_BayesTraits_", currentlabel, "_", nsim, "sims.csv"))
means <- apply(X = simplebtwV1Out,MARGIN = 2,FUN = mean)
meansdf <- as.data.frame(rbind(means,means))
#meansdf <- as.data.frame(as.matrix(means))
pvalMed <- median(simplebtwV1Out$pval)
pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
#plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwV1Out, nsim = nsim, newpdf = T)
dev.off()

simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, nsim = nsim, )



#### Cycle Simple BT Discrete Jackknifed ----
# Simple bayestraits discrete tests for Female Song and CoopBreed

library(devtools)
install_github("rgriff23/btw")
library(btw)

newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"

source("btwDiscreteKTS.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4"
#.BayesTraitsPath = "~/Documents/BayesTraitsV3"
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumn = "MeanCoopTie2Noncoop"
columns = c(CBcolumn, "FemaleSong_Agg01")
subset1 <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf
subsetdf1[,columns[1]] <- as.character(subsetdf1[,columns[1]])
subsetdf1[,columns[2]] <- as.character(subsetdf1[,columns[2]])

# unique(subsetdf1$Order3_BirdtreeMatchSpecies2)
# familyvec = unique(subsetdf1$Family3_BirdtreeMatchSpecies2)
# familyvec = na.omit(familyvec)

require(dplyr)
familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2[which(familycounts$n > 2)]

for (i in 1:length(familyvec)) {
  
  familyToRemove = familyvec[i]
  currentlabel = paste0("remove",familyToRemove)
  print(currentlabel)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 != familyToRemove),]

subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
subsettree <- subsetbtw$subsettree
subsetdf <- subsetbtw$subsetdf
subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
numSpecies = length(subsetdf$species)

## Using my altered btw::bayestraits function
source("btwV2bayestraitsKTS.R")
seeds = 121:140
Version = "V4"
Method = "ML"
MLtries = 100

outputdf = set.seed(10)
# Do Independent first
for (i in seeds) {
  print(paste("Independent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Model = "Independent"
  commandVector = c("2", "1", paste("mlt", MLtries), paste("Se", i)) # mlt number of tries
  
  outInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, BTversionNum = Version, remove_files = T, BTdirpath = "~/Documents")
  resultsInd = outInd$Log$results
  temprow = cbind(Seed, Version, Model, Method, MLtries, numSpecies, resultsInd)
  tempdf = rbind(tempdf, temprow)
  
  outdf = tempdf[,c("Seed","Version","Model","Method","MLtries", "numSpecies", "Tree.No", "Lh")]
  outdf$q12 = tempdf$alpha2
  outdf$q13 = tempdf$alpha1
  outdf$q21 = tempdf$beta2
  outdf$q24 = tempdf$alpha1
  outdf$q31 = tempdf$beta1
  outdf$q34 = tempdf$alpha2
  outdf$q42 = tempdf$beta1
  outdf$q43 = tempdf$beta2
  outdf = cbind(outdf, tempdf[,c("Root...P.0.0.", "Root...P.0.1.", "Root...P.1.0.", "Root...P.1.1.")])
  
  outputdf = rbind(outputdf,outdf)
}
write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", currentlabel, ".csv"))

# Then do dependent
for (i in seeds[1:length(seeds)]) {
  print(paste("Dependent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Model = "Dependent"
  
  commandVector = c("3", "1", paste("mlt", MLtries), paste("Se", i))
  
  outDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, BTversionNum = Version, remove_files = T, BTdirpath = "~/Documents")
  resultsDep = outDep$Log$results
  temprow = cbind(Seed, Version, Model, Method, MLtries, numSpecies, resultsDep)
  tempdf = rbind(tempdf, temprow)
  
  outputdf = rbind(outputdf, tempdf)
  
  if (Seed %in% c(110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 240, 260, 280, 300, 350, 400, 450, 500)) {
    write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", currentlabel, ".csv"))
  }
}
write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", currentlabel, ".csv"))

} # end cycle through familyvec




#### Cycle Simple BT Discrete Single Family ----
# Simple bayestraits discrete tests for Female Song and CoopBreed

library(devtools)
install_github("rgriff23/btw")
library(btw)

newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"

source("btwDiscreteKTS.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4"
#.BayesTraitsPath = "~/Documents/BayesTraitsV3"
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumn = "MeanCoopTie2Noncoop"
columns = c(CBcolumn, "FemaleSong_Agg01")
subset1 <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf
subsetdf1[,columns[1]] <- as.character(subsetdf1[,columns[1]])
subsetdf1[,columns[2]] <- as.character(subsetdf1[,columns[2]])

# unique(subsetdf1$Order3_BirdtreeMatchSpecies2)
# familyvec = unique(subsetdf1$Family3_BirdtreeMatchSpecies2)
# familyvec = na.omit(familyvec)

require(dplyr)
familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2[which(familycounts$n > 4)]

for (j in 1:length(familyvec)) {
  
  familyToKeep = familyvec[j]
  templabel = paste0(columns[1], "_", columns[2], "_", "only",familyToKeep)
  print(paste(j, templabel))
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 == familyToKeep),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  if ( length(unique(subsetdf[,columns[1]])) == 2 & length(unique(subsetdf[,columns[2]])) == 2 ) {
    print(paste(familyToKeep, "has both FS/noFS and CB/nonCB"))

  
  ## Using my altered btw::bayestraits function
  source("btwV2bayestraitsKTS.R")
  seeds = 141:190
  Version = "V4"
  Method = "ML"
  MLtries = 100
  
  outputdf = set.seed(10)
  # Do Independent first
  for (i in seeds) {
    print(paste("Independent, Seed:", i))
    tempdf = set.seed(i)
    Seed = i
    Model = "Independent"
    commandVector = c("2", "1", paste("mlt", MLtries), paste("Se", i)) # mlt number of tries
    
    outInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, BTversionNum = Version, remove_files = T, BTdirpath = "~/Documents")
    resultsInd = outInd$Log$results
    temprow = cbind(Seed, Version, Model, Method, MLtries, numSpecies, resultsInd)
    tempdf = rbind(tempdf, temprow)
    
    outdf = tempdf[,c("Seed","Version","Model","Method","MLtries", "numSpecies", "Tree.No", "Lh")]
    outdf$q12 = tempdf$alpha2
    outdf$q13 = tempdf$alpha1
    outdf$q21 = tempdf$beta2
    outdf$q24 = tempdf$alpha1
    outdf$q31 = tempdf$beta1
    outdf$q34 = tempdf$alpha2
    outdf$q42 = tempdf$beta1
    outdf$q43 = tempdf$beta2
    outdf = cbind(outdf, tempdf[,c("Root...P.0.0.", "Root...P.0.1.", "Root...P.1.0.", "Root...P.1.1.")])
    
    outputdf = rbind(outputdf,outdf)
  }
  write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", templabel, ".csv"))
  
  # Then do dependent
  for (i in seeds[1:length(seeds)]) {
    print(paste("Dependent, Seed:", i))
    tempdf = set.seed(i)
    Seed = i
    Model = "Dependent"
    
    commandVector = c("3", "1", paste("mlt", MLtries), paste("Se", i))
    
    outDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, BTversionNum = Version, remove_files = T, BTdirpath = "~/Documents")
    resultsDep = outDep$Log$results
    temprow = cbind(Seed, Version, Model, Method, MLtries, numSpecies, resultsDep)
    tempdf = rbind(tempdf, temprow)
    
    outputdf = rbind(outputdf, tempdf)
    
    if (Seed %in% c(110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 240, 260, 280, 300, 350, 400, 450, 500)) {
      write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", templabel, ".csv"))
    }
  }
  write.csv(outputdf, paste0("BayesTraitsDiscreteML_", columns[1], "-", columns[2], "_", templabel, ".csv"))
  
  } # end if family has instances of CB = 0,1 and FS = 0,1
  
} # end cycle through familyvec for single family tests


#### Cycle plot simple bayestraits ----
# use plotSimpleDiscreteBayes from plotSimpleDiscreteBayes.R

trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"

# EITHER
tempfolder = paste0("BayesTraitsDiscreteML_", trait1, "-", trait2, " Jackknife Outputs")
filelist = list.files(tempfolder)
filelist = filelist[which(str_detect(filelist, "remove") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
familytreatment = "removed"
pdf(file = paste0("jackknifed BayesTraitsDiscrete ", trait1, " ", trait2, ".pdf"), width = 12, height = 10)

# OR
# (single family)
tempfolder = paste0("BayesTraitsDiscreteML_", trait1, "-", trait2, " SingleFamily Outputs")
filelist = list.files(tempfolder)
filelist = filelist[which(str_detect(filelist, "only") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
familytreatment = "only"
pdf(file = paste0("SingleFamily BayesTraitsDiscrete ", trait1, " ", trait2, ".pdf"), width = 12, height = 10)

par(mfrow=c(3,3), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
    font.main=1, cex.main=1.25, cex.lab=1.1, cex.axis=1.1)

for (j in 1:length(filelist)) {
tempfile = filelist[j]

if (str_detect(tempfile, "remove")) {
  tempFamily <- unique(gsub(".*remove([A-Za-z]+)\\.csv", "\\1", tempfile))
} else if (str_detect(tempfile, "only")) {
  tempFamily <- unique(gsub(".*only([A-Za-z]+)\\.csv", "\\1", tempfile))
}

tempfilepath = paste(tempfolder, tempfile, sep = "/")

dfIn = read.csv(tempfilepath)
dfInDep = dfIn[which(dfIn$Model == "Dependent"),]
dfInInd = dfIn[which(dfIn$Model == "Independent"),]

Nspecies = dfIn[1,"numSpecies"]

templabel = paste0(familytreatment, tempFamily, " Nspecies=", Nspecies)

tryCatch(
  {
plotSimpleDiscreteBayes(columns = c(trait1,trait2), df = dfInDep, nocorrDdf = dfInInd, LhCol = "Lh", nsim = NULL, treelabel = NULL, newpdf = FALSE, otherlabel = "", ylabel = templabel, arrowmod = 1, roundDigits = 3)
  }, error=function(e) {
    message(paste("An Error occurred, could not plot", templabel))
  }) # end trycatch

}
dev.off()


#### Cycle plot simple bayestraits AND simmap overlap ----
# For each family, plot 1) Simmap overlap Real D hist, 2) dummy D hist, 3) Observed state boxplots, 4) bayestraits Dependent plot, 5) BayesTraits Pval plot, 6) family tree with CB tips, 7) family tree with FS tips, 8?) whole tree with family indicated w simmap

newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
newdatadf = read.csv(newdata)
newdatadf$X = NULL
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"

# get which tests were done per family
alltestdf = set.seed(10)

# BayesTraits jackknife
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"
trait1 = "MeanCoopTie2Coop"
trait2 = "FemaleSong_Agg01"
tempfolder = paste0("BayesTraitsDiscreteML_", trait1, "-", trait2, " Jackknife Outputs")
filelist = list.files(tempfolder)
filelist = filelist[which(str_detect(filelist, "remove") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempFamily <- unique(gsub(".*remove([A-Za-z]+)\\.csv", "\\1", filelist))
test = "BayesTraits Jackknife"
famdf <- cbind(test, trait1, trait2, tempFamily)
alltestdf <- rbind(alltestdf, famdf)
alltestdf <- as.data.frame(alltestdf)

# BayesTraits single family
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"
tempfolder = paste0("BayesTraitsDiscreteML_", trait1, "-", trait2, " SingleFamily Outputs")
filelist = list.files(tempfolder)
filelist = filelist[which(str_detect(filelist, "only") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempFamily <- unique(gsub(".*only([A-Za-z]+)\\.csv", "\\1", filelist))
test = "BayesTraits SingleFamily"
famdf <- cbind(test, trait1, trait2, tempFamily)
alltestdf <- rbind(alltestdf, famdf)

# Simmap Overlap jackknife
trait1 = "MeanCoopTie2Coop"
trait2 = "FemaleSong_Agg01"
require(reshape2)
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"

trait1 = "MeanCoopTie2Coop"
trait2 = "HighConfidence_FemaleSong"
filelist = list.files(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/Jackknife ", trait1, " ", trait2), recursive = T, pattern = ".csv", full.names = T)
tempFiles = jackknifeFiles = filelist[which(str_detect(filelist, "remove") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempFamily <- tempUniqueFamilies <- unique_remove_strings <- unique(gsub(".*remove([^ ]*) .*", "\\1", jackknifeFiles))
test = "SimmapOverlap Jackknife"
famdf <- cbind(test, trait1, trait2, tempFamily)
alltestdf <- rbind(alltestdf, famdf)
familytreatment = "removed"

## Simmap overlap single family
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"

filelist = list.files(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/Single Family ", trait1, " ", trait2), recursive = T, pattern = ".csv", full.names = T)
tempFiles = onefamilyFiles = filelist[which(str_detect(filelist, "only") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempFamily <- tempUniqueFamilies <- unique_only_strings <- unique(gsub(".*only([^ ]*) .*", "\\1", onefamilyFiles))
test = "SimmapOverlap SingleFamily"
famdf <- cbind(test, trait1, trait2, tempFamily)
alltestdf <- rbind(alltestdf, famdf)
familytreatment = "only"
#write.csv(alltestdf, "summary of jackknife singlefamily BayesTraits SimmapOverlap tests.csv")

famcounts = alltestdf %>% group_by(tempFamily) %>% count
famcounts %>% group_by(n) %>% count
testcounts = alltestdf %>% group_by(test, trait1, trait2) %>% count


require(reshape2)
tempcount = 1
if (tempcount == 7) {
  pdf(file = "plots by family - SimmapOverlap BayesTraits Jackknife SingleFamily - 7tests.pdf", height = 12, width = 30)
  par(mfcol=c(3,8), mai=c(1,1,0.8,0.3), oma=c(3,2,3,2),
    font.main=1, cex.main=1.25, cex.lab=1.1, cex.axis=1.1)
  families = famcounts$tempFamily[which(famcounts$n == tempcount)]
  
  for (i in 1:length(families)) {
    tempfamily = families[i]
    
    # plot single-family trees
    subsetout = subsettreedata(columns = c(temprow$trait1, temprow$trait2), cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladesubsetvalue = tempfamily, newdata = newdatadf, newtree = treefile)
    subsettree = subsetout$subsettree
    subsetdf = subsetout$subsetdf
    
    # no tips
    plot.phylo(subsettree, type = "f", show.tip.label = TRUE, cex = 0.2)
    title(main = paste(tempfamily))
    
    # trait1 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py = c("blue","red")
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait1] == 1)]
    tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=1)
    legend("bottomleft", legend = c("Noncooperative", "Cooperative"), cex = 0.9, fill=py, bty="n")
    
    # trait2 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py2 = c("black","orange")
    treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait2] == 1)]
    tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=1)
    legend("bottomleft", legend = c("Female Song Absent", "Female Song Present"), cex = 0.9, fill=py2, bty="n")
    
    
    for (j in 1:7) {
      temprow = testcounts[j,]
      
      if (str_detect(temprow$test, "BayesTraits")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " Jackknife Outputs")
          familytreatment = "remove"
        } else {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " SingleFamily Outputs")
          familytreatment = "only"
        }
        tempfile = list.files(tempfolder, pattern = tempfamily)
        # bayestraits plots
        dfIn = read.csv(paste0(tempfolder,"/",tempfile))
        dfInDep = dfIn[which(dfIn$Model == "Dependent"),]
        dfInInd = dfIn[which(dfIn$Model == "Independent"),]
        
        Nspecies = dfIn[1,"numSpecies"]
        
        templabel = paste0(familytreatment, tempfamily, " Nspecies=", Nspecies)
        
        tryCatch(
          {
            plotSimpleDiscreteBayes(columns = c(temprow$trait1, temprow$trait2), df = dfInDep, nocorrDdf = dfInInd, LhCol = "Lh", nsim = NULL, treelabel = NULL, newpdf = FALSE, otherlabel = "", ylabel = templabel, arrowmod = 1, roundDigits = 3)
          }, error=function(e) {
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            message(paste("An Error occurred, could not plot", templabel))
          }) # end trycatch
      } # end if BayesTraits
      
      
      
      if (str_detect(temprow$test, "SimmapOverlap")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste("Simmap Overlap Outputs/Jackknife", temprow$trait1, temprow$trait2)
          familytreatment = "remove"
        } else {
          tempfolder = paste("Simmap Overlap Outputs/Single Family", temprow$trait1, temprow$trait2)
          familytreatment = "only"
        }
        tempfiles = list.files(tempfolder, pattern = tempfamily)
        tempRealFile = tempfiles[which(str_detect(tempfiles, "REAL"))]
        tempDummyFile = tempfiles[which(str_detect(tempfiles, "DUMMY"))]
        RealDF = read.csv(paste0(tempfolder,"/",tempRealFile))
        DummyDF = read.csv(paste0(tempfolder,"/",tempDummyFile))
        
        # plot simmap overlap plots
        templabel = paste0(familytreatment, tempfamily)
        
        calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = templabel)
        
        # Third plot: observed state proportions
        dfDummyMelt <- melt(RealDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
        
        # Create the basic boxplot without x-axis labels (xaxt = "n") and without outlines
        boxplot(ObservedState.prop ~ ObservedState, data = dfDummyMelt, xlab = "", ylab = "Observed State Proportion", xaxt = "n", outline = FALSE, boxwex = 0.5, col = "lightgray")
        
        # Add points to the plot
        points(jitter(as.numeric(dfDummyMelt$ObservedState), amount = 0.05), dfDummyMelt$ObservedState.prop, pch = 16, col = "darkred", cex = 0.6)
        
        # Add rotated x-axis labels
        axis(1, at = 1:length(unique(dfDummyMelt$ObservedState)), 
             labels = FALSE)
        text(1:length(unique(dfDummyMelt$ObservedState)), 
             par("usr")[3] - 0.001, srt = 45, adj = 1.2, 
             labels = as.character(unique(dfDummyMelt$ObservedState)), xpd = TRUE, cex=0.8)
      } # end if SimmapOverlap
      
    } # end for j in 1:7
    
  } # end for i in 1: length(families)
  dev.off()
} else if (tempcount == 5) {
  pdf(file = "plots by family - SimmapOverlap BayesTraits Jackknife SingleFamily - 5tests.pdf", height = 12, width = 25)
  par(mfcol=c(3,6), mai=c(1,1,0.8,0.3), oma=c(3,2,3,2),
      font.main=1, cex.main=1.25, cex.lab=1.1, cex.axis=1.1)
  families = famcounts$tempFamily[which(famcounts$n == tempcount)]
  for (i in 1:length(families)) {
    tempfamily = families[i]
    
    # plot single-family trees
    subsetout = subsettreedata(columns = c(temprow$trait1, temprow$trait2), cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladesubsetvalue = tempfamily, newdata = newdatadf, newtree = treefile)
    subsettree = subsetout$subsettree
    subsetdf = subsetout$subsetdf
    
    # no tips
    plot.phylo(subsettree, type = "f", show.tip.label = TRUE, cex = 0.2)
    title(main = paste(tempfamily))
    
    # trait1 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py = c("blue","red")
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait1] == 1)]
    tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=1)
    legend("bottomleft", legend = c("Noncooperative", "Cooperative"), cex = 0.9, fill=py, bty="n")
    
    # trait2 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py2 = c("black","orange")
    treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait2] == 1)]
    tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=1)
    legend("bottomleft", legend = c("Female Song Absent", "Female Song Present"), cex = 0.9, fill=py2, bty="n")
    
    
    for (j in 1:5) {
      testcounts5 = testcounts[which(testcounts$n > 30),]
      temprow = testcounts5[j,]
      
      if (str_detect(temprow$test, "BayesTraits")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " Jackknife Outputs")
          familytreatment = "remove"
        } else {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " SingleFamily Outputs")
          familytreatment = "only"
        }
        tempfile = list.files(tempfolder, pattern = tempfamily)
        # bayestraits plots
        dfIn = read.csv(paste0(tempfolder,"/",tempfile))
        dfInDep = dfIn[which(dfIn$Model == "Dependent"),]
        dfInInd = dfIn[which(dfIn$Model == "Independent"),]
        
        Nspecies = dfIn[1,"numSpecies"]
        
        templabel = paste0(familytreatment, tempfamily, " Nspecies=", Nspecies)
        
        tryCatch(
          {
            plotSimpleDiscreteBayes(columns = c(temprow$trait1, temprow$trait2), df = dfInDep, nocorrDdf = dfInInd, LhCol = "Lh", nsim = NULL, treelabel = NULL, newpdf = FALSE, otherlabel = "", ylabel = templabel, arrowmod = 1, roundDigits = 3)
          }, error=function(e) {
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            message(paste("An Error occurred, could not plot", templabel))
          }) # end trycatch
      } # end if BayesTraits
      
      
      
      if (str_detect(temprow$test, "SimmapOverlap")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste("Simmap Overlap Outputs/Jackknife", temprow$trait1, temprow$trait2)
          familytreatment = "remove"
        } else {
          tempfolder = paste("Simmap Overlap Outputs/Single Family", temprow$trait1, temprow$trait2)
          familytreatment = "only"
        }
        tempfiles = list.files(tempfolder, pattern = tempfamily)
        tempRealFile = tempfiles[which(str_detect(tempfiles, "REAL"))]
        tempDummyFile = tempfiles[which(str_detect(tempfiles, "DUMMY"))]
        RealDF = read.csv(paste0(tempfolder,"/",tempRealFile))
        DummyDF = read.csv(paste0(tempfolder,"/",tempDummyFile))
        
        # plot simmap overlap plots
        templabel = paste0(familytreatment, tempfamily)
        
        calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = templabel)
        
        # Third plot: observed state proportions
        dfDummyMelt <- melt(RealDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
        
        # Create the basic boxplot without x-axis labels (xaxt = "n") and without outlines
        boxplot(ObservedState.prop ~ ObservedState, data = dfDummyMelt, xlab = "", ylab = "Observed State Proportion", xaxt = "n", outline = FALSE, boxwex = 0.5, col = "lightgray")
        
        # Add points to the plot
        points(jitter(as.numeric(dfDummyMelt$ObservedState), amount = 0.05), dfDummyMelt$ObservedState.prop, pch = 16, col = "darkred", cex = 0.6)
        
        # Add rotated x-axis labels
        axis(1, at = 1:length(unique(dfDummyMelt$ObservedState)), 
             labels = FALSE)
        text(1:length(unique(dfDummyMelt$ObservedState)), 
             par("usr")[3] - 0.001, srt = 45, adj = 1.2, 
             labels = as.character(unique(dfDummyMelt$ObservedState)), xpd = TRUE, cex=0.8)
      } # end if SimmapOverlap
      
    } # end for j in 1:7
    
  } # end for i in 1: length(families)
  dev.off()
  
} else if (tempcount == 1) {
  pdf(file = "plots by family - SimmapOverlap BayesTraits Jackknife SingleFamily - 1test.pdf", height = 12, width = 8)
  par(mfcol=c(3,2), mai=c(1,1,0.8,0.3), oma=c(3,2,3,2),
      font.main=1, cex.main=1.25, cex.lab=1.1, cex.axis=1.1)
  families = famcounts$tempFamily[which(famcounts$n == tempcount)]
  for (i in 1:length(families)) {
    tempfamily = families[i]
    
    # plot single-family trees
    subsetout = subsettreedata(columns = c(temprow$trait1, temprow$trait2), cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladesubsetvalue = tempfamily, newdata = newdatadf, newtree = treefile)
    subsettree = subsetout$subsettree
    subsetdf = subsetout$subsetdf
    
    # no tips
    plot.phylo(subsettree, type = "f", show.tip.label = TRUE, cex = 0.2)
    title(main = paste(tempfamily))
    
    # trait1 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py = c("blue","red")
    treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait1] == 1)]
    tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=1)
    legend("bottomleft", legend = c("Noncooperative", "Cooperative"), cex = 0.9, fill=py, bty="n")
    
    # trait2 tips
    plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
    py2 = c("black","orange")
    treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,temprow$trait2] == 1)]
    tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=1)
    legend("bottomleft", legend = c("Female Song Absent", "Female Song Present"), cex = 0.9, fill=py2, bty="n")
    
    
    for (j in 1:1) {
      testcounts1 = testcounts[which(testcounts$n > 60),]
      temprow = testcounts1[j,]
      
      if (str_detect(temprow$test, "BayesTraits")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " Jackknife Outputs")
          familytreatment = "remove"
        } else {
          tempfolder = paste0("BayesTraitsDiscreteML_", temprow$trait1, "-", temprow$trait2, " SingleFamily Outputs")
          familytreatment = "only"
        }
        tempfile = list.files(tempfolder, pattern = tempfamily)
        # bayestraits plots
        dfIn = read.csv(paste0(tempfolder,"/",tempfile))
        dfInDep = dfIn[which(dfIn$Model == "Dependent"),]
        dfInInd = dfIn[which(dfIn$Model == "Independent"),]
        
        Nspecies = dfIn[1,"numSpecies"]
        
        templabel = paste0(familytreatment, tempfamily, " Nspecies=", Nspecies)
        
        tryCatch(
          {
            plotSimpleDiscreteBayes(columns = c(temprow$trait1, temprow$trait2), df = dfInDep, nocorrDdf = dfInInd, LhCol = "Lh", nsim = NULL, treelabel = NULL, newpdf = FALSE, otherlabel = "", ylabel = templabel, arrowmod = 1, roundDigits = 3)
          }, error=function(e) {
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            plot(x = 1, y = 1)
            message(paste("An Error occurred, could not plot", templabel))
          }) # end trycatch
      } # end if BayesTraits
      
      
      
      if (str_detect(temprow$test, "SimmapOverlap")) {
        if (str_detect(temprow$test, "Jackknife")) {
          tempfolder = paste("Simmap Overlap Outputs/Jackknife", temprow$trait1, temprow$trait2)
          familytreatment = "remove"
        } else {
          tempfolder = paste("Simmap Overlap Outputs/Single Family", temprow$trait1, temprow$trait2)
          familytreatment = "only"
        }
        tempfiles = list.files(tempfolder, pattern = tempfamily)
        tempRealFile = tempfiles[which(str_detect(tempfiles, "REAL"))]
        tempDummyFile = tempfiles[which(str_detect(tempfiles, "DUMMY"))]
        RealDF = read.csv(paste0(tempfolder,"/",tempRealFile))
        DummyDF = read.csv(paste0(tempfolder,"/",tempDummyFile))
        
        # plot simmap overlap plots
        templabel = paste0(familytreatment, tempfamily)
        
        calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = templabel)
        
        # Third plot: observed state proportions
        dfDummyMelt <- melt(RealDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
        
        # Create the basic boxplot without x-axis labels (xaxt = "n") and without outlines
        boxplot(ObservedState.prop ~ ObservedState, data = dfDummyMelt, xlab = "", ylab = "Observed State Proportion", xaxt = "n", outline = FALSE, boxwex = 0.5, col = "lightgray")
        
        # Add points to the plot
        points(jitter(as.numeric(dfDummyMelt$ObservedState), amount = 0.05), dfDummyMelt$ObservedState.prop, pch = 16, col = "darkred", cex = 0.6)
        
        # Add rotated x-axis labels
        axis(1, at = 1:length(unique(dfDummyMelt$ObservedState)), 
             labels = FALSE)
        text(1:length(unique(dfDummyMelt$ObservedState)), 
             par("usr")[3] - 0.001, srt = 45, adj = 1.2, 
             labels = as.character(unique(dfDummyMelt$ObservedState)), xpd = TRUE, cex=0.8)
      } # end if SimmapOverlap
      
    } # end for j in 1:1
    
  } # end for i in 1: length(families)
  dev.off()
}
