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
newdata = "2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv"
dataNoSongless = read.csv(newdata)
dataNoSongless$X = NULL
#treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
currentlabel <- "Hackett-Tie2Noncoop-FSHighConf"
CBcolumn = "MeanCoopTie2Noncoop"
#columns = c("HighConfidence_FemaleSong", CBcolumn)
columns = c("HighConfidence_FemaleSong", CBcolumn)


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
nsim = 500
for (k in 1:17) {
  print(Sys.time())
  newdata = newdata
  treefile = treefile
  currentlabel <- currentlabel
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c(CBcolumn, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  
  if (file.exists(paste0("OutputFiles/",Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"))) {
    print("file exists")
    plotbrownie(data = paste0(Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
  } else if (file.exists(paste0("OutputFiles/",Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"))) {
    plotbrownie(data = paste0(Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
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
columns = c("MeanCoopTie2Noncoop", "HighConfidence_FemaleSong")
currentlabel <- c("Tie2Noncoop FSHighConf Hackett")
subsetout <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile)
subsetdf = subsetout$subsetdf
subsettree = subsetout$subsettree
findQrates(columns = columns[1], plot=TRUE, newtree = subsettree, newdata = subsetdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
findQrates(columns = columns[2], plot=TRUE, newtree = subsettree, newdata = subsetdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
plotACEtree(columns = "MeanCoopTie2Noncoop", cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = subsetdf, newtree = subsettree, discretelabels = c("Non-cooperative","Cooperative"), discretemodel = "ARD", otherlabel = currentlabel)
plotACEtree(columns = "HighConfidence_FemaleSong", cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = subsetdf, newtree = subsettree, discretelabels = c("Female Song Absent","Female Song Present"), discretemodel = "ARD", otherlabel = currentlabel)



#### Simple bayestraits discrete tests for Female Song and CoopBreed ----
source("btwDiscreteKTS.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4"
.BayesTraitsPath = "~/Documents/BayesTraitsV3"
dataIn = read.csv(newdata)
dataIn$X = NULL
#columns = c(CBcolumn, "FemaleSong_Agg01")
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
for (i in seeds[118:length(seeds)]) {
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
songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate")
#newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
#newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
#treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
currentlabel <- "PasserineTreeHackett-Tie2Noncoop"
nsim = 100
for (k in 3) { 
  feature <- songfeatures[k]
  btwfunction(columns = c(CBcolumn, feature), plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = dataNoSongless)
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"_Bayes_Tie2Noncoop_",feature, nsim, "reps.csv") 
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  #colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = "CoopBreed",SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
} 


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