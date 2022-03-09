## Execute analysis code
## Coded by Kate T Snyder
## Created 8/18/2021
## Last Edited 3/8/2022

setwd("~/Desktop/CooperativeBreedingEvolution")

songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate", "Continuity")
newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
currentlabel <- "PasserineTreeEricson-AllCoop"

# test ER/ARD brownie
source("findQrates")
#findQrates(columns = "CoopBreed", plot=TRUE, newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = "PasserineTreeEricson-MeanCoop")
findQrates(columns = c("CoopBreed","Syll.song.final"), plot=TRUE, newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = currentlabel)
findQrates(columns = "CoopBreed", plot=FALSE, newtree = "birdzillatreeMaybeConsensus.nex", newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
## getting these warnings: Warning messages: 1: In rstate(p/sum(p)) : Some probabilities (slightly?) < 0. Setting p < 0 to zero. - only when doing AnyCoop though? Not MeanCoop

source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
#below run with "2020-10-11ConsensusPasserineTreeHack100.nex" 8/19/21
#run with "birdzillatreeMaybeConsensus.nex" 8/24/21
songfeatures <- c("Syll.song.min", "Syll.song.max", "Syll.song.final", "Song.rep.final", "Song.rep.min", "Song.rep.max")
nsim = 500
for (k in 1:6) {
  newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
  treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
  currentlabel <- "PasserineTreeEricson-AnyCoop"
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c("CoopBreed", feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  
  plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
  
  newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
  currentlabel <- "PasserineTreeEricson-MeanCoop"
  browniefunction(columns = c("CoopBreed", feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  
  plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
  
}



source("scatterboxes.R")
scatterboxes(DiscreteTrait = "CoopBreed", newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", otherlabel = "AnyCoop")
scatterboxes(DiscreteTrait = "CoopBreed", newdata = "2022-03-08CoopSong_MeanCoop_All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", otherlabel = "MeanCoop")


source("plotACEtree.R")
for (k in 1:7) {
  feature <- songfeatures[k]
plotACEtree(columns = c("CoopBreed", feature), cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = newdata, newtree = newtree, islog = feature, discretelabels = c("Non-cooperative","Cooperative"), discretemodel = "ARD", otherlabel = currentlabel)
} 


# Change each "Cooperative" species to "Noncooperative" iteratively to test robustness of results
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
  thisdf <- read.csv("2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
  thisdf[which(thisdf$species == coopbird), "CoopBreed"] <- 0
  browniefunction(columns = c("CoopBreed", feature), newdata = thisdf, newtree = "birdzillatreeMaybeConsensus.nex", nsim = 500, islog = feature, plotsimmaps = FALSE)
  plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, "brownie500sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = paste0(coopbird,"NonCoop"), newpdf = FALSE, nsim = 500, islog = TRUE)
}
dev.off()

## Simple bayestraits discrete test for Female Song and CoopBreed
## Using btw V1 - run 3/9/2022
require("btw")
subsetbtw <- subsettreedata(columns = c("CoopBreed","FemaleSong"), newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", skinnydata = TRUE)
currentlabel <- "PasserineTreeEricson-AnyCoop"
subsetdf <- subsetbtw$subsetdf
subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- "1"
subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- "0"
subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
subsettree <- subsetbtw$subsettree
subsetdf %>% group_by(CoopBreed,FemaleSong) %>% summarise(n=n())
simplebtwOut <- set.seed(10)
nsim = 2000
for (n in 1:nsim) {
  nocorrD <- Discrete(subsettree, subsetdf)
  corrD <- Discrete(subsettree, subsetdf, dependent=TRUE)
  lrtestresults <- lrtest(corrD, nocorrD)
  tempRow <- cbind(corrD, lrtestresults)
  simplebtwOut <- rbind(simplebtwOut, tempRow)
}
simplebtwOut <- as.data.frame(simplebtwOut)
means <- apply(X = simplebtwOut,MARGIN = 2,FUN = mean)
meansdf <- as.data.frame(rbind(means,means))
#meansdf <- as.data.frame(as.matrix(means))
pvalMed <- median(simplebtwOut$pval)
pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
dev.off()

subsetbtw <- subsettreedata(columns = c("CoopBreed","FemaleSong"), newdata = "2022-03-08CoopSong_MeanCoop_All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", skinnydata = TRUE)
currentlabel <- "PasserineTreeEricson-MeanCoop"
subsetdf <- subsetbtw$subsetdf
subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- "1"
subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- "0"
subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
subsettree <- subsetbtw$subsettree
subsetdf %>% group_by(CoopBreed,FemaleSong) %>% summarise(n=n())
simplebtwOut <- set.seed(10)
nsim = 2000
for (n in 1:nsim) {
  nocorrD <- Discrete(subsettree, subsetdf)
  corrD <- Discrete(subsettree, subsetdf, dependent=TRUE)
  lrtestresults <- lrtest(corrD, nocorrD)
  tempRow <- cbind(corrD, lrtestresults)
  simplebtwOut <- rbind(simplebtwOut, tempRow)
}
simplebtwOut <- as.data.frame(simplebtwOut)
means <- apply(X = simplebtwOut,MARGIN = 2,FUN = mean)
meansdf <- as.data.frame(rbind(means,means))
#meansdf <- as.data.frame(as.matrix(means))
pvalMed <- median(simplebtwOut$pval)
pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
dev.off()


# BayesTraits song features
source("btwfunction.R")
source("BayesPlots_choosebin.R")
songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate")
newdata = "2022-03-08CoopSong_AnyCoopEqualsCoop_All.csv"
#newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
currentlabel <- "PasserineTreeEricson-AnyCoop"
nsim = 100
for (k in 1:6) { 
  feature <- songfeatures[k]
  btwfunction(MateParam = "CoopBreed",SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"BayesCoopBreed",feature, nsim, "reps.csv") 
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = "CoopBreed",SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
} 


# Do sylls/song brownie with fake CoopBreed data (randomly assign 16 species as Coop)
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