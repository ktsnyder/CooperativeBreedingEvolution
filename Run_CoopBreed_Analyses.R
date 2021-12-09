## Execute analysis code
## Coded by Kate T Snyder
## Created 8/18/2021

setwd("~/Desktop/CooperativeBreedingEvolution")

songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate", "Continuity")
newdata = "~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2021-08-18CoopSong_MeanCoop_NatCommsSubset.csv"

source("subsettreedata.R")
source("browniefunction.R")
#below run with "2020-10-11ConsensusPasserineTreeHack100.nex" 8/19/21
#run with "birdzillatreeMaybeConsensus.nex" 8/24/21
songfeatures <- c("Syll.song.min", "Syll.song.max")
for (k in 1:2) {
feature <- songfeatures[k]
print(feature)
browniefunction(columns = c("CoopBreed", feature), newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", newtree = "birdzillatreeMaybeConsensus.nex", nsim = 101, islog = feature, plotsimmaps = TRUE)
} 

source("plotbrownie.R")
for (i in 1:2) {
  feature <- songfeatures[i]
plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, "_brownie101sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = NULL, newpdf = TRUE, nsim = 101, islog = TRUE)
}


source("scatterboxes.R")
scatterboxes(DiscreteTrait = "CoopBreed", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", newtree = "birdzillatreeMaybeConsensus.nex", otherlabel = "_AnyCoopEqualsCoop_birdzillaMaybeConsensus")



source("plotACEtree.R")
for (k in 1:7) {
  feature <- songfeatures[k]
plotACEtree(columns = c("CoopBreed", feature), cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", newtree = "birdzillatreeMaybeConsensus.nex", islog = feature, discretelabels = c("Non-cooperative","Cooperative"), discretemodel = "ARD", otherlabel = "birdzillatreeMaybeConsensus")
} 

source("btwfunction")
source("BayesPlots_choosebin.R")
for (k in 1:2) { 
  feature <- songfeatures[k]
btwfunction(MateParam = "CoopBreed",SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = 100, newtreefile = "birdzillatreeMaybeConsensus.nex", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
}
for (k in 1:2) { # need to run the rest of the features
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"BayesCoopBreed",feature, "100reps.csv")
BTdf <- read.csv(filename)
transitionBinplots(MateParam = "CoopBreed",SongParam = feature, df = BTdf,newpdf = TRUE, nsim = 100)
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