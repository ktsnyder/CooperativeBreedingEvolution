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
browniefunction(columns = c("CoopBreed", feature), newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv", newtree = "birdzillatreeMaybeConsensus.nex", nsim = 101, islog = feature, plotsimmaps = FALSE)
} 

source("plotbrownie.R")
for (i in 1:2) {
  feature <- songfeatures[i]
plotbrownie(data = paste0(Sys.Date(),"CoopBreed",feature, "brownie101sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = NULL, newpdf = TRUE, nsim = 101, islog = TRUE)
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
source("plotbrownie.R")
feature <- "Syll.song.final"
subset <- subsettreedata(columns = c("CoopBreed", "Syll.song.final"), newtree = "birdzillatreeMaybeConsensus.nex", newdata = "2021-08-18CoopSong_AnyCoopEqualsCoop_NatCommsSubset.csv")
pdf("Brownie_Syllsongfinal_fakeCoopData.pdf", width = 8, height = 12)
par(mfrow = c(4,2))
for (j in 1:20) {
  df <- subset$subsetdf
  fakecoops <- sample(1:length(df$species), 16)
  df$CoopBreed <- rep(0, times = length(df$species))
  df$CoopBreed[fakecoops] <- 1
  browniefunction(columns = c("CoopBreed", feature), newdata = df, newtree = "birdzillatreeMaybeConsensus.nex", nsim = 100, islog = feature, plotsimmaps = FALSE, otherlabel = paste0("fakeCoopData", j, "_"))

plotbrownie(data = paste0(Sys.Date(),columns[1], columns[2],"fakeCoopData",j, "_brownie",nsim,"sim.csv"), columns = c("CoopBreed",feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = "fakedata", newpdf = FALSE, nsim = 100, islog = TRUE)
}