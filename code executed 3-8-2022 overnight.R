## Temp - copied from Run_CoopBreed_Analyses for running overnight 3/8/2022-3/9/2022

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
  print(feature)
  print(currentlabel)
  btwfunction(MateParam = "CoopBreed",SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"BayesCoopBreed",feature, nsim, "reps.csv") 
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = "CoopBreed",SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
} 


newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
currentlabel <- "PasserineTreeEricson-MeanCoop"
nsim = 150
for (k in 1:6) { 
  feature <- songfeatures[k]
  print(feature)
  print(currentlabel)
  btwfunction(MateParam = "CoopBreed",SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
  feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"BayesCoopBreed",feature, nsim, "reps.csv") 
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = "CoopBreed",SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
} 