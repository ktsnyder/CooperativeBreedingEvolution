## Temp - copied from Run_CoopBreed_Analyses for running 3/10/2022-3/10/2022
# 12:30pm - re-ran CoopBreed_species_summary and 
# 
setwd("~/Desktop/CooperativeBreedingEvolution/")
setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")
require(dplyr)

source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Syll.song.final", "Song.rep.final", "Syll.rep.final")
classmethods = c("MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop")
newdata = "2022-03-10CoopSong_All.csv"
olddata = "2022-03-09CoopSong__All_preJetzCorrection.csv"
treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" 
treelabel = "PasserineTreeEricson-"

dfnew <- read.csv(newdata)


jackspecies <- dfnew$species[which(dfnew$numSourcesNonCoop == dfnew$numSourcesCoop & dfnew$numSourcesNonCoop != 0 & !is.na(dfnew$MeanCoopTie2Noncoop & !is.na(dfnew$FemaleSong)))][1:4]

dfnew %>% group_by(MeanCoopTie2Noncoop ,Griesser2017FamilialLiving) %>% summarize(n = n())
dfnew %>% group_by(MeanCoopTie2Noncoop, FemaleSong) %>% summarize(n = n())



dotTree(tree,cbind(x,y),fsize=0.7,colors=colors)



# plot BayesTraits song rep with more bins
BToutdf <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-15BayesMeanCoopTie2NoncoopSong.rep.final20reps.csv")
source("BayesPlots_choosebin.R")
transitionBinplots(MateParam = "MeanCoopTie2Noncoop", SongParam = "Song.rep.final", df = BToutdf, newpdf = TRUE, nsim = 20, binnum = 5)
transitionBinplots(MateParam = "MeanCoopTie2Noncoop", SongParam = "Song.rep.final", df = BToutdf, newpdf = TRUE, nsim = 20, binnum = 1, sigonly = TRUE)


# ## ACE tree
 treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex"
# treefile <- "2020-10-11ConsensusPasserineTreeHack100.nex"
 source("plotACEtree.R")
 for (feature in songfeatures) {  # repeated Syllrep, Songrep, Interval with "ER"
  currentclassmethod <- classmethods[2]
  currentlabel <- paste0("PasserTreeEric-ER-",currentclassmethod)
  print(feature)
  print(currentlabel)
plotACEtree(columns = c(currentclassmethod,feature), newdata = newdata, newtree = treefile, islog = feature, discretelabels = c("NonCoop","Coop"), discretemodel = "ER", otherlabel = currentlabel)

  currentclassmethod <- classmethods[3]
  currentlabel <- paste0("PasserTreeEric-ER-",currentclassmethod)
  print(feature)
  print(currentlabel)
  plotACEtree(columns = c(currentclassmethod,feature), newdata = newdata, newtree = treefile, islog = feature, discretelabels = c("NonCoop","Coop"), discretemodel = "ER", otherlabel = currentlabel)
}
# 

newdata = "2022-03-10CoopSong_All.csv"
treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" 
currentclassmethod = classmethods[2]
currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
print(currentlabel)
nsim = 11
source("btwfunction.R")
source("BayesPlots_choosebin.R")
for (k in 1:3) {
  feature <- songfeatures[k]
  print(feature)
  print(currentlabel)
  btwfunction(MateParam = currentclassmethod,SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
  #feature <- songfeatures[k]
  filename <- paste0(Sys.Date(),"Bayes", currentclassmethod,feature, nsim, "reps.csv")
  BTdf <- read.csv(filename)
  #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
  colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
  transitionBinplots(MateParam = currentclassmethod,SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
}



## Brownie jacks by family

source("findQrates.R")
source("plotbrowniejacks.R")
source("jackknifingbrownie.R")

features <- c("Syll.song.final", "Song.rep.final") #,"Syllable.rep.final")
currentclassmethod <- "MeanCoopTie2Noncoop"
for (i in 1:length(features)) {
 # i=2
  feature = features[i]
  print(feature)
  nsim = 100
  columns = c(currentclassmethod,feature)
  #jackspecies <- dfnew$species[which(dfnew$numSourcesNonCoop == dfnew$numSourcesCoop & dfnew$numSourcesNonCoop != 0 & !is.na(dfnew[,feature]))]  # nsim = 100
  #jackspecies <- dfnew$species[which(!is.na(dfnew[,currentclassmethod]) & !is.na(dfnew[,feature]))]  # nsim = 20
  
  browniejackout <- jackbrowniefunction(columns = columns, islog = feature, matemodel = "ARD", matensim = nsim, allcsvs = TRUE, plotsimmaps = FALSE, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family.x", cladeJackvalues = NULL, otherlabel = NULL)
  
  browniejacks = browniejackout$brownielist
  familyvec = browniejackout$familyvec
  plotbrowniejacks(MateParam = currentclassmethod, SongParam = feature, browniejacks = browniejacks, allcsvs = "2022-03-17", nsim = nsim)
}


# Merge in Webb female song data
df <- read.csv(newdata)
WebbFSdf <- read.csv("Webb et al 2016 Female Song Plumage Data.csv")
dfnew <- merge(df, WebbFSdf, by.x = "species", by.y = "TipLabel")
dfnew %>% group_by(MeanCoopTie2Noncoop, Female_song_score) %>% summarise(n=n())


require("btw")
currentclassmethod = "MeanCoopTie2Noncoop"
songcol = "Female_song_score"
songcol = "FemaleSong"
dfnew <- dfnew[which(dfnew$Female_song_score %in% c("Present","Absent")),]
subsetbtw <- subsettreedata(columns = c(currentclassmethod,songcol), newdata = dfnew, newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", skinnydata = TRUE)
currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
#subsetbtw <- subsettreedata(columns = c(currentclassmethod,songcol), newdata = dfnew, newtree = "2022-03-16ConsensusPasserineTreeHackett4_1000.nex", skinnydata = TRUE)
#currentlabel <- paste0("PasserTreeHack-",currentclassmethod)
#subsetdfOdom <- subsetbtw$subsetdf
#subsettreeOdom <- subsetbtw$subsettree
colnames(subsetdfOdom)[which(colnames(subsetdfOdom) == currentclassmethod)] <- "CoopBreed"
subsetdfOdom[,"FemaleSong"][which(subsetdfOdom[,"FemaleSong"] == "Present")] <- "1"
subsetdfOdom[,"FemaleSong"][which(subsetdfOdom[,"FemaleSong"] == "Absent")] <- "0"
subsetdfOdom$CoopBreed <- as.character(subsetdfOdom$CoopBreed)
#subsetdfWebb <- subsetbtw$subsetdf
#subsettreeWebb <- subsetbtw$subsettree
colnames(subsetdfWebb)[which(colnames(subsetdfWebb) == currentclassmethod)] <- "CoopBreed"
subsetdfWebb[,"Female_song_score"][which(subsetdfWebb[,"Female_song_score"] == "Present")] <- "1"
subsetdfWebb[,"Female_song_score"][which(subsetdfWebb[,"Female_song_score"] == "Absent")] <- "0"
subsetdfWebb$CoopBreed <- as.character(subsetdfWebb$CoopBreed)

subsetdf <- subsetbtw$subsetdf
subsettree <- subsetbtw$subsettree
colnames(subsetdf)[which(colnames(subsetdf) == currentclassmethod)] <- "CoopBreed"
subsetdf[,"Female_song_score"][which(subsetdf[,"Female_song_score"] == "Present")] <- "1"
subsetdf[,"Female_song_score"][which(subsetdf[,"Female_song_score"] == "Absent")] <- "0"
subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)

subsetdf %>% group_by(CoopBreed,Female_song_score) %>% summarise(n=n())

simplebtwOut <- set.seed(10)
nsim = 50
for (n in 1:nsim) {
  if (n %in% c(1, 2, 5,10,25,50,100,150,200, 500, 1000, 1500, 2000)) {
    print(n)
    print(Sys.time())
  }
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
pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSongWebb ", nsim, "sims.pdf"))
plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSongWebb, \nnsims =",nsim, "median pval =", pvalMed))
dev.off()


### Brownie with Hackett tree
### 
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
treefile <- "2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
currentclassmethod = "MeanCoopTie2Noncoop"
currentlabel <- paste0("PasserTreeHack-",currentclassmethod)
feature <- "Syll.song.final"
nsim = 100
newdata <- "2022-03-10CoopSong_All.csv"
print(currentlabel)
browniefunction(columns = c(currentclassmethod, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)

plotbrownie(data = paste0(Sys.Date(),currentclassmethod,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(currentclassmethod,feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)




# Attempt to use btw with BayesTraitsV4
.BayesTraitsPath <- "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"
