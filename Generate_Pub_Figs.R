## Run analyses to generate figures in manuscript
## 
## 

# Table of jackknifed BayesTraits results - Jackknifed Coop/FS (8 median dependent rates + bayesfactor)
# Table of # of species with each combo Coop/FS
# Phylanova for Griesser sociality metrics - just a table

# Large multipanel figure of simmap overlap boxplots with pvals - Jackknifed Coop/FS
# 8 multipanel of simmap overlap boxplots for Each combo Coop/FS
# 8 multipanel of Bayestraits dependent outputs for Each combo Coop/FS
# Brownie for MeanCoopTie2Noncoop + song features
### Brownie for Song rep min/max + MeanCoopTie2Noncoop
# Jackknife Brownie SongRep + CoopBreed
# Brownie for Song rep + Griesser sociality metrics
# Brownie for Song rep + Other coop breed classification methods
# Scatterbox for Griesser sociality metrics+CoopBreed classification methods + song repertoire only



newdata = "2023-11-15_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"


#### Brownie ---- 
# Brownie - Cooperative Breeding and all song features, with song repertoire min/max
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Song.rep.min", "Song.rep.max", "Song.rep.final", "Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")
CBcolumn <- "MeanCoopTie2Noncoop"
discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 500
currentlabel <- "HackettOscine"

for (k in 1:length(songfeatures)) {
  print(Sys.time())
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
} # end for k


# Brownie - All other cooperative breeding classification methods with Song repertoire
CBcolumns <- c("MeanCoopTie2Coop", "AnyCoopEqualsCoop", "MeanCoopOmitTies")
discreteCatLabels = c("Non-cooperative", "Cooperative")
feature <- "Song.rep.final"
nsim = 500
currentlabel <- "HackettOscine"

for (k in 1:length(CBcolumns)) {
  print(Sys.time())
  CBcolumn = CBcolumns[k]
  print(CBcolumn)
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
} # end for k

# Brownie - Song.rep.final and other sociality metrics
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving")
discLabelList = list(c("Non-Colonial", "Colonial"), c("Two or fewer caretakers", "More than two caretakers"), c("Season or shorter social bonds", "Longest social bonds"), c("Groups Pair or Smaller", "Groups Larger than Pair"), c("One caretaker", "Two or more caretakers"), c("Asocial", "Social"), c("Smaller groups", "Largest group sizes"), c("Shortest social bonds", "Season or longer social bonds"), c("Non-familial", "Familial") )
feature <- "Song.rep.final"
nsim = 500
currentlabel <- "HackettOscine"

for (k in 5:length(SocialColumns)) {
  print(Sys.time())
  CBcolumn = SocialColumns[k]
  discreteCatLabels = discLabelList[[k]]
  print(CBcolumn)
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
} # end for k


# Jackknife by family - Brownie MeanCoopTie2Noncoop and Song.rep.final
source("jackknifingbrownie.R")

jackbrowniefunction(columns = c("MeanCoopTie2Noncoop", "Song.rep.final"), islog = TRUE, matemodel = "ARD", matensim = 100, allcsvs = TRUE, plotsimmaps = FALSE, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladeJackvalues = "Mimidae", otherlabel = "HackettOscine")
