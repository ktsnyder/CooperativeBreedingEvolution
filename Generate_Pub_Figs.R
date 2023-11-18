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
### Jackknife Brownie SongRep + CoopBreed
### Brownie for Song rep + Griesser sociality metrics
### Brownie for Song rep + Other coop breed classification methods
# Scatterbox for Griesser sociality metrics+CoopBreed classification methods + song repertoire only
setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")


#newdata = "2023-11-15_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
#df = read.csv(newdata)
#df[,c("Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds")] <- sapply(df[,c("Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds")], FUN = as.integer)
#write.csv(df, file = "2023-11-17_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv", row.names = F)
newdata = "2023-11-17_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
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

for (k in 1:length(SocialColumns)) {
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

jackbrowniefunction(columns = c("MeanCoopTie2Noncoop", "Song.rep.final"), islog = TRUE, matemodel = "ARD", matensim = 100, allcsvs = TRUE, plotsimmaps = FALSE, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladeJackvalues = NULL, otherlabel = "HackettOscine")

source("plotbrowniejacks.R")
plotbrowniejacks(columns = c("MeanCoopTie2Noncoop", "Song.rep.final"), allcsvs = TRUE, islog = TRUE, otherlabel = "HackettOscine", csvFolder = "BrownieJackknifeOutputs")

#### phylANOVA and ScatterBoxes ----
source("scatterboxes.R")
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving")
discLabelList = list(c("Non-Colonial", "Colonial"), c("Two or fewer caretakers", "More than two caretakers"), c("Season or shorter social bonds", "Longest social bonds"), c("Groups Pair or Smaller", "Groups Larger than Pair"), c("One caretaker", "Two or more caretakers"), c("Asocial", "Social"), c("Smaller groups", "Largest group sizes"), c("Shortest social bonds", "Season or longer social bonds"), c("Non-familial", "Familial") )

for (i in 5:length(SocialColumns)) {
  CBcolumn = SocialColumns[i]
  tempLabels = discLabelList[[i]]
  scatterboxes(DiscreteTrait = CBcolumn, newdata = newdata, newtree = treefile, otherlabel = currentlabel, discreteCategoryLabels = tempLabels)
}

# phylANOVA for multi-group traits - Table
multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
songtrait = "Song.rep.final"
phynovaDF = set.seed(10)
for (i in 1:length(multigrouptraits)) {
tempgrouptrait = multigrouptraits[i]
subsets = subsettreedata(columns = c(tempgrouptrait, songtrait), newdata = newdata, newtree = treefile)
subsetdf = subsets$subsetdf
discvec = subsetdf[,tempgrouptrait]
names(discvec) = subsetdf$species
contvec = subsetdf[,songtrait]
contvec = log(contvec)
names(contvec) = subsetdf$species
subsettree = subsets$subsettree
phylANOVAout= phylANOVA(subsettree, x = discvec, y = contvec, nsim = 50000, posthoc = TRUE)
phylANOVAp=phylANOVAout$Pf
temprow = c(tempgrouptrait, songtrait, phylANOVAp)
phynovaDF = rbind(phynovaDF, temprow)
phynovaDF = as.data.frame(phynovaDF)
colnames(phynovaDF) <- c("DiscreteTrait", "ContinuousTrait", "PhylANOVApval")
write.csv(phynovaDF, file = "phylANOVA outputs Song.rep.final w Coops.csv", row.names = F)
print(paste(tempgrouptrait, songtrait))
print(phylANOVAout)
}


#### Jackknife CoopBreed/FS ----

# BayesTraits
source("simplebtwDiscrete.R")
source(file = "subsettreedata.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"

trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"
subset = subsettreedata(columns = c(trait1, trait2), newdata = newdata, newtree = treefile)
subsetdf = subset$subsetdf
subsettree = subset$subsettree
familyvec = unique(subsetdf$Family3_BirdtreeMatchSpecies2)
familyvec = c("None", familyvec)

for (i in 1:length(familyvec)) {
  tempfam = familyvec[i]
  jackedsongdf <- subsetdf  #have to do this so songdf doesn't get whittled down every time the for loop loops
  jackedsongdf <- jackedsongdf[which(jackedsongdf[, cladesubsetcolumn] != tempfam),]
  
  numSpecies = length(jackedsongdf[,1])
  
  print(paste("BayesTraits", trait1, trait2, tempfam, i, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  notjackedvec <- subsettree$tip.label %in% as.character(jackedsongdf$species)
  dropforjack <- which(notjackedvec == FALSE)
  
  jacktree <- drop.tip(subsettree,tip = dropforjack)
  
  fulllabel = paste0("HackettOscine_removed", tempfam)
  
simplebtwDiscrete(columns = c(trait1, trait2), newdata = jackedsongdf, newtree = jacktree, treelabel = fulllabel, nsim = 100, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, savecsvs = TRUE, KeepBTInputFiles = FALSE)
}


# Simmap Overlap
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"
columns = c(trait1, trait2)

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2[which(familycounts$n > 2)]
familyvec = c("None", familyvec)

nsims_real = 50
nsims_dummy = 200

for (j in 1:length(familyvec)) {
  
  familyToRemove = familyvec[j]
  templabel = paste0(columns[1], " ", columns[2], " ", "removed",familyToRemove)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  print(paste("Simmap Overlap", trait1, trait2, familyToRemove, j, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap")
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap")
} # end cycle through families for jackknife