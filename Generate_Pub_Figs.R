## Run analyses to generate figures in manuscript
## 
## 

# Table of jackknifed BayesTraits results - Jackknifed Coop/FS (8 median dependent rates + bayesfactor)
# Table of # of species with each combo Coop/FS
### Phylanova for Griesser sociality metrics - just a table

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
#newdata = "2023-11-17_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
#newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
#newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
subsetout =subsettreedata(newdata = newdata, newtree = treefile)
subsetdf= subsetout$subsetdf
subsetdf$AnyNoncoopEqualsNoncoop = subsetdf$MeanCoopTie2Noncoop
subsetdf$AnyNoncoopEqualsNoncoop[which(subsetdf$SourceDiscrepancy == 1)] = 0
subsetdf %>% group_by(AnyNoncoopEqualsNoncoop) %>% count

#### Brownie ---- 
# Brownie - Cooperative Breeding and all song features, with song repertoire min/max
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Song.rep.final", "Syllable.rep.final") #, "Syll.song.final", "Duration.final", "Song.rep.min", "Song.rep.max", "Syllable.rep.min", "Syllable.rep.max") # "Interval.final",
CBcolumn <- "HighConfidence_Coop"
discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 500
currentlabel <- "UpdatedSongData"

treefile = "2024-05-26ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex"
currentlabel <- "UpdatedSongData_Hackett4Passerine-MeanEdge-IgnoreAbsent"

treeIn = read.nexus(treefile)
tree = drop.tip(treeIn, tip = which(!treeIn$tip.label %in% OscineTree$tip.label))

for (k in 1:length(songfeatures)) {
  print(Sys.time())
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c(CBcolumn, feature), newdata = newdata, newtree = tree, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  
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
dfIn = read.csv(newdata)
dfIn$AnyNoncoopEqualsNoncoop = dfIn$MeanCoopTie2Noncoop
dfIn$AnyNoncoopEqualsNoncoop[which(dfIn$SourceDiscrepancy == 1)] = 0
newdata = dfIn
CBcolumns <- c("MeanCoopTie2Noncoop","MeanCoopTie2Coop", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "MeanCoopOmitTies")
discreteCatLabels = c("Non-cooperative", "Cooperative")
feature <- "Syll.song.final"
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
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0VsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving")
discLabelList = list(c("Non-Colonial", "Colonial"), c("Two or fewer caretakers", "More than two caretakers"), c("Season or shorter social bonds", "Longest social bonds"), c("Groups Pair or Smaller", "Groups Larger than Pair"), c("One caretaker", "Two or more caretakers"), c("Asocial", "Social"), c("Smaller groups", "Largest group sizes"), c("Shortest social bonds", "Season or longer social bonds"), c("Non-familial", "Familial") )
SocialColumns <- c("Griesser2023.MoreThanTwoCaretakers", "Griesser2023.TwoOrMoreCaretakers")
discLabelList = list(c("Two or fewer caretakers", "More than two caretakers"), c("One caretaker", "Two or more caretakers") )
feature <- "Song.rep.final"
feature <- "Syllable.rep.final"
nsim = 500
currentlabel <- "HackettOscineER"

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
source("test_trait_overlap_simmaps.R")
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving")
discLabelList = list(c("Non-Colonial", "Colonial"), c("Two or fewer caretakers", "More than two caretakers"), c("Season or shorter social bonds", "Longest social bonds"), c("Groups Pair or Smaller", "Groups Larger than Pair"), c("One caretaker", "Two or more caretakers"), c("Asocial", "Social"), c("Smaller groups", "Largest group sizes"), c("Shortest social bonds", "Season or longer social bonds"), c("Non-familial", "Familial") )

for (i in 1:length(SocialColumns)) {
  CBcolumn = SocialColumns[i]
  tempLabels = discLabelList[[i]]
  scatterboxes(DiscreteTrait = CBcolumn, newdata = newdata, newtree = treefile, otherlabel = currentlabel, discreteCategoryLabels = tempLabels)
}

# phylANOVA for multi-group traits - Table
#multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
multigrouptraits = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" , "HighConfidence_Coop","MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "grouping_Griesser2023", "social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "social_bonds_Griesser2023")
songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final", "Song.rep.max", "Song.rep.min","Syllable.rep.max", "Syllable.rep.min")
#songtraits = c("Song.rep.max", "Song.rep.min","Syllable.rep.max", "Syllable.rep.min", "Song.rep.final","Syllable.rep.final")
#multigrouptraits = c("HighConfidence_Coop","HighConfidence_Coop", "HighConfidence_Coop", "HighConfidence_Coop", "AnyNoncoopEqualsNoncoop", "AnyNoncoopEqualsNoncoop")
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
dfIn = read.csv(newdata)
dfIn$AnyNoncoopEqualsNoncoop = dfIn$MeanCoopTie2Noncoop
dfIn$AnyNoncoopEqualsNoncoop[which(dfIn$SourceDiscrepancy == 1)] = 0
newdata = dfIn

phynovaDF = set.seed(10)
for (j in 1:length(songtraits)) {
  songtrait = songtraits[j]
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
    Nspecies = length(subsettree$tip.label)
    Ngroups = length(unique(discvec))
    
    tryCatch({
      # Code that might fail
      phylANOVAout = phylANOVA(subsettree, x = discvec, y = contvec, nsim = 50000, posthoc = TRUE)
      phylANOVAp = phylANOVAout$Pf
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, phylANOVAp)
    }, error = function(e) {
      # Code to run in case of an error
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, NA)
      message("Error in phylANOVA computation: ", e$message)
    })
    
    phynovaDF = rbind(phynovaDF, temprow)
    phynovaDF = as.data.frame(phynovaDF)
    colnames(phynovaDF) <- c("DiscreteTrait", "DiscreteNumGroups", "ContinuousTrait", "n_Species", "PhylANOVApval")
    print(paste(tempgrouptrait, songtrait))
    print(phylANOVAout)
  }
  if (j == 5) {
    write.csv(phynovaDF, file = paste(Sys.Date(), "phylANOVA outputs FinalSongs_Coops.csv"), row.names = F)
  }
}
write.csv(phynovaDF, file = paste(Sys.Date(), "phylANOVA outputs Songs_Coops.csv"), row.names = F)


#### Simmap Overlap CoopBreed/FS ----
source("test_trait_overlap_simmaps.R")
dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = 100, treelabel = "HackettOscine", datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = 200, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
calcHuel(dfout4, dfDummy4)

dfout4 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ MeanCoopTie2Noncoop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 1000 HackettOscine .csv")
dfDummy4 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ MeanCoopTie2Noncoop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 2500 HackettOscine .csv")

tempdfGather = tempdfSub %>% gather("Rate", "RateValue", c(q12:q43, q12.1:q43.1)) 

# With new consensus tree 5/26/2024
treeIn = read.nexus("2024-05-26ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex")
newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
tree = drop.tip(treeIn, which(!treeIn$tip.label %in% OscineTree$tip.label ))
tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 0.000000000000001
sum(OscineTree$edge.length)
sum(tree$edge.length)
dfout4 <- CharacterSimmaps(columns = c("HighConfidence_Coop","FemaleSong_Agg01"), df = newdata, tree =  tree, dummy = FALSE, nsims = 500, treelabel = "HackettPasserineMeanEdgeIgnoreAbsent", datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c("HighConfidence_Coop","FemaleSong_Agg01"), df = newdata, tree =  tree, dummy = TRUE, nsims = 500, treelabel = "HackettPasserineMeanEdgeIgnoreAbsent", datalabel = NULL, dummyMethod = "makeSimmap")
calcHuelout = calcHuel(dfout4, dfDummy4, otherlabel = "HackettPasserineMeanEdgeIgnoreAbsent")


#### Simmap Overlap Sociality metrics ----

dfIn = read.csv(newdata)
dfIn$NonkinNoncoop0_FamAndOrCoop1 = NA
dfIn$NonkinNoncoop0_FamAndOrCoop1[dfIn$social_system_incl_nk_coop_Griesser2017 == "no_fam"] <- 0
dfIn$NonkinNoncoop0_FamAndOrCoop1[dfIn$social_system_incl_nk_coop_Griesser2017 %in% c("coop_families", "family", "nk-coop")] = 1

unique(dfIn$grouping)
dfIn$AsocPairLarge0_SmallGroup1 = NA
dfIn$AsocPairLarge0_SmallGroup1[which(dfIn$grouping %in% c("asocial","pair","large_groups"))] = 0
dfIn$AsocPairLarge0_SmallGroup1[which(dfIn$grouping == "small_groups")] = 1
dfIn$AsocSmallLarge0_Pair1 = NA
dfIn$AsocSmallLarge0_Pair1[which(dfIn$grouping %in% c("asocial","small_groups","large_groups"))] = 0
dfIn$AsocSmallLarge0_Pair1[which(dfIn$grouping == "pair")] = 1
newdata = dfIn
sum(!is.na(dfIn$AsocPairLarge0_SmallGroup1))
sum(!is.na(dfIn$AsocSmallLarge0_Pair1))

socialityMetrics = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.MoreThanTwoCaretakers")#, "MeanCoopTie2Noncoop", "MeanCoopTie2Coop")
socialityMetrics = c("MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
socialityMetrics = c("NonkinNoncoop0_FamAndOrCoop1", "AsocPairLarge0_SmallGroup1", "AsocSmallLarge0_Pair1")

treelabel = "HackettOscineER"
nsims_real = 500
nsims_dummy = 500

source("test_trait_overlap_simmaps.R")
socialityPlots = list()
socialityStats = list()
for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
dfout4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
calcHuelout = calcHuel(dfout4, dfDummy4)
require(gridExtra)
plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")
#nPlots = length(calcHuelout)-5
#m3 <- marrangeGrob(calcHuelout, ncol = 1, nrow = nPlots)
#ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
nPlots = 6
m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
#socialityPlots[[i]] <- calcHuelout
}

#### Simmap overlap and counts ----
# Using already-generated data from CharacterSimmaps

for (i in 1:length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(tempMetric)
  
  columns = c(tempMetric, "FemaleSong_Agg01")
  trait1 = tempMetric
  trait2 = columns[2]
  
  filelist = list.files(path = "Simmap Overlap Outputs", pattern = "overlap_counts", full.names = T)
  filelist = filelist[which(str_detect(filelist,paste(tempMetric, "FemaleSong_Agg01")))]
  filelist = filelist[which(str_detect(filelist,".csv"))]
  IndFile = filelist[which(str_detect(filelist,"DUMMY"))]
  DepFile = filelist[which(str_detect(filelist,"REAL"))]
  
  tempdfDep = read.csv(DepFile)
  tempdfInd = read.csv(IndFile)
  
  nsims_real = length(tempdfDep[,1])
    nsims_dummy = length(tempdfInd[,1])
  
  calcHuelout = calcHuel(dfout = tempdfDep, dfDummy = tempdfInd, otherlabel = "ExpObsCompare")
  require(gridExtra)
  plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts-ExpObsCompare.pdf")
  nPlots = length(calcHuelout)-4
  print(calcHuelout$mediansRow)
  m3 <- marrangeGrob(calcHuelout[1:nPlots], ncol = 1, nrow = nPlots) # nrow can be nPlots if no "filename" or "TransitionStats" in calcHuelout
  ggsave(plotname, m3, width = 7.5, height = 4.2*nPlots, units = "in")
}


# use previously-generated data to make specific panels for figure (boxplots simmap overlap - HighConf Coop and Fem Song)
DepFile = "Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine .csv"
IndFile = "Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine .csv"
tempdfDep = read.csv(DepFile)
tempdfInd = read.csv(IndFile)
calcHuelout = calcHuel(dfout = tempdfDep, dfDummy = tempdfInd)
filename1 = "simmap overlap HighConfidence_Coop FemaleSong_Agg01 ObservedStates boxplot.pdf"
calcHuelout$p3
#ggsave(filename1, calcHuelout$p3, width = 6, height = 5, units = "in", device = "pdf")

DepFile = "Simmap Overlap Outputs/ Griesser2017FamilialLiving FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine .csv"
IndFile = "Simmap Overlap Outputs/ Griesser2017FamilialLiving FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine .csv"
tempdfDep = read.csv(DepFile)
tempdfInd = read.csv(IndFile)
calcHuelout = calcHuel(dfout = tempdfDep, dfDummy = tempdfInd)
filename2 = "simmap overlap Griesser2017FamilialLiving FemaleSong_Agg01 ObservedStates boxplot.pdf"
#ggsave(filename2, calcHuelout$p3, width = 6, height = 5, units = "in", device = "pdf")

# multistate 
multiFile = "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Multistate/Sim Output CSVs/simmap overlap social_system_incl_nk_coop_Griesser2017 FemaleSong_Agg01 1500 sims.csv"
overlapdf = read.csv(multiFile)
calcHuelout3 = calcHuelflex(overlapdf) # from scratch simmap overlap processing
calcHuelout3$boxplotStates
filename3 = "simmap overlap social_system_incl_nk_coop_Griesser2017 FemaleSong_Agg01 ObservedStates boxplot.pdf"
#ggsave(filename3, calcHuelout3$boxplotStates, width = 9, height = 5, units = "in", device = "pdf")


# plot transition counts as if BayesTraits
source("plotSimpleDiscreteBayes.R")
source("test_trait_overlap_simmaps.R")
source("transition_plot.R")
require(stringr)
plotlist = list()
plotlistGray = list()
transitionList = list()
outlist = list()
for (i in 1:length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(tempMetric)
  
  columns = c(tempMetric, "FemaleSong_Agg01")
  trait1 = tempMetric
  trait2 = columns[2]
  
  filelist = list.files(path = "Simmap Overlap Outputs", pattern = "overlap_counts", full.names = T)
  filelist = filelist[which(str_detect(filelist,paste(tempMetric, "FemaleSong_Agg01")))]
  filelist = filelist[which(str_detect(filelist,".csv"))]
  IndFile = filelist[which(str_detect(filelist,"DUMMY"))]
  DepFile = filelist[which(str_detect(filelist,"REAL"))]
  
  tempdfDep = read.csv(DepFile)
  tempdfInd = read.csv(IndFile)
  nsims = nsims_real = length(tempdfDep[,1])
  nsims_dummy = length(tempdfInd[,1])
  
  transitioncols = colnames(tempdfDep)[18:25]
  ExpectedCols = paste0(transitioncols, "Expected")
  ChiSqStat = rowSums((tempdfDep[,transitioncols] - tempdfDep[,ExpectedCols])^2/tempdfDep[,ExpectedCols])
  ChiSqPvals = pchisq(ChiSqStat, df = 7)
  tempdfDep = cbind(tempdfDep, ChiSqStat, ChiSqPvals)
  
  #plot(density(tempdfDep$ChiSqPvals, na.rm = TRUE), main = tempMetric)
  hist(tempdfDep$ChiSqPvals, main = tempMetric, breaks = 20)
  abline(v = 0.05, col = "red")
  #next
  
  calcHuelout = calcHuel(tempdfDep, tempdfInd)
  
  Transitions = calcHuelout$TransitionStats$logPairwisePostHoc[2]
  pvals = calcHuelout$TransitionStats$logPairwisePostHoc[7]
  pvaldf = as.data.frame(cbind(Transitions, pvals))
  pvaldf$SignificanceLabel = "n.s."
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.05)] = "p < 0.05"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.01)] = "p < 0.01"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.001)] = "p < 0.001"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.0001)] = "p < 0.0001"
  
  NtimesActualGreaterThanExpected = calcHuelout$TransitionStats$NtimesActualGreaterThanExpected
  NtimesActualGreaterThanExpecteddf = cbind(names(NtimesActualGreaterThanExpected), NtimesActualGreaterThanExpected)
  NtimesActualGreaterThanExpecteddf = as.data.frame(NtimesActualGreaterThanExpecteddf)
  colnames(NtimesActualGreaterThanExpecteddf) = c("Transition", "CountNActualGreaterThanExpected")
  NtimesActualGreaterThanExpecteddf$CountNActualGreaterThanExpected = as.integer(NtimesActualGreaterThanExpecteddf$CountNActualGreaterThanExpected)
  pvaldf = merge(pvaldf, NtimesActualGreaterThanExpecteddf, by = "Transition")
  pvaldf$FractionActualGreaterThanExpected = pvaldf$CountNActualGreaterThanExpected/nsims
  
  RateRef = as.data.frame(rbind(c("FS0to1inCoop0", "q12"),c("Coop0to1inFS0", "q13"),c("FS1to0inCoop0", "q21"), c("Coop0to1inFS1", "q24"),c("Coop1to0inFS0", "q31"), c("FS0to1inCoop1", "q34"), c("Coop1to0inFS1", "q42"),c("FS1to0inCoop1", "q43")))
  colnames(RateRef) <- c("Transitions", "qRate")
  ratePvals = merge(RateRef, pvaldf, by.x = "Transitions", by.y = "Transition")
  
  # remove these lines if want to use ANOVA pval labels instead
  ratePvals$PercentTrendingLabel = "<80%"
  ratePvals$PercentTrendingLabel[which(ratePvals$FractionActualGreaterThanExpected > .8 | ratePvals$FractionActualGreaterThanExpected < .2)] <- ">80%"
  ratePvals$PercentTrendingLabel[which(ratePvals$FractionActualGreaterThanExpected > .9 | ratePvals$FractionActualGreaterThanExpected < .1)] <- ">90%"
  ratePvals$PercentTrendingLabel[which(ratePvals$FractionActualGreaterThanExpected > .95 | ratePvals$FractionActualGreaterThanExpected < .05)] <- ">95%"
  ratePvals$PercentTrendingLabel[which(ratePvals$FractionActualGreaterThanExpected > .99 | ratePvals$FractionActualGreaterThanExpected < .01)] <- ">99%"
  
  tempdfDep[,paste0(colnames(tempdfDep)[18:25],"DifferenceFromExpected")] = tempdfDep[,colnames(tempdfDep)[18:25]] - tempdfDep[,paste0(colnames(tempdfDep)[18:25],"Expected")]
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop0DifferenceFromExpected")] <- "q12"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS0DifferenceFromExpected")] <- "q13"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop0DifferenceFromExpected")] <- "q21"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS1DifferenceFromExpected")] <- "q24"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS0DifferenceFromExpected")] <- "q31"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop1DifferenceFromExpected")] <- "q34"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS1DifferenceFromExpected")] <- "q42"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop1DifferenceFromExpected")] <- "q43"
  
  labx0 = paste("Female Song Absent")
  labx1 = paste("Female Song Present")
  
  if (trait1 == "Final.polygyny") {
    lab0x = paste("Monogamy")
    lab1x = paste("Polygyny")
  } else if (grepl("coop", trait1, ignore.case = T)) {  #(str_detect(trait1, "coop")) {
    lab0x = paste("Non-Cooperative")
    lab1x = paste("Cooperative")
  } else if (str_detect(trait1, "Kin")) {
    lab0x = paste("Non-kin")
    lab1x = paste("Kin")
  } else if (str_detect(trait1, "Familial")) {
    lab0x = paste("Non-Familial Living")
    lab1x = paste("Familial Living")
  } else if (str_detect(trait1, "Colonial")) {
    lab0x = paste("Non-Colonial")
    lab1x = paste("Colonial") 
  } else if (str_detect(trait1, "GroupsLargerThanPair")) {
    lab0x = paste("Asocial or pair")
    lab1x = paste("Small or large groups") 
  } else if (str_detect(trait1, "LongSocialBonds")) {
    lab0x = paste("Bonds last one season or less")
    lab1x = paste("Multi-year bonds") 
  } else if (str_detect(trait1, "MoreThanTwoCaretakers")) {
    lab0x = paste("Two or fewer caretakers")
    lab1x = paste("More than two caretakers") 
  } else if (str_detect(trait1, "TwoOrMoreCaretakers")) {
    lab0x = paste("Fewer than two caretakers")
    lab1x = paste("Two or more caretakers") 
  } else if (str_detect(trait1, "Asocial")) {
    lab0x = paste("Asocial")
    lab1x = paste("Pair or group sociality") 
  } else if (str_detect(trait1,"SeasonOrLonger")) {
    lab0x = paste("Bonds last less than one season")
    lab1x = paste("Season or longer social bonds") 
  } else if (str_detect(trait1, "LargestGroupSizes")) {
    lab0x = paste("Asocial, pair, or small groups")
    lab1x = paste("Large groups") 
  }
  
  #plotSimpleDiscreteBayes(columns = columns, df = tempdfDep, nocorrDdf = NULL, LhCol = "Lh", nsim = nsims, treelabel = "HackettOscine", newpdf = TRUE, cladesubsetvalue = NULL, ylabel = "Transition Counts", arrowmod = 0.5, otherlabel = "TransitionCountArrows_halfTotalTransIndependent", roundDigits = 3)
  trait1StateLabels = c(lab0x, lab1x)
  trait2StateLabels = c(labx0, labx1)
  plottitle = paste(trait1, trait2, "nsims:", nsims_real)
  outlist[[i]] = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle, center="median")
  names(outlist)[i] <- paste(trait1, trait2, sep = "_")
  templist = outlist[[i]]
  plotlist[[i]] = templist$transition_plot
  plotlistGray[[i]] = templist$transitionplot_GrayNS
  
  # calcHuelout = calcHuel(tempdfDep, tempdfInd)
  # require(gridExtra)
  # plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCountsStatsSig.pdf")
  # nPlots = length(grep("^p", names(calcHuelout)))
  # m3 <- marrangeGrob(calcHuelout[1:nPlots], ncol = 1, nrow = nPlots)
  # #ggsave(plotname, m3, width = 7.5, height = 3.9*nPlots, units = "in")
  # socialityPlots[[i]] <- calcHuelout
  # socialityStats[[i]] <- calcHuelout$TransitionStats
  # names(socialityStats)[i] <- calcHuelout$filename
  # 
  # print(calcHuelout$filename)
  # print(calcHuelout$TransitionStats)
}
#plotlist[[i]]$transition_plot
#plotlist[[i]]$transition_df

require(gridExtra)
require(grid)
nTransPlots= length(plotlist)
#mSoc <- marrangeGrob(plotlist[1:nTransPlots], ncol = 2, nrow = 2)
grobs_with_margins <- lapply(plotlist, function(plot) {
  plot_with_margin <- plot + 
    theme(plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) # Adjust margins as needed
  ggplotGrob(plot_with_margin)
})
# Arrange the grobs on a single page
single_page_plot <- grid.arrange(grobs = grobs_with_margins, ncol = 2, nrow = 6) # Adjust ncol and nrow as needed

# Save the arranged plot to a file
ggsave("transition plots sociality FS_onepage_median_percentStates_percentTrending.pdf", single_page_plot, width = 18, height = 38, units = "in")

# With gray nonsig arrows
nTransPlots= length(plotlistGray)
#mSoc <- marrangeGrob(plotlist[1:nTransPlots], ncol = 2, nrow = 2)
grobs_with_marginsGray <- lapply(plotlistGray, function(plot) {
  plot_with_margin <- plot + 
    theme(plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) # Adjust margins as needed
  ggplotGrob(plot_with_margin)
})
# Arrange the grobs on a single page
single_page_plotGray <- grid.arrange(grobs = grobs_with_marginsGray, ncol = 2, nrow = 6) # Adjust ncol and nrow as needed

# Save the arranged plot to a file
ggsave("transition plots sociality FS_onepage_median_percentStates_weightsNumTransitions_GrayNonsig.pdf", single_page_plotGray, width = 18, height = 38, units = "in")

# just Coop-FS plot
ggsave("transition plot Tie2Noncoop FSAgg median_percentStates_weightsNumTransitions_labsPercentTrending_GrayNonsig.png", plotlistGray[[11]], width = 9, height = 5, units = "in", device = "png")

# just Coop-FS plot
ggsave("transition plot HighConf_Coop FSAgg median_percentStates_weightsNumTransitions_labsPercentTrending_GrayNonsig.pdf", plotlistGray[[1]], width = 9, height = 5, units = "in", device = "pdf")

# just target 4 plots
grobs_with_marginsGray4 <- lapply(plotlistGray[c(3,5,6,10)], function(plot) {
  plot_with_margin <- plot + 
    theme(plot.margin = margin(t = 40, r = 30, b = 10, l = 20, unit = "pt")) # Adjust margins as needed
  ggplotGrob(plot_with_margin)
})
single_page_plotGray4 <- grid.arrange(grobs = grobs_with_marginsGray4, ncol = 2, nrow = 2)
ggsave("transition plots 4TargetSociality FS_median_percentStates_weightsNumTransitions_labsPercentTrending_GrayNonsig.png", single_page_plotGray4, width = 20, height = 13, units = "in", device = "png")



#### Jackknife BayesTraits CoopBreed/FS ----

# BayesTraits
source("simplebtwDiscrete.R")
source(file = "subsettreedata.R")
.BayesTraitsPath = "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"

trait1 = "MeanCoopTie2Noncoop"
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
subset = subsettreedata(columns = c(trait1, trait2), newdata = newdata, newtree = treefile)
subsetdf = subset$subsetdf
subsettree = subset$subsettree
familyvec = unique(subsetdf$Family3_BirdtreeMatchSpecies2)
cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2"
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

filelist = list.files(path = "/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/BayesTraitsDiscrete", pattern = "removed", full.names = T)
mediandf = set.seed(10)
bayesbarplots = list()
for (i in 1:length(filelist)) {
  tempfile = filelist[i]
  splitfile = str_split(tempfile, "_removed", simplify = F)
  splitfile2 = str_split(splitfile[[1]][2], "Coop")
  removedFam = splitfile2[[1]][1]
  splitfile3 = str_split(splitfile2[[1]][2], " ")
  trait1 = paste0("Coop", splitfile3[[1]][1])
  trait2 = splitfile3[[1]][2]
  nSims = str_remove(splitfile3[[1]][3], "sims.csv")
  tempdf = read.csv(tempfile)
  tempdf$CalcBFLh = 2*(tempdf$Lh - tempdf$Lh.1)
  tempdf$CalcBFModelLhs = 2*(tempdf$model1.Lh - tempdf$model2.Lh)
  write.csv(tempdf, paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/BayesTraitsDiscrete/", Sys.Date(), "_BayesTraitsDiscrete_HackettOscine_removed", removedFam, " Coop FemaleSong_Agg01 100sims.csv"), row.names = FALSE)
  
  tempdfSub = tempdf[which(tempdf$CalcBFLh > 2),]
  nSig = length(tempdfSub[,1])
  
  tempdfGather = tempdfSub %>% gather("Rate", "RateValue", c(q12:q43, q12.1:q43.1)) 
  tempdfGather$Model = NA
  tempdfGather$Model[which(tempdfGather$Rate %in% c("q12", "q13", "q21", "q31", "q24", "q42", "q34", "q43"))] <- "Dependent"
  tempdfGather$Model[which(tempdfGather$Rate %in% c("q12.1", "q13.1", "q21.1", "q31.1", "q24.1", "q42.1", "q34.1", "q43.1"))] <- "Independent"
  bayesbarplots[[i]] <- ggplot(tempdfGather, aes(x = Rate, y = RateValue, fill = Model)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    geom_point(position = position_jitter(width = 0.2), size = 0.7) +
    theme_minimal() +
    labs(y = "Rate", x = paste("removed", removedFam, "nSignificant:", nSig), fill = "Model") +
    scale_fill_manual(values = c("Independent" = "blue", "Dependent" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10)) +
    ggtitle(paste(trait1, trait2))
  
  medianvec = sapply(tempdfSub, FUN=median, simplify = T)
  meanvec = sapply(tempdfSub, FUN=mean, simplify = T)
  temprow = c(medianvec, meanvec, removedFam, trait1, trait2, nSims, nSig)
  names(temprow) <- c(paste0(names(medianvec), ".median"), paste0(names(meanvec), ".mean"), "removedFam", "trait1", "trait2", "nSims", "nSignificant")
  mediandf = rbind(mediandf, temprow)
  mediandf = as.data.frame(mediandf)
}
mediandf$X = NULL
colnames(mediandf)
write.csv(mediandf, paste("jackknifed BayesTraitsDiscrete MeanMedianResults SigOnly", trait1, trait2, ".csv"), row.names = FALSE)

# Plot jackknife BayesTraitsDiscrete output rates as boxplots
require(gridExtra)
mBayes <- marrangeGrob(bayesbarplots, ncol = 2, nrow = 4)
ggsave("jackknifed BayesTraitsDiscrete RateScatterBarplots Coop FemaleSong_Agg01 HackettOscine.pdf", mBayes, width = 8, height = 10, units = "in")


#### Jackknife Coop/FS Simmap Overlap ----
trait1 = "MeanCoopTie2Noncoop"
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
columns = c(trait1, trait2)

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 2)]
familyvec = c("None", familyvec)

nsims_real = 50
nsims_dummy = 200

plotlist = list()
for (j in 1:length(familyvec)) {
  
  familyToRemove = familyvec[j]
  templabel = paste0(columns[1], " ", columns[2], " ", "removed",familyToRemove)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2_AVONET != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  print(paste("Simmap Overlap", trait1, trait2, familyToRemove, j, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap")
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap")
  
  HuelOut = calcHuel(dfout = dfout4, dfDummy = dfDummy4, newplot = FALSE, otherlabel = templabel, plot_ggplots_pdf = FALSE)
  plotlist[[j]] <- HuelOut
  plot(length(plotlist))
} # end cycle through families for jackknife

# Plot jackknife output
require(gridExtra)
index = 0
plotlist2 = list()
for (i in 1:54) {
  tempplots = plotlist[[i]]
  for (k in 1:3) {
    index = index+1
    plotlist2[[index]] = tempplots[[k]]
  }
}
m1 <- marrangeGrob(plotlist2, ncol = 1, nrow = 3)
ggsave(paste(Sys.Date(), "jackknifed Simmap Overlaps Coop FemaleSong_Agg01 HackettOscine.pdf"), m1, width = 8, height = 9, units = "in")



#### Jackknife CoopBreed/Song features Brownie ----
source("jackknifingbrownie.R")
songfeatures = c("Song.rep.final", "Syllable.rep.final")
currentlabel = "UpdatedSongData"

for (i in 1:2) {
  tempfeature = songfeatures[i]
  print(tempfeature)
  subset = subsettreedata(columns = c("HighConfidence_Coop", tempfeature), newtree = treefile, newdata = newdata)
  familycounts = subset$subsetdf %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
  familysubset = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 1)]
jackbrowniefunction(columns = c("HighConfidence_Coop", tempfeature), islog = T, matemodel = "ARD", matensim = 100, allcsvs = T, plotsimmaps = F, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2_AVONET", otherlabel = currentlabel, cladeJackvalues = familysubset)
}

#### Generate Counts Table ----

dfIn = read.csv(newdata)
OscineSubset = subsettreedata(newdata = dfIn, newtree = treefile)
df = OscineSubset$subsetdf
columnsToSummarize = c("Griesser2023.Colonial01", "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.TwoOrMoreCaretakers", "Griesser2017FamilialLiving", "Final.polygyny", "MeanCoopTie2Noncoop") 
otheraxiscolumns = c("FemaleSong_Agg01", "HighConfidence_Coop") # "MeanCoopTie2Noncoop")

fulltable = c(0,1)
for (i in 1:length(columnsToSummarize)) {
  tablesegment = c(0,1)
  for (j in 1: length(otheraxiscolumns)) {
  tempCol = columnsToSummarize[i]
  othertempCol = otheraxiscolumns[j]
  
  if (tempCol != othertempCol) {
grouped_df <- df %>%
  group_by(.data[[othertempCol]], .data[[tempCol]]) %>%
  summarize(Count = n(), .groups='drop')
filtered_df = grouped_df %>% filter(!is.na(.data[[othertempCol]]) & !is.na(.data[[tempCol]]))

wide_df <- filtered_df %>%
  spread(key = .data[[tempCol]], value = Count)
wide_df[is.na(wide_df)] <- 0
wide_df = as.data.frame(wide_df)
  } else {
    wide_df = as.data.frame(matrix(data = NA, nrow = 2, ncol = 2))
  }
rownames(wide_df) = paste(othertempCol, c(0,1), sep = "_")
wide_df = wide_df[,which(colnames(wide_df) != othertempCol)]
colnames(wide_df) = paste(tempCol, c(0,1), sep = "_")
tablesegment = rbind(tablesegment, wide_df)
  }
  tablesegment = tablesegment[which(rownames(tablesegment) != 1),]
  fulltable = cbind(fulltable, tablesegment)
}
fulltable = fulltable[,which(colnames(fulltable) != "fulltable")]
fulltable
#write.csv(fulltable, "SuppTable_binary traits state intersections NumSpecies - HighConfidence_Coop.csv", row.names = T)


# Multistate counts table
dfIn = read.csv(newdata)
OscineSubset = subsettreedata(columns = c("social_system_incl_nk_coop_Griesser2017", "FemaleSong_Agg01"), newdata = dfIn, newtree = treefile)
df = OscineSubset$subsetdf
MultistateCBxFS = df %>% group_by(social_system_incl_nk_coop_Griesser2017, FemaleSong_Agg01) %>% count
write.csv(MultistateCBxFS, file = "SuppTable_multistate CB-Fam FSintersections NumSpecies.csv", row.names = FALSE)


# male/female-biased helping
bigdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-10-26_Aggregate_CBSource_Data_AllColumns.csv")
bigdataRiehl = bigdata[, c("BirdtreeSpecies", "Dispersal_Riehl", "Frequency", "O.F", "Kin", "Social_breeding_system_when_cooperative", "Category")]
bigdataRiehl$HelperSexBias = NA
bigdataRiehl$HelperSexBias[which(str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "female", negate = T) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "male", negate = T))] <- "NotSpecified"
bigdataRiehl$HelperSexBias[which(str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "female", negate = T) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "male", negate = F))] <- "MaleOnly"
bigdataRiehl$HelperSexBias[which(str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "female", negate = F) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, " male", negate = F))] <- "BothMaleAndFemale"
bigdataRiehl$HelperSexBias[which(str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "female", negate = F) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, " males", negate = F))] <- "BothMaleAndFemale"
bigdataRiehl$HelperSexBias[which(str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, "female", negate = F) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, " male", negate = T) & str_detect(bigdataRiehl$Social_breeding_system_when_cooperative, " males", negate = T))] <- "FemaleOnly"
View(bigdataRiehl[which(!is.na(bigdataRiehl$Social_breeding_system_when_cooperative)),])

df = read.csv(newdata)
df2 = merge(df, bigdataRiehl, by.x = "species", by.y = "BirdtreeSpecies", all.x = T, suffixes = c("","_Riehl"))
df2 %>% group_by(FemaleSong_Agg01, HelperSexBias) %>% count

df = read.csv(newdata)
length(which(df$SourceDiscrepancy == 1 & !is.na(df$FemaleSong_Agg01)))

