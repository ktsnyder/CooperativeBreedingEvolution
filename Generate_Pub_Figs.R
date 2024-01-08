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

#### Simmap Overlap CoopBreed/FS ----
dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = 1000, treelabel = "HackettOscine", datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = 2500, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
calcHuel(dfout4, dfDummy4)

dfout4 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ MeanCoopTie2Noncoop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 1000 HackettOscine .csv")
dfDummy4 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ MeanCoopTie2Noncoop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 2500 HackettOscine .csv")

tempdfGather = tempdfSub %>% gather("Rate", "RateValue", c(q12:q43, q12.1:q43.1)) 


#### Simmap Overlap Sociality metrics ----

socialityMetrics = c("Griesser2023.Colonial01", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LongSocialBonds", "Griesser2023.MoreThanTwoCaretakers", "Griesser2017FamilialLiving", "Final.polygyny", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes",   "Griesser2023.SeasonOrLongerSocialBonds", "MeanCoopTie2Noncoop")
treelabel = "HackettOscine"
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
#dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
#calcHuelout = calcHuel(dfout4, dfDummy4)
#require(gridExtra)
plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")
nPlots = length(calcHuelout)
#m3 <- marrangeGrob(calcHuelout, ncol = 1, nrow = nPlots)
#ggsave(plotname, m3, width = 7.5, height = 3.5*nPlots, units = "in")
#socialityPlots[[i]] <- calcHuelout
}


# plot transition counts as if BayesTraits
source("plotSimpleDiscreteBayes.R")
source("test_trait_overlap_simmaps.R")
source("transition_plot.R")
require(stringr)
plotlist = list()
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
  next
  
  calcHuelout = calcHuel(tempdfDep, tempdfInd)
  Transitions = calcHuelout$TransitionStats$logPairwisePostHoc[2]
  pvals = calcHuelout$TransitionStats$logPairwisePostHoc[7]
  pvaldf = as.data.frame(cbind(Transitions, pvals))
  pvaldf$SignificanceLabel = "n.s."
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.05)] = "p < 0.05"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.01)] = "p < 0.01"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.001)] = "p < 0.001"
  pvaldf$SignificanceLabel[which(pvaldf$p.value < 0.0001)] = "p < 0.0001"
  
  RateRef = as.data.frame(rbind(c("FS0to1inCoop0", "q12"),c("Coop0to1inFS0", "q13"),c("FS1to0inCoop0", "q21"), c("Coop0to1inFS1", "q24"),c("Coop1to0inFS0", "q31"), c("FS0to1inCoop1", "q34"), c("Coop1to0inFS1", "q42"),c("FS1to0inCoop1", "q43")))
  colnames(RateRef) <- c("Transitions", "qRate")
  ratePvals = merge(RateRef, pvaldf, by.x = "Transitions", by.y = "Transition")
  
  
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
  plottitle = paste(trait1, trait2, "nsims:", nsims)
  outlist[[i]] = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle)
  names(outlist)[i] <- paste(trait1, trait2, sep = "_")
  templist = outlist[[i]]
  plotlist[[i]] = templist$transition_plot
  
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
ggsave("transition plots sociality FS_onepage.pdf", single_page_plot, width = 18, height = 38, units = "in")
#ggsave("transition plots sociality FS.pdf", mSoc, width = 24.5, height = 7.9*(round(nTransPlots/2)), units = "in", limitsize = FALSE)





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

filelist = list.files(path = "/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/BayesTraitsDiscrete", pattern = "removed", full.names = T)
mediandf = set.seed(10)
bayesbarplots = list()
for (i in 1:length(filelist)) {
  tempfile = filelist[i]
  splitfile = str_split(tempfile, "_removed", simplify = F)
  splitfile2 = str_split(splitfile[[1]][2], "MeanCoop")
  removedFam = splitfile2[[1]][1]
  splitfile3 = str_split(splitfile2[[1]][2], " ")
  trait1 = paste0("MeanCoop", splitfile3[[1]][1])
  trait2 = splitfile3[[1]][2]
  nSims = str_remove(splitfile3[[1]][3], "sims.csv")
  tempdf = read.csv(tempfile)
  tempdf$CalcBFLh = 2*(tempdf$Lh - tempdf$Lh.1)
  tempdf$CalcBFModelLhs = 2*(tempdf$model1.Lh - tempdf$model2.Lh)
  write.csv(tempdf, paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/BayesTraitsDiscrete/", Sys.Date(), "_2023-11-18BayesTraitsDiscrete_HackettOscine_removed", removedFam, "MeanCoopTie2Noncoop FemaleSong_Agg01 100sims_WithCalcBFs.csv"), row.names = FALSE)
  
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
ggsave("jackknifed BayesTraitsDiscrete RateScatterBarplots MeanCoopTie2Noncoop FemaleSong_Agg01 HackettOscine.pdf", mBayes, width = 8, height = 10, units = "in")


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

plotlist = list()
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
ggsave("jackknifed Simmap Overlaps MeanCoopTie2Noncoop FemaleSong_Agg01 HackettOscine.pdf", m1, width = 8, height = 9, units = "in")


#### Generate Counts Table ----

dfIn = read.csv(newdata)
OscineSubset = subsettreedata(newdata = dfIn, newtree = treefile)
df = OscineSubset$subsetdf
columnsToSummarize = c("Griesser2023.Colonial01", "Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.TwoOrMoreCaretakers", "Griesser2017FamilialLiving", "Final.polygyny", "MeanCoopTie2Noncoop") 
otheraxiscolumns = c("FemaleSong_Agg01", "MeanCoopTie2Noncoop")

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
#write.csv(fulltable, "SuppTable_binary traits state intersections NumSpecies.csv", row.names = T)
