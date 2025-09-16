## Run analyses to generate figures in manuscript 
## Kate T Snyder
## Last edited 6/19/2024

# Source code files ----
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")
source("jackknifingbrownie.R")
source("brownie relative rates.R")
source("BrownieMultistate.R")
source("test_trait_overlap_simmaps.R")
source("find transition counts by state for 2 Discrete traits.R")
source("transition_plot.R")


# Set data, tree, and variable vectors ----
newdata = "Data_R.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex" # to perform test using the alternative consensus tree for any given analysis, replace this with "ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex"
# treefile can also be set to any individual tree extracted from a BirdTree multiphylo object in order to test across many individual trees
songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")
socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "grouping_Griesser2023", "social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "social_bonds_Griesser2023")
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"

# Supplemental Table 3 - phylANOVA ----
phynovaDF = set.seed(10)
for (j in 1:length(songtraits)) {
  songtrait = songtraits[j]
    tempgrouptrait = trait1
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
      phylANOVAout = phylANOVA(subsettree, x = discvec, y = contvec, nsim = 500, posthoc = TRUE)
      phylANOVAp = phylANOVAout$Pf
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, phylANOVAp)
    }, error = function(e) {
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, NA)
      message("Error in phylANOVA computation: ", e$message)
    })
    
    phynovaDF = rbind(phynovaDF, temprow)
    phynovaDF = as.data.frame(phynovaDF)
    colnames(phynovaDF) <- c("DiscreteTrait", "DiscreteNumGroups", "ContinuousTrait", "n_Species", "PhylANOVApval")
    print(paste(tempgrouptrait, songtrait))
    print(phylANOVAout)
}
write.csv(phynovaDF, file = paste(Sys.Date(), "phylANOVA outputs Songs.csv"), row.names = F)






# Figure 1, Supplemental Figure 1, Supplemental Table 4 - Brownie ----
source("brownie relative rates.R")
discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 50
currentlabel = "test"
CBcolumn = "HighConfidence_Coop"

for (k in 1:length(songtraits)) {
  print(Sys.time())
  feature <- songtraits[k]
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
# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = "OutputFiles", otherlabel = currentlabel)


# Supplemental Table 5 - Brownie Jackknife ----
require(dplyr)
source("jackknifingbrownie.R")
source("brownie relative rates.R")

currentlabel = "test jackknife"
subset = subsettreedata(columns = c("HighConfidence_Coop", "Song.rep.final"), newtree = treefile, newdata = newdata)
familycounts = subset$subsetdf %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familysubset = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 20)]

jackbrowniefunction(columns = c("HighConfidence_Coop", "Song.rep.final"), islog = T, matemodel = "ARD", matensim = 10, allcsvs = T, plotsimmaps = F, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2_AVONET", otherlabel = currentlabel, cladeJackvalues = familysubset)

# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = "BrownieJackknifeOutputs", otherlabel = currentlabel)


# Supplemental Figures 2 & 3; Supplemental Tables 6 & 7 - Brownie with multi-state categorical traits ----
BrownieMultistate(DiscreteTrait = "grouping_Griesser2023", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = 5, plotsimmaps = F, plotResults = T, otherlabel = "test")
BrownieMultistate(DiscreteTrait = "social_system_incl_nk_coop_Griesser2017", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = 10, plotsimmaps = F, plotResults = T, otherlabel = "test")




# Figure 3A & B; Supplemental Table 9 - Co-occurrance of cooperative breeding/familial living and female song  ----
source("test_trait_overlap_simmaps.R")
nsims_real = 10
nsims_dummy = 10

targetMetrics = c("HighConfidence_Coop", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "JetzCoopInclCockburn", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop")
temptrait2 = "FemaleSong_Agg01" 
for (i in 1: length(targetMetrics)) {
  tempMetric = targetMetrics[i]
  print(i)
  print(tempMetric)
  dfout <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
  dfDummy <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  calcHuelout = calcHuel(dfout, dfDummy)
  require(gridExtra)
  plotname = file.path("Simmap Overlap Outputs", paste(tempMetric, temptrait2, nsims_real, nsims_dummy, treelabel, "withTransCounts.pdf"))
  
  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
}





# Tables 1 & 2 - Co-occurrance of cooperative breeding (or female song) with sociality traits ----
source("test_trait_overlap_simmaps.R")
nsims_real = 20
nsims_dummy = 20

temptrait2 = "HighConfidence_Coop"
# temptrait2 = "FemaleSong_Agg01" # uncomment to run analyses found in Table 2
for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  dfout <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
  dfDummy <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  calcHuelout = calcHuel(dfout, dfDummy)
  require(gridExtra)
  plotname = file.path("Simmap Overlap Outputs", paste(tempMetric, temptrait2, nsims_real, nsims_dummy, treelabel, "withTransCounts.pdf"))
  
  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
}




# Figures 4 & 5; Supplemental Figures 5 & 7 - Transitions between discrete state combinations ---- 
# Uses files generated during steps associated with Figure 3 and Table 2
source("test_trait_overlap_simmaps.R")
source("transition_plot.R")
require(stringr)
plotlist = list()
plotlistGray = list()
transitionList = list()
outlist = list()
socialityMetrics = c("HighConfidence_Coop", "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop")
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
  
  print(DepFile)
  print(IndFile)
  tempdfDep = read.csv(DepFile)
  tempdfInd = read.csv(IndFile)
  nsims = nsims_real = length(tempdfDep[,1])
  nsims_dummy = length(tempdfInd[,1])
  
  transitioncols = colnames(tempdfDep)[18:25]
  ExpectedCols = paste0(transitioncols, "Expected")
  ChiSqStat = rowSums((tempdfDep[,transitioncols] - tempdfDep[,ExpectedCols])^2/tempdfDep[,ExpectedCols])
  ChiSqPvals = pchisq(ChiSqStat, df = 7)
  tempdfDep = cbind(tempdfDep, ChiSqStat, ChiSqPvals)
  
  hist(tempdfDep$ChiSqPvals, main = tempMetric, breaks = 20)
  abline(v = 0.05, col = "red")
  
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
  
  trait1StateLabels = c(lab0x, lab1x)
  trait2StateLabels = c(labx0, labx1)
  plottitle = paste(trait1, trait2, "nsims:", nsims_real)
  outlist[[i]] = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle, center="median")
  names(outlist)[i] <- paste(trait1, trait2, sep = "_")
  templist = outlist[[i]]
  plotlist[[i]] = templist$transition_plot
  plotlistGray[[i]] = templist$transitionplot_GrayNS
}

require(gridExtra)
require(grid)
nTransPlots= length(plotlistGray)

grobs_with_margins <- lapply(plotlistGray, function(plot) {
  plot_with_margin <- plot + 
    theme(plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) # Adjust margins as needed
  ggplotGrob(plot_with_margin)
})
# Arrange the grobs on a single page
single_page_plot <- grid.arrange(grobs = grobs_with_margins, ncol = 2, nrow = 3) # Adjust ncol and nrow as needed

# Save the arranged plot to a file
ggsave("transition plots sociality FS_onepage_median_percentStates_percentTrending.pdf", single_page_plot, width = 18, height = 38, units = "in")




# Supplemental Table 8 - Co-occurrance jackknife  ----
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
templabel = "test"
columns = c(trait1, trait2)

source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = treefile)
qrates= Qout$qrates

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 65)]
familyvec = c("None", familyvec)

nsims_real = 5
nsims_dummy = 10

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
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = qrates)
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = qrates)
  
  HuelOut = calcHuel(dfout = dfout4, dfDummy = dfDummy4, newplot = FALSE, otherlabel = templabel)
  plotlist[[j]] <- HuelOut
} # end cycle through families for jackknife

# Plot jackknife output
require(gridExtra)
index = 0
plotlist2 = list()
for (i in 1:length(plotlist)) {
  tempplots = plotlist[[i]]
  for (k in 1:3) {
    index = index+1
    plotlist2[[index]] = tempplots[[k]]
  }
}
m1 <- marrangeGrob(plotlist2, ncol = 1, nrow = 3)
ggsave(paste(Sys.Date(), "jackknifed Simmap Overlaps Coop FemaleSong_Agg01.pdf"), m1, width = 8, height = 9, units = "in")


# Supplemental Table 12 - ARD vs ER rates for binary traits ----
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0VsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "Final.polygyny", "HighConfidence_Coop", "FemaleSong_Agg01", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop")
allQout = set.seed(10)

for ( i in 1:length(SocialColumns)) {
  temptrait = SocialColumns[i]
  
  Qoutput = findQrates(columns = temptrait, newtree = treefile, newdata = newdata, plot = F)
  Qoutput
  
  browniedata = set.seed(10)
  
  browniedata$trait = temptrait
  browniedata$ERsimmapQ = gsub("ERrates ", "", Qoutput$ERrates)
  browniedata$ARDsimmapQ0to1 = gsub("ARDrates ", "", Qoutput$ARDrates[2])
  browniedata$ARDsimmapQ1to0 = gsub("ARDrates ", "", Qoutput$ARDrates[1])
  browniedata$ERsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[1]
  browniedata$ARDsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[2]
  browniedata$ARDvERsimmapQ.LRtestPval = Qoutput$anovaERARD$`Pr(>|Chi|)`[2]
  
  browniedata = as.data.frame(browniedata)
  allQout = rbind(allQout, browniedata)
}
allQout$ARDvERsimmapQ.LRtestPval.abbr = as.numeric(allQout$ARDvERsimmapQ.LRtestPval)
allQout$ARDvERsimmapQ.LRtestPval.abbr = round(allQout$ARDvERsimmapQ.LRtestPval.abbr, digits = 3)
allQout$ARDvERsimmapQ.LRtestPval.abbr[which(allQout$ARDvERsimmapQ.LRtestPval.abbr < 0.001)] <- "<0.001"
write.csv(allQout,"binary trait Qrates_rounded.csv")



# Figure 3C; Supplemental Table 11 - Co-occurrence between female song and multistate traits ----
multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "grouping_Griesser2023")
for (i in 1:length(multistateTraits)) {
  trait1 = multistateTrait = multistateTraits[i]
  trait2 = othertrait = "FemaleSong_Agg01"
  columns = c(multistateTrait, othertrait)
  subsetout = subsettreedata(columns = multistateTrait, newdata = newdata, newtree = treefile)
  subsetDisctree = subsetout$subsettree
  subsetDiscdf = subsetout$subsetdf
  
  discretetraitvecDisc = subsetDiscdf[,multistateTrait]
  names(discretetraitvecDisc) = subsetDiscdf$species
  
  ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
  
  # get rates from ace() output
  aceARDratesVec = ARDmodel$rates
  aceARDrates = cbind(1:length(aceARDratesVec), aceARDratesVec)
  aceARDrates = as.data.frame(aceARDrates)
  colnames(aceARDrates) <- c("rate_index", "rates")
  rate_index_matrix = ARDmodel$index.matrix
  groupnames = colnames(ARDmodel$lik.anc)
  nGroups = length(groupnames)
  rate_matrix = matrix(rep(NA,nGroups^2), nrow = nGroups)
  rownames(rate_matrix) = colnames(rate_matrix) = groupnames
  
  for (i in 1:nGroups) {
    for (j in 1:nGroups) {
      index = rate_index_matrix[i,j]
      if (!is.na(index)) {
        rate_matrix[i,j] = aceARDrates$rates[which(aceARDrates$rate_index == index)]
      }
    }
  }
  diagvals = rowSums(rate_matrix, na.rm = T)*-1
  diag(rate_matrix) <- diagvals
  print(rate_matrix)
  
  source("findQrates.R")
  FSrates <- findQrates(columns = othertrait, newdata = newdata, newtree = treefile)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  
  # Make simmaps from data subsetted to those with both
  subsetSong = subsettreedata(columns = c(multistateTrait, othertrait), newdata = newdata, newtree = treefile)
  subsetdf = subsetSong$subsetdf
  subsettree = subsetSong$subsettree
  
  discretetraitvec = subsetdf[,multistateTrait]
  names(discretetraitvec) = subsetdf$species
  othertraitvec = subsetdf[,othertrait]
  names(othertraitvec) = subsetdf$species
  
  simmapMultistate = make.simmap(subsettree, discretetraitvec, nsim = nsims, Q= rate_matrix, type = "discrete") 
  realDiscreteTraitVecList = list(discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec)
  
  simmapTrait2 = make.simmap(subsettree, othertraitvec, nsim = nsims, Q= FSQ, type = "discrete")
  realTrait2vecList = list(othertraitvec,othertraitvec,othertraitvec,othertraitvec,othertraitvec)
  
  ## DUMMY
  
  CoopsimtreesRand <- list()
  RanddiscretetraitvecList <- list()
  for (j in 1:nsims) { 
    Coopvec <- subsetdf[,columns[1]]
    CoopvecRandom <- sample(Coopvec)
    names(CoopvecRandom) <- subsetdf$species
    Coopsimtree <- make.simmap(tree = subsettree, x = CoopvecRandom, model = "ARD", nsim = 1, Q = rate_matrix)
    CoopsimtreesRand[[j]] <- Coopsimtree
    RanddiscretetraitvecList[[j]] <- CoopvecRandom
    if (j == 1) {
      CoopsimtreesMulti = Coopsimtree
    } else {
      CoopsimtreesMulti = c(CoopsimtreesMulti, Coopsimtree)
    }
    print(paste(j, Sys.time()))
  } # end for j in 1:nsims (Coop)
  Coopsimtrees <- CoopsimtreesRand
  
  # Make randomized versions of FemaleSong simmaps / DUMMY data
  FSsimtreesRand <- list()
  RandTrait2vecList <- list()
  print(paste("starting Dummy FemSong simmaps", Sys.time()))
  for (j in 1:nsims) {
    FSvec <- subsetdf[,columns[2]]
    FSvecRandom <- sample(FSvec)
    names(FSvecRandom) <- subsetdf$species
    FSsimtree <- make.simmap(tree = subsettree, x = FSvecRandom, model = "ARD", nsim = 1, Q = FSQ)
    FSsimtreesRand[[j]] <- FSsimtree
    RandTrait2vecList[[j]] <- FSvecRandom
    if (j == 1) {
      FSsimtreesMulti = FSsimtree
    } else {
      FSsimtreesMulti = c(FSsimtreesMulti, FSsimtree)
    }
    print(j)
    
  } # end for j in 1:nsims (FS)
  FSsimtrees<- FSsimtreesRand
  
  overlapdf = set.seed(10)
  for (k in 1:nsims) {
    # calculate overlap - real
    realOverlap = Map.Overlap(simmapMultistate[[k]], simmapTrait2[[k]])
    overlapVec = as.vector(realOverlap)
    new_names <- outer(rownames(realOverlap), colnames(realOverlap), paste, sep = "_FS")
    new_names <- as.vector(new_names)
    names(overlapVec) = paste0(new_names, "_REAL")
    overlapVec
    
    # calculate overlap - dummy
    dummyOverlap = Map.Overlap(Coopsimtrees[[k]], FSsimtrees[[k]])
    overlapVecDummy = as.vector(dummyOverlap)
    new_names2 <- outer(rownames(dummyOverlap), colnames(dummyOverlap), paste, sep = "_FS")
    new_names2 <- as.vector(new_names2)
    names(overlapVecDummy) = paste0(new_names2, "_DUMMY")
    overlapVecDummy
    
    temprow = c(k, multistateTrait, othertrait, overlapVec, overlapVecDummy)
    
    overlapdf = rbind(overlapdf, temprow)
    overlapdf = as.data.frame(overlapdf)
    colnames(overlapdf) = c("tree", "trait1", "trait2", names(overlapVec), names(overlapVecDummy))
  }
  write.csv(overlapdf, file = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.csv"))
  
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.pdf")
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:3], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 12, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = paste("simmap overlap FractionDummyLessThanMedianReal", multistateTrait, othertrait, nsims,"sims.csv"))
}
