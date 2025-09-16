
length(socialityPlots)
overlaphists = list()
for (i in 1:length(socialityPlots)) {
  sociality1 = socialityPlots[[i]]
  sociality1filename = sociality1$filename
  traits = str_remove(sociality1filename, "Simmap Overlap Counts ")
  traits = str_remove(traits, " 500 500 ")
  traits = str_split(traits," ")
  trait1 = traits[[1]][1]
  trait2 = traits[[1]][2]
  plotname = paste0("Simmap Overlap Outputs/", sociality1filename, ".pdf")
  nPlots = 6
  calcHuelout = sociality1[1:6]
  #m3 <- marrangeGrob(calcHuelout, ncol = 1, nrow = nPlots)
  #ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
  
  overlaphists[[3*(i-1)+1]] = sociality1$p1
  overlaphists[[3*(i-1)+2]] = sociality1$p2
  overlaphists[[3*(i-1)+3]] = sociality1$p3
  print(3*(i-1)+1)
  print(3*(i-1)+2)
  print(3*(i-1)+3)
}
histsArranged = marrangeGrob(overlaphists, ncol = length(socialityPlots), nrow = 3)
plotname = paste0("Simmap Overlap Outputs/sociality Coop - Simmap Overlap Hists.pdf")
ggsave(plotname, histsArranged, height = 7.5, width = 4.8*length(socialityPlots), units = "in")


#### transition plots ----
source("test_trait_overlap_simmaps.R")
source("transition_plot.R")
for (i in 1:length(socialityPlots)) {
  calcHuelout = socialityPlots[[i]]
  sociality1filename = calcHuelout$filename
  traits = str_remove(sociality1filename, "Simmap Overlap Counts ")
  traits = str_remove(traits, " 500 500 ")
  traits = str_split(traits," ")
  trait1 = traits[[1]][1]
  trait2 = traits[[1]][2]
  nsims = 500
  
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

# calculate stats on the difference from expected - 4/30/2024
SummarizeCountsColumns = c(paste0(colnames(tempdfDep)[18:25],"DifferenceFromExpected"), colnames(tempdfDep)[18:25], paste0(colnames(tempdfDep)[18:25],"Expected"))
RateCountSummary <- tempdfDep %>%
  summarise(across(all_of(SummarizeCountsColumns), 
                   list(mean = ~mean(.x, na.rm = TRUE),
                        p2_5 = ~quantile(.x, probs = 0.025, na.rm = TRUE),
                        p97_5 = ~quantile(.x, probs = 0.975, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE),
                        ci_lower = ~mean(.x, na.rm = TRUE) - 1.96 * (sd(.x, na.rm = TRUE) / sqrt(n())),
                        ci_upper = ~mean(.x, na.rm = TRUE) + 1.96 * (sd(.x, na.rm = TRUE) / sqrt(n())),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(mean|p2_5|p97_5|sd|ci_lower|ci_upper|min|max|median)"
  )
write.csv(RateCountSummary, paste0("Simmap Overlap Counts Summary Stats ", tempdfDep$column1[1], " ", tempdfDep$column2[1], " ", nsims, "sims.csv"))

# plot histgrams with CI marks
DifferenceFromExpectedcolumns <- c("Coop0to1inFS0DifferenceFromExpected", "Coop1to0inFS0DifferenceFromExpected", 
             "Coop0to1inFS1DifferenceFromExpected", "Coop1to0inFS1DifferenceFromExpected", 
             "FS0to1inCoop0DifferenceFromExpected", "FS1to0inCoop0DifferenceFromExpected", 
             "FS0to1inCoop1DifferenceFromExpected", "FS1to0inCoop1DifferenceFromExpected")

# Reshape the dataframe to long format
longDifferenceFromExpected_df <- tempdfDep %>%
  pivot_longer(cols = all_of(DifferenceFromExpectedcolumns),
               names_to = "variable",
               values_to = "value")

longDifferenceFromExpected_df <- left_join(longDifferenceFromExpected_df, RateCountSummary, by = "variable")

# Plotting histograms with a vertical line at x = 0
ggplot(longDifferenceFromExpected_df, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Histograms of Differences from Expected",
       x = "Difference from Expected",
       y = "Frequency")



# prep for transition_plots.R
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop0DifferenceFromExpected")] <- "q12"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS0DifferenceFromExpected")] <- "q13"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop0DifferenceFromExpected")] <- "q21"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS1DifferenceFromExpected")] <- "q24"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS0DifferenceFromExpected")] <- "q31"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop1DifferenceFromExpected")] <- "q34"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS1DifferenceFromExpected")] <- "q42"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop1DifferenceFromExpected")] <- "q43"

labx0 = paste("Non-Cooperative")
labx1 = paste("Cooperative")

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
outlist[[i]] = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle, center="median")
names(outlist)[i] <- paste(trait1, trait2, sep = "_")
templist = outlist[[i]]
plotlist[[i]] = templist$transition_plot
plotlistGray[[i]] = templist$transitionplot_GrayNS
}

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
ggsave("transition plots sociality Coop_onepage_median_percentStates_weightsNumTransitions_GrayNonsig.pdf", single_page_plotGray, width = 18, height = 38, units = "in")


## make tables of simmap overlap AND make state boxplots and transition plots
source("transition_plot.R")
socTraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.MoreThanTwoCaretakers", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop", "AnyNoncoopEqualsNoncoop", "HighConfidence_Coop") # , "Griesser2023.TwoOrMoreCaretakers" - removed because too few species in lower group
# transition plots not in main text: socTraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.Colonial01", "Final.polygyny", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
socTraits = c("HighConfidence_Coop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "MeanCoopTie2Noncoop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "JetzCoop", "DaleCoop", "DowningCoop", "BiagoliniCoop")
otherTraits = c("FemaleSong_Agg01")

# dfout1 = read.csv("Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettPasserineMeanEdgeIgnoreAbsent .csv")
# dfDummy1 = read.csv("Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettPasserineMeanEdgeIgnoreAbsent .csv")
# calcHuelout = calcHuel(dfout1, dfDummy1, otherlabel = "HackettPasserineMeanEdgeIgnoreAbsent")

dfDummy1 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/HighConfidence_Coop FemaleSong_Agg01 200trees 20simsPerTree_multitree_setQtreeCorrected_All DUMMY.csv")
dfout1 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/HighConfidence_Coop FemaleSong_Agg01 200trees 20simsPerTree_multitree_setQtreeCorrected_All REAL.csv")
calcHuelout = calcHuel(dfout1, dfDummy1, otherlabel = "Multitree")


outlist = list()
boxplotlist = list()
transplotlist = list()
outdf = set.seed(10)
allcsvfiles = list.files("Simmap Overlap Outputs", pattern = ".csv")  
for (i in 1:length(socTraits)) {
  for (j in 1:length(otherTraits)) {
  temptrait1 = socTraits[i]
  temptrait2 = otherTraits[j]
  
  
  filelist = allcsvfiles[which(str_detect(allcsvfiles, temptrait1) & str_detect(allcsvfiles, temptrait2))]
  
  if (length(filelist > 0)) {
    dummyfile = filelist[which(str_detect(filelist, "DUMMY"))]
    realfile = filelist[which(str_detect(filelist, "DUMMY", negate = T))]
    if (length(dummyfile) == 1 & length(realfile) == 1) {
      dfDummy1 = read.csv(paste0("Simmap Overlap Outputs/",dummyfile))
      dfout1 = read.csv(paste0("Simmap Overlap Outputs/",realfile))
    } else {
      dummyfile = dummyfile[which(str_detect(dummyfile, "00") & str_detect(dummyfile, "counts") & str_detect(dummyfile, "removed", negate = T))]
      realfile = realfile[which(str_detect(realfile, "00") & str_detect(realfile, "counts") & str_detect(realfile, "removed", negate = T))]
      print(dummyfile)
      print(realfile)
      
      if (length(dummyfile) > 1 | length(realfile) > 1) {
        print("skipped")
        next
      } else {
        dfDummy1 = read.csv(paste0("Simmap Overlap Outputs/",dummyfile))
        dfout1 = read.csv(paste0("Simmap Overlap Outputs/",realfile))
      }
    }
  } else {
    print(paste("no files found", temptrait1, temptrait2))
    next
  }
  
  nsims = length(dfDummy1[,1])
  calcHuelout = calcHuel(dfout = dfout1, dfDummy = dfDummy1)
  outlist[[length(outlist)+1]] <- calcHuelout
  boxplotlist[[length(boxplotlist)+1]] <- calcHuelout$p3
  
  ObsProp0Absent_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp0Absent <= median(dfout1$ObsProp0Absent))/length(dfDummy1$treenum)
  ObsProp0Present_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp0Present <= median(dfout1$ObsProp0Present))/length(dfDummy1$treenum)
  ObsProp1Absent_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp1Absent <= median(dfout1$ObsProp1Absent))/length(dfDummy1$treenum)
  ObsProp1Present_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp1Present <= median(dfout1$ObsProp1Present))/length(dfDummy1$treenum)
  
  temprow = cbind(calcHuelout$mediansRow, ObsProp0Absent_FractionDummyLessThanMedianReal, ObsProp0Present_FractionDummyLessThanMedianReal, ObsProp1Absent_FractionDummyLessThanMedianReal, ObsProp1Present_FractionDummyLessThanMedianReal)
  outdf = rbind(outdf, temprow)
  
  # Transition plots
  trait1 = temptrait1
  trait2 = temptrait2
  
  tempdfDep = dfout1
  
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
  colnames(tempdfDep)
  
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop0DifferenceFromExpected")] <- "q12"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS0DifferenceFromExpected")] <- "q13"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop0DifferenceFromExpected")] <- "q21"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS1DifferenceFromExpected")] <- "q24"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS0DifferenceFromExpected")] <- "q31"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop1DifferenceFromExpected")] <- "q34"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS1DifferenceFromExpected")] <- "q42"
  colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop1DifferenceFromExpected")] <- "q43"
  
  meanObservedRateCounts = sapply(tempdfDep[,colnames(tempdfDep)[18:25]], mean) # added 4/26/2024
  meanExpectedRateCounts = sapply(tempdfDep[,paste0(colnames(tempdfDep)[18:25],"Expected")], mean)
  meanRateCounts = cbind(names(meanObservedRateCounts), meanObservedRateCounts, meanExpectedRateCounts)
  colnames(meanRateCounts) = c("Transitions", "meanObservedRateCounts", "meanExpectedRateCounts")
  meanRateCounts = as.data.frame(meanRateCounts)
  ratePvals = merge(ratePvals, meanRateCounts)
  
  
  if (trait2 == "MeanCoopTie2Noncoop") {
    labx0 = paste("Non-Cooperative")
    labx1 = paste("Cooperative")
  } else if (trait2 == "FemaleSong_Agg01") {
    labx0 = paste("Female Song Absent")
    labx1 = paste("Female Song Present")
  } else if (trait2 == "HighConfidence_Coop") {
    labx0 = paste("Non-Cooperative")
    labx1 = paste("Cooperative")
  }
  
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

  print("ratePvals:")
  print(ratePvals)
  write.csv(ratePvals, paste0("simmap overlap counts output - ratePvals - ", trait1, trait2, nsims, ".csv"), row.names = F) # added 4/26/2024
  
  trait1StateLabels = c(lab0x, lab1x)
  trait2StateLabels = c(labx0, labx1)
  plottitle = paste(trait1, trait2, "nsims:", nsims)
  transplotout = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle, center="median")
  transplotlist[[length(transplotlist) + 1]] = transplotout$transitionplot_GrayNS
  }
}
#write.csv(outdf, file = "AllCoops simmap overlap summary table_StateCompare. .csv", row.names = FALSE)
nplotrows = ceiling(length(boxplotlist)/2)
single_page_plotBox <- grid.arrange(grobs = boxplotlist, ncol = 2, nrow = nplotrows)
ggsave("simmap overlap all boxplots_onepage_500dummySims_testFactorLabs.pdf", single_page_plotBox, width = 14, height = 56, units = "in", limitsize = FALSE)

nplotrows = ceiling(length(transplotlist)/2)
single_page_plotTrans <- grid.arrange(grobs = transplotlist, ncol = 2, nrow = nplotrows)
ggsave("simmap overlap counts All transition plots_onepage.pdf", single_page_plotTrans, width = 18, height = 45, units = "in")

length(transplotlist)
length(boxplotlist)
bothlist = list()
for (i in 1:20) {
  bothlist[[length(bothlist)+1]] = boxplotlist[[i]]
  bothlist[[length(bothlist)+1]] = transplotlist[[i]]
}

single_page_plotTransBox = grid.arrange(grobs = bothlist, ncol = 4, nrow = 10)
ggsave("simmap overlap all boxplots-transitionplots_onepage_500dummySims_testFactorLabs.pdf", single_page_plotTransBox, width = 30, height = 50, units = "in", limitsize = FALSE)



#### multistate ----
multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "grouping_Griesser2023")
#df =  read.csv("2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv")
df =  read.csv("2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK = df$social_system_incl_nk_coop_Griesser2017
df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK[which(is.na(df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK) & df$MeanCoopTie2Noncoop == 1 & df$Kin_NK == "Kin")] = "coop_families"
df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK[which(is.na(df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK) & df$MeanCoopTie2Noncoop == 1 & df$Kin_NK == "NonKin")] = "nk-coop"
df$Griesser2017SocialSystem_w_Tie2Noncoop_KinNK[which(df$species %in% c("Prunella_modularis", "Parus_caeruleus"))] = NA # two species with contradicting Griesser 2017 data
# make species counts of these coop/fam classifications
SocCounts= df %>% group_by(social_system_incl_nk_coop_Griesser2017, social_system_Griesser2017, Griesser2017SocialSystem_w_Tie2Noncoop_KinNK, Griesser2017FamilialLiving, FemaleSong_Agg01) %>% count
write.csv(SocCounts, file = "Intersecting species counts - FemaleSong and Griesser multistates CoopFam.csv", row.names = FALSE)

df$Griesser2017SocialSystem_w_Tie2Coop_KinNK = df$social_system_incl_nk_coop_Griesser2017
df$Griesser2017SocialSystem_w_Tie2Coop_KinNK[which(is.na(df$Griesser2017SocialSystem_w_Tie2Coop_KinNK) & df$MeanCoopTie2Coop == 1 & df$Kin_NK == "Kin")] = "coop_families"
df$Griesser2017SocialSystem_w_Tie2Coop_KinNK[which(is.na(df$Griesser2017SocialSystem_w_Tie2Coop_KinNK) & df$MeanCoopTie2Coop == 1 & df$Kin_NK == "NonKin")] = "nk-coop"
df$Griesser2017SocialSystem_w_Tie2Coop_KinNK[which(df$species %in% c("Prunella_modularis", "Parus_caeruleus"))] = NA # two species with contradicting Griesser 2017 data

df$social_system_incl_nk_coop_Griesser2017_noBlTiAlpAcc = df$social_system_incl_nk_coop_Griesser2017
df$social_system_incl_nk_coop_Griesser2017_noBlTiAlpAcc[which(df$species %in% c("Prunella_modularis", "Parus_caeruleus"))] = NA # two species with contradicting Griesser 2017 data
multistateTraits = "social_system_incl_nk_coop_Griesser2017_noBlTiAlpAcc"

multistateTraits = "Griesser2017SocialSystem_w_Tie2Noncoop_KinNK"
multistateTraits = "Griesser2017SocialSystem_w_Tie2Coop_KinNK"
multistateTraits = "grouping"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
nsims = 500

for (i in 1:length(multistateTraits)) {
  trait1 = multistateTrait = multistateTraits[i]
  trait2 = othertrait = "MeanCoopTie2Noncoop" #"FemaleSong_Agg01"
  columns = c(multistateTrait, othertrait)
  subsetout = subsettreedata(columns = multistateTrait, newdata = df, newtree = treefile)
  subsetDisctree = subsetout$subsettree
  subsetDiscdf = subsetout$subsetdf
  #subsetDiscdf %>% group_by(social_system_incl_nk_coop_Griesser2017) %>% count
  
  discretetraitvecDisc = subsetDiscdf[,multistateTrait]
  names(discretetraitvecDisc) = subsetDiscdf$species
  
  #simmapER <- make.simmap(subsetDisctree,discretetraitvecDisc,nsim=1,model = "ER") 
  
  #ERmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ER")
  ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
  #anovaERARD <- anova(ERmodel,ARDmodel)
  #SYMmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "SYM")
  #anovaERSYM = anova(ERmodel,SYMmodel)
  #anovaSYMARD <- anova(SYMmodel,ARDmodel)
  
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
  FSrates <- findQrates(columns = othertrait, newdata = df, newtree = treefile)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  
  # get rates from make.simmap
  #simmapARD = make.simmap(subsettree,discretetraitvec, nsim = 1, model = "ARD")
  #makesimmapARDrates = simmapARD$Q
  
  # Make simmaps from data subsetted to those with song reps
  subsetSong = subsettreedata(columns = c(multistateTrait, othertrait), newdata = df, newtree = treefile)
  subsetdf = subsetSong$subsetdf
  subsettree = subsetSong$subsettree
  
  discretetraitvec = subsetdf[,multistateTrait]
  names(discretetraitvec) = subsetdf$species
  othertraitvec = subsetdf[,othertrait]
  names(othertraitvec) = subsetdf$species
  
  simmapMultistate = make.simmap(subsettree, discretetraitvec, nsim = nsims, Q= rate_matrix, type = "discrete") 
  realDiscreteTraitVecList = list(discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec)
  plot_my_simmaps(thistrait = trait1, subsettrait = trait2, simmappy = simmapMultistate, otherlabel = "REAL", discretetraitvecList = realDiscreteTraitVecList)
  write.simmap(simmapMultistate, file = paste("multistateTrait simmaps", multistateTrait, othertrait, nsims, "sims REAL.nex"), format = "nexus", version = 1.5)
  
  simmapTrait2 = make.simmap(subsettree, othertraitvec, nsim = nsims, Q= FSQ, type = "discrete")
  realTrait2vecList = list(othertraitvec,othertraitvec,othertraitvec,othertraitvec,othertraitvec)
  
  plot_my_simmaps(thistrait = trait2, subsettrait = trait1, simmappy = simmapTrait2, otherlabel = "REAL", discretetraitvecList = realTrait2vecList)
  write.simmap(simmapTrait2, file = paste("FemaleSong simmaps", multistateTrait, othertrait, nsims, "sims REAL.nex"), format = "nexus", version = 1.5)
  
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
  
  plot_my_simmaps(thistrait = trait1, subsettrait = trait2, simmappy = CoopsimtreesMulti, otherlabel = "DUMMY", discretetraitvecList = RanddiscretetraitvecList)
  write.simmap(CoopsimtreesMulti, file = paste("multistateTrait simmaps", multistateTrait, othertrait, nsims, "sims DUMMY.nex"), format = "nexus", version = 1.5)
  
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
    #write.simmap(FSsimtree, file = paste("FemaleSong simmaps", multistateTrait, othertrait, "500sims DUMMY.nex"), format = "nexus", version = 1.5, append = TRUE)
    print(j)
    
  } # end for j in 1:nsims (FS)
  FSsimtrees<- FSsimtreesRand
  plot_my_simmaps(thistrait = trait2, subsettrait = trait1, simmappy = FSsimtreesMulti, otherlabel = "DUMMY", discretetraitvecList = RandTrait2vecList)
  write.simmap(FSsimtreesMulti, file = paste("FemaleSong simmaps", multistateTrait, othertrait, nsims, "sims DUMMY.nex"), format = "nexus", version = 1.5)
  
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


## Using prior run data
overlapdf = read.csv("simmap overlap social_system_incl_nk_coop_Griesser2017 FemaleSong_Agg01 500sims.csv")
overlapdf = read.csv("simmap overlap social_system_Griesser2017 FemaleSong_Agg01 500sims.csv")
overlapdf = read.csv("simmap overlap grouping FemaleSong_Agg01 500sims.csv")

files = c("simmap overlap grouping FemaleSong_Agg01 1000 sims.csv", "simmap overlap social_system_Griesser2017 FemaleSong_Agg01 1000 sims.csv", "simmap overlap social_system_incl_nk_coop_Griesser2017 FemaleSong_Agg01 1000 sims.csv")

for (i in 1:length(files)) {
  tempfile = files[i]
  overlapdf = read.csv(tempfile)
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = gsub(".csv", ".pdf", tempfile)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout, ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 12, units = "in", limitsize = FALSE)
}

overlapdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/simmap overlap social_system_incl_nk_coop_Griesser2017 FemaleSong_Agg01 500sims.csv")
calcHuelout = calcHuelflex(overlapdf)
#

calcHuelflex = function(overlapdf) { # currently not good bc it requires
  require(dplyr)
overlapdf$X=NULL
colnames(overlapdf)

nsims = length(overlapdf[,1])

overlapdf[,4:length(colnames(overlapdf))] = apply(overlapdf[,4:length(colnames(overlapdf))], MARGIN = 2, FUN = as.numeric) 
overlaplonger = overlapdf %>%   pivot_longer(
  cols = !c(tree, trait1, trait2),  
  names_to = "state",  
  values_to = "proportion"          # The name of the new column for the values
)
overlaplonger$Which = NA
overlaplonger$Which[which(str_detect(overlaplonger$state, "DUMMY"))] = "Dummy"
overlaplonger$Which[which(str_detect(overlaplonger$state, "REAL"))] = "Real"
overlaplonger$state <- gsub( "_DUMMY", "", overlaplonger$state)
overlaplonger$state <- gsub( "_REAL", "", overlaplonger$state)

df_long <- overlaplonger %>%
  separate(state, into = c("trait_state", "FS"), sep = "_FS") 

total_times_trait <- df_long %>%
  group_by(tree, trait1, trait2, trait_state, Which) %>%
  summarize(total_trait = sum(proportion), .groups = "drop")

total_times_FS <- df_long %>%
  group_by(tree, trait1, trait2, FS, Which) %>%
  summarize(total_FS = sum(proportion), .groups = "drop")

df_long <- df_long %>%
  left_join(total_times_trait, by = c("tree", "trait1", "trait2", "trait_state", "Which")) %>%
  left_join(total_times_FS, by = c("tree", "trait1", "trait2", "FS", "Which"))

df_long <- df_long %>%
  mutate(expected_proportion = total_trait * total_FS)

df_real <- df_long %>%
  filter(Which == "Real")
df_real <- df_real %>%
  mutate(abs_diff = abs(proportion - expected_proportion))

D_real <- df_real %>%
  summarize(total_abs_diff = sum(abs_diff)) %>%
  pull(total_abs_diff) / nsims

Real_dsims <- df_real %>%
  select(tree, trait1, trait2, trait_state, FS, abs_diff) %>%
  group_by(tree, trait1, trait2) %>%
  summarize(Real_dsim = sum(abs_diff), .groups = "drop")

# Dummy data
# Step 1: Filter for "Dummy" data
df_dummy <- df_long %>%
  filter(Which == "Dummy")

# Step 2: Calculate the absolute differences
df_dummy <- df_dummy %>%
  mutate(abs_diff = abs(proportion - expected_proportion))

# Step 3: Calculate Dummy_dsums (row-wise sums of the absolute differences)
Dummy_dsums <- df_dummy %>%
  group_by(tree, trait1, trait2) %>%
  summarize(Dummy_dsum = sum(abs_diff), .groups = "drop")

#hist(Real_dsims$Real_dsim)
#hist(Dummy_dsums$Dummy_dsum)
#abline(v = D_real)
numGreater = sum(Dummy_dsums$Dummy_dsum > D_real)
pval = sum(Dummy_dsums$Dummy_dsum > D_real)/nsims

## Get num Dummy greater than median real
# Calculate medians for "Real" data
medians_real <- df_long %>%
  filter(Which == "Real") %>%
  group_by(trait_state, FS) %>%
  summarize(median_real = median(proportion), .groups = "drop")

# Filter for "Dummy" data
df_dummy <- df_long %>%
  filter(Which == "Dummy")

# Join the medians back to the "Dummy" data
df_dummy <- df_dummy %>%
  left_join(medians_real, by = c("trait_state", "FS"))

# Calculate the fraction for each state in "Dummy" data
fraction_dummy_less_than_median_real <- df_dummy %>%
  group_by(trait_state, FS) %>%
  summarize(fraction = sum(proportion <= median_real) / n(), .groups = "drop")


trait1 = overlaplonger$trait1[1]
trait2 = overlaplonger$trait2[1]

plotlabel = paste(trait1, trait2)
dummytitle = paste("Nsims =", nsims, "\nnum Dummy dsums > D_real:", numGreater, ", pval =", pval)

xmax = max(c(Real_dsims$Real_dsim, Dummy_dsums$Dummy_dsum))*1.1

# ggplot histograms to return
p1 <- ggplot(Real_dsims, aes(x = Real_dsim)) +
  geom_histogram(binwidth = xmax/20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
  geom_vline(xintercept = D_real, color = "red") +
  xlim(c(0, xmax)) +
  labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
  theme_minimal(base_size = 10)

# Create the histogram for Dummy_dsums
p2 <- ggplot(Dummy_dsums, aes(x = Dummy_dsum)) +
  geom_histogram(binwidth = xmax/20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
  xlim(c(0, xmax)) +
  labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
  theme_minimal(base_size = 10)

#pdf(paste("simmap overlap Real Dummy multiState Proportions Boxplot -", trait1, trait2, nsims, "sims.pdf"))
boxplotStates <- ggplot(overlaplonger, aes(x = state, y = proportion, fill = Which)) +
  geom_boxplot(outlier.shape = NA) + # Exclude outliers
  theme_minimal() +
  labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
  scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste(trait1, trait2, "p =", pval))
#dev.off()
  
  return(list(p1=p1, p2=p2, boxplotStates = boxplotStates, fraction_dummy_less_than_median_real = fraction_dummy_less_than_median_real))

}


plot_my_simmaps <- function(thistrait, subsettrait, simmappy, otherlabel, discretetraitvecList) {
    simmapFileName = paste0(thistrait, " multistate simmap plots ", subsettrait, " subset ", otherlabel, Sys.Date(),".pdf")
    pdf(simmapFileName, height = 9, width = 12)
    par(mfrow = c(2,3))
    par(mar = c(3.8,3.8,3,1))
    for (simmapNum in 1:5) {
      simmap = simmappy[[simmapNum]]
      discretetraitvec = discretetraitvecList[[simmapNum]]
      simmapQ = simmap$Q
      SimmapStates = colnames(simmap$Q)
      py = c("orange", "#009E73", "blue", "#CC79A7")
      pynamed <- py[1:length(SimmapStates)]
      names(pynamed) <- SimmapStates
      tipcols = rep(NA, length(simmap$tip.label))
      for (stateNum in 1:length(SimmapStates)) { # get color vector of tips
        tempstate = SimmapStates[stateNum]
        tipcols[which(simmap$tip.label %in% names(which(discretetraitvec==tempstate)))] <- pynamed[tempstate]
      }
      # Plot simmap 
      plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
      tiplabels(pch=21,bg=tipcols, col = tipcols, cex=0.3)
      legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
    }
    # plot to add rate matrix
    # Create an empty plot
    par(mar = c(5,5,4,1))
    plot(1, type = "n", xlim = c(0, ncol(simmapQ)+1), ylim = c(0, nrow(simmapQ)+1), 
         xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    
    # Add column and row names
    axis(1, at = 1:ncol(simmapQ), labels = colnames(simmapQ), las = 2, tick = FALSE)
    axis(2, at = 1:nrow(simmapQ), labels = rev(rownames(simmapQ)), las = 2, tick = FALSE)
    # Add the matrix values
    for (i in 1:nrow(simmapQ)) {
      for (j in 1:ncol(simmapQ)) {
        text(j, nrow(simmapQ) - i + 1, round(simmapQ[i, j], 6))
      }
    }
    # Add label re transitions
    axis(1, at = 0.1, labels = "To:", las = 1, tick = FALSE, font = 2)
    axis(2, at = nrow(simmapQ)+0.9, labels = "From:", las = 2, tick = FALSE, font = 2)
    title("Transition Rates")
    dev.off()
  
}



#### misc ----
dfout = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine .csv")
dfDummy = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ HighConfidence_Coop FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine .csv")
calcHuelout = calcHuel(dfout, dfDummy)
calcHuelout$TransitionStats
calcHuelout$mediansRow
dfCounts = calcHuelout$dfMeltCounts

colnames(dfCounts)

