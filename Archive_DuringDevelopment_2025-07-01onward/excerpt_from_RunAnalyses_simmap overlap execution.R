## excerpted from Run_Analyses.R 
## 



# Figure 3A & B; Supplemental Table 9 - Co-occurrance of cooperative breeding/familial living and female song  ----
source("test_trait_overlap_simmaps.R")
nsims_real = 10
nsims_dummy = 10

targetMetrics = socialityMetrics = c("HighConfidence_Coop", "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "TerritorialityWeakVsStrong", "Territory_12vs3")
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