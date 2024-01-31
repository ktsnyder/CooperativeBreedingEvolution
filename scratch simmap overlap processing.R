
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


## make tables of simmap overlap
socTraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.MoreThanTwoCaretakers", "MeanCoopTie2Coop") 
otherTraits = c("FemaleSong_Agg01","MeanCoopTie2Noncoop")

outlist = list()
boxplotlist = list()
outdf = set.seed(10)
for (i in 1:length(socTraits)) {
  for (j in 1:2) {
  temptrait1 = socTraits[i]
  temptrait2 = otherTraits[j]
  
  filelist = list.files("Simmap Overlap Outputs", pattern = ".csv")  
  filelist = filelist[which(str_detect(filelist, temptrait1) & str_detect(filelist, temptrait2))]
  
  if (length(filelist > 0)) {
    dummyfile = filelist[which(str_detect(filelist, "DUMMY"))]
    realfile = filelist[which(str_detect(filelist, "DUMMY", negate = T))]
    if (length(dummyfile) == 1 & length(realfile) == 1) {
      dfDummy1 = read.csv(paste0("Simmap Overlap Outputs/",dummyfile))
      dfout1 = read.csv(paste0("Simmap Overlap Outputs/",realfile))
    } else {
      print(dummyfile)
      print(realfile)
      dummyfile = dummyfile[which(str_detect(dummyfile, "2000"))]
      realfile = realfile[which(str_detect(realfile, "counts"))]
      dfDummy1 = read.csv(paste0("Simmap Overlap Outputs/",dummyfile))
      dfout1 = read.csv(paste0("Simmap Overlap Outputs/",realfile))
      #next
    }
  } else {
    print(paste("no files found", temptrait1, temptrait2))
    next
  }
  
  calcHuelout = calcHuel(dfout = dfout1, dfDummy = dfDummy1)
  outlist[[length(outlist)+1]] <- calcHuelout
  boxplotlist[[length(boxplotlist)+1]] <- calcHuelout$p3
  
  temprow = calcHuelout$mediansRow
  outdf = rbind(outdf, temprow)

  }
}
write.csv(outdf, file = "simmap overlap summary table.csv", row.names = FALSE)
nplotrows = ceiling(length(boxplotlist)/2)
single_page_plotBox <- grid.arrange(grobs = boxplotlist, ncol = 2, nrow = nplotrows)
ggsave("simmap overlap all boxplots_onepage.pdf", single_page_plotBox, width = 14, height = 38, units = "in")
