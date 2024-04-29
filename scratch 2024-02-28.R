setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
library(stringr)
library(dplyr)
source("test_trait_overlap_simmaps.R")
source("transition_plot.R")


#### table ----
socTraits = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.MoreThanTwoCaretakers") # , "Griesser2023.TwoOrMoreCaretakers" - removed because too few species in lower group
# transition plots not in main text: socTraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.Colonial01", "Final.polygyny", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
otherTraits = c("HighConfidence_Coop")
groupLabel = "HighConfidence-Coop_Sociality"

socTraits = c("HighConfidence_Coop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "JetzCoop", "DaleCoop", "DowningCoop", "BiagoliniCoop", "RubensteinCoop")
otherTraits = "FemaleSong_Agg01"
groupLabel = "SingleSourceCoop_FemaleSong"

socTraits = c("HighConfidence_Coop")
otherTraits = "FemaleSong_Agg01"
groupLabel = "Jackknifed_HighConfCoop_FemaleSong"
jackknife = TRUE

outlist = list()
boxplotlist = list()
transplotlist = list()
outdf = set.seed(10)
for (i in 1:length(socTraits)) {
  for (j in 1:length(otherTraits)) {
    temptrait1 = socTraits[i]
    temptrait2 = otherTraits[j]
    
    filelist = list.files("Simmap Overlap Outputs", pattern = ".csv")  
    
    if (jackknife == TRUE) {
      filelist = filelist[which(str_detect(filelist, temptrait1) & str_detect(filelist, temptrait2) & str_detect(filelist, "remove"))]
      tempfamilies <- sapply(filelist, function(name) {
        matches <- regmatches(name, gregexpr("(?<=removed)[A-Za-z]+", name, perl = TRUE))
        matches[[1]]
      })
      familyvec = unique(tempfamilies)
    } else { # if jackknife != TRUE
      filelist = filelist[which(str_detect(filelist, temptrait1) & str_detect(filelist, temptrait2) & !str_detect(filelist, "remove"))]
      
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
      familyvec = "None"
    } # end if jackknife != TRUE
    
    for (tempfam in familyvec) {
      
      if (jackknife == TRUE) {
        tempotherlabel =  paste0("removed", tempfam)
        print(tempotherlabel)
        filelistFam = filelist[which(str_detect(filelist, tempotherlabel))]
        dummyfile = filelistFam[which(str_detect(filelistFam, "DUMMY"))]
        realfile = filelistFam[which(str_detect(filelistFam, "DUMMY", negate = T))]
        dfout1 = read.csv(paste0("Simmap Overlap Outputs/",realfile))
        dfDummy1 = read.csv(paste0("Simmap Overlap Outputs/",dummyfile))
      } else {
        tempotherlabel = ""
      }
      
      nsims = length(dfDummy1$treenum)
      
      calcHuelout = calcHuel(dfout = dfout1, dfDummy = dfDummy1, otherlabel = tempotherlabel)
      outlist[[length(outlist)+1]] <- calcHuelout
      boxplotlist[[length(boxplotlist)+1]] <- calcHuelout$p3
      SimmapOverlapPval = calcHuelout$mediansRow$pval
      
      ObsProp0Absent_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp0Absent <= median(dfout1$ObsProp0Absent))/length(dfDummy1$treenum)
      ObsProp0Present_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp0Present <= median(dfout1$ObsProp0Present))/length(dfDummy1$treenum)
      ObsProp1Absent_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp1Absent <= median(dfout1$ObsProp1Absent))/length(dfDummy1$treenum)
      ObsProp1Present_FractionDummyLessThanMedianReal = sum(dfDummy1$ObsProp1Present <= median(dfout1$ObsProp1Present))/length(dfDummy1$treenum)
      
      temprow = cbind(calcHuelout$mediansRow, ObsProp0Absent_FractionDummyLessThanMedianReal, ObsProp0Present_FractionDummyLessThanMedianReal, ObsProp1Absent_FractionDummyLessThanMedianReal, ObsProp1Present_FractionDummyLessThanMedianReal, tempfam)
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
      colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop0DifferenceFromExpected")] <- "q12"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS0DifferenceFromExpected")] <- "q13"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop0DifferenceFromExpected")] <- "q21"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS1DifferenceFromExpected")] <- "q24"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS0DifferenceFromExpected")] <- "q31"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop1DifferenceFromExpected")] <- "q34"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS1DifferenceFromExpected")] <- "q42"
      colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop1DifferenceFromExpected")] <- "q43"
      
      if (trait2 == "MeanCoopTie2Noncoop") {
        labx0 = paste("Non-Cooperative")
        labx1 = paste("Cooperative")
      } else if (trait2 == "FemaleSong_Agg01") {
        labx0 = paste("Female Song Absent")
        labx1 = paste("Female Song Present")
      } else if (trait2 == "HighConfidence_Coop") {
        labx0 = paste("Non-Cooperative")
        labx1 = paste("Cooperative")
      } else if (trait2 == "Final.polygyny") {
        labx0 = paste("Monogamy")
        labx1 = paste("Polygyny")
      } else if (grepl("coop", trait2, ignore.case = T)) {  #(str_detect(trait2, "coop")) {
        labx0 = paste("Non-Cooperative")
        labx1 = paste("Cooperative")
      } else if (str_detect(trait2, "Kin")) {
        labx0 = paste("Non-kin")
        labx1 = paste("Kin")
      } else if (str_detect(trait2, "Familial")) {
        labx0 = paste("Non-Familial Living")
        labx1 = paste("Familial Living")
      } else if (str_detect(trait2, "Colonial")) {
        labx0 = paste("Non-Colonial")
        labx1 = paste("Colonial") 
      } else if (str_detect(trait2, "GroupsLargerThanPair")) {
        labx0 = paste("Asocial or pair")
        labx1 = paste("Small or large groups") 
      } else if (str_detect(trait2, "LongSocialBonds")) {
        labx0 = paste("Bonds last one season or less")
        labx1 = paste("Multi-year bonds") 
      } else if (str_detect(trait2, "MoreThanTwoCaretakers")) {
        labx0 = paste("Two or fewer caretakers")
        labx1 = paste("More than two caretakers") 
      } else if (str_detect(trait2, "TwoOrMoreCaretakers")) {
        labx0 = paste("Fewer than two caretakers")
        labx1 = paste("Two or more caretakers") 
      } else if (str_detect(trait2, "Asocial")) {
        labx0 = paste("Asocial")
        labx1 = paste("Pair or group sociality") 
      } else if (str_detect(trait2,"SeasonOrLonger")) {
        labx0 = paste("Bonds last less than one season")
        labx1 = paste("Season or longer social bonds") 
      } else if (str_detect(trait2, "LargestGroupSizes")) {
        labx0 = paste("Asocial, pair, or small groups")
        labx1 = paste("Large groups") 
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
      } else if (str_detect(trait1, "FemaleSong")) {
        lab0x = paste("Female Song Absent")
        lab1x = paste("Female Song Present")
      }
      
     # print("ratePvals:")
    #  print(ratePvals)
      
      trait1StateLabels = c(lab0x, lab1x)
      trait2StateLabels = c(labx0, labx1)
      plottitle = paste(trait1, trait2, "\nnsims:", nsims, tempotherlabel, "pval =", SimmapOverlapPval)
      transplotout = transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 1, offset = 0.2, lengthen = 0.2, ratePvals = ratePvals, plottitle = plottitle, center="median")
      transplotlist[[length(transplotlist) + 1]] = transplotout$transitionplot_GrayNS
    } # end for tempfam in familyvec
  }
}
write.csv(outdf, file = paste(groupLabel,"simmap overlap summary table_StateCompare.csv"), row.names = FALSE)
nplotrows = ceiling(length(boxplotlist)/2)
require(gridExtra)
single_page_plotBox <- grid.arrange(grobs = boxplotlist, ncol = 2, nrow = nplotrows)
ggsave(paste(groupLabel,"simmap overlap all boxplots_onepage_500dummySims_testFactorLabs.pdf"), single_page_plotBox, width = 14, height = 4.5*nplotrows, units = "in", limitsize = FALSE)

nplotrows = ceiling(length(transplotlist)/2)
single_page_plotTrans <- grid.arrange(grobs = transplotlist, ncol = 2, nrow = nplotrows)
ggsave(paste(groupLabel,"simmap overlap counts All transition plots_onepage.pdf"), single_page_plotTrans, width = 18, height = 4.5*nplotrows, units = "in", limitsize = FALSE)





# 
library(phytools)
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
tree = read.nexus(treefile)
df = read.csv(newdata)
dfOscine = df[which(df$species %in% tree$tip.label),]
sort(unique(dfOscine$Family3_BirdtreeMatchSpecies2_AVONET[which(!is.na(dfOscine$FemaleSong_Agg01))]))
