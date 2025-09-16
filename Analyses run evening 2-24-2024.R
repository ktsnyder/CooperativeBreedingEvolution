# run late night 2/24/2024

setwd("~/Desktop/CooperativeBreedingEvolution")
newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

newdf = read.csv(newdata)
newdf$System01_Jetz2011 = NA
newdf$System01_Jetz2011[which(newdf$System_Jetz2011 == "Non-cooperative")] <- 0
newdf$System01_Jetz2011[which(newdf$System_Jetz2011 == "Cooperative")] <- 1
tempMetric = "System01_Jetz2011"

# Coop and female song simmap overlap
source("test_trait_overlap_simmaps.R")
dfout4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdf, tree =  treefile, dummy = FALSE, nsims = 500, treelabel = "HackettOscine", datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdf, tree =  treefile, dummy = TRUE, nsims = 500, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
calcHuelout = calcHuel(dfout4, dfDummy4)

nsims_real = length(dfout4[,1])    
nsims_dummy = length(dfDummy4[,1])

require(gridExtra)  
treelabel = "HackettOscine"
tempMetric = "System_Jetz2011" #"HighConfidence_Coop"
plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts-ExpObsCompare.pdf")  
nPlots = 6 #length(calcHuelout)-4  
print(calcHuelout$mediansRow)  
m3 <- marrangeGrob(calcHuelout[c("p1", "p2", "p3", "p4", "p5", "p6")], ncol = 1, nrow = nPlots) # nrow can be nPlots if no "filename" or "TransitionStats" in calcHuelout  
ggsave(plotname, m3, width = 7.5, height = 4.2*nPlots, units = "in")


# Coop and sociality metrics

socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny",  "Griesser2023.MoreThanTwoCaretakers")
nsims_real = 500
nsims_dummy = 500

source("test_trait_overlap_simmaps.R")
for (i in 2: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  dfout1 <- CharacterSimmaps(columns = c(tempMetric,"HighConfidence_Coop"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = FALSE)
  dfDummy1 <- CharacterSimmaps(columns = c(tempMetric,"HighConfidence_Coop"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  calcHuelout = calcHuel(dfout1, dfDummy1)
  require(gridExtra)
  plotname = paste0("Simmap Overlap Outputs/",tempMetric, " HighConfidence_Coop ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts-ExpObsCompare.pdf")
  nPlots = 6# length(calcHuelout)-5
  m3 <- marrangeGrob(calcHuelout[c("p1", "p2", "p3", "p4", "p5", "p6")], ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 4.2*nPlots, units = "in")
}

#### Brownie - coop and song features ----
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Song.rep.final", "Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final", "Song.rep.min", "Song.rep.max")
CBcolumn <- "HighConfidence_Coop"
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

# process brownie outputs
filelist = list.files("OutputFiles/Brownie", pattern = "2024-02",full.names = T)
filelist = filelist[which(str_detect(filelist, ".csv"))]

browniesummarydf = set.seed(10)
for (i in 1:length(filelist)) {
  tempfile = filelist[i]
  df = read.csv(tempfile) 
  nsims = length(df$Pval)
  percentSig = sum(df$Pval < 0.05)/nsims
  DiscreteTrait = df$DiscreteTrait[1]
  SongFeature = df$ContinuousTrait[1]
  PercentSimsRate0GreaterThanRate1 = sum(df$ARDRate0 > df$ARDRate1)/nsims
  MedianRate0 = median(df$ARDRate0)
  MedianRate1 = median(df$ARDRate1)
  MedianPval = median(df$Pval)
  temprow = c(tempfile, DiscreteTrait, SongFeature, nsims, MedianPval, percentSig, PercentSimsRate0GreaterThanRate1, MedianRate0, MedianRate1)
  names(temprow) <- sapply(substitute(list(tempfile, DiscreteTrait, SongFeature, nsims, MedianPval, percentSig, PercentSimsRate0GreaterThanRate1, MedianRate0, MedianRate1))[-1], deparse)
  browniesummarydf = rbind(browniesummarydf, temprow)
}
browniesummarydf = as.data.frame(browniesummarydf)
write.csv(browniesummarydf, "HighConfidence_Coop brownie output summary.csv", row.names = F)


#### phylANOVA and ScatterBoxes ----
source("scatterboxes.R")
socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny",  "Griesser2023.MoreThanTwoCaretakers", "HighConfidence_Coop")
discLabelList = list(c("Asocial", "Social"), c("Groups Pair or Smaller", "Groups Larger than Pair"), c("Smaller groups", "Largest group sizes"), c("Season or shorter social bonds", "Longest social bonds"), c("Non-familial", "Familial"), c("Non-Colonial", "Colonial"), c("Monogamy", "Polygyny"), c("Two or fewer caretakers", "More than two caretakers"), c("Noncooperative", "Cooperative"))

for (i in 9:length(socialityMetrics)) {
  CBcolumn = socialityMetrics[i]
  tempLabels = discLabelList[[i]]
  scatterboxes(DiscreteTrait = CBcolumn, newdata = newdata, newtree = treefile, otherlabel = currentlabel, discreteCategoryLabels = tempLabels)
}


#### Generate counts tables ----

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
dfIn = read.csv(newdata)
dfIn$AnyNoncoopEqualsNoncoop = dfIn$MeanCoopTie2Noncoop
dfIn$AnyNoncoopEqualsNoncoop[which(dfIn$SourceDiscrepancy == 1)] = 0
OscineSubset = subsettreedata(newdata = dfIn, newtree = treefile)
df = OscineSubset$subsetdf
columnsToSummarize = c("HighConfidence_Coop", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop", "MeanCoopOmitTies", "AnyNoncoopEqualsNoncoop", "Griesser2023.Colonial01", "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.TwoOrMoreCaretakers", "Griesser2017FamilialLiving", "Final.polygyny")
otheraxiscolumns = c("FemaleSong_Agg01")

 otheraxiscolumns= columnsToSummarize = c("FemaleSong_Agg01", "HighConfidence_Coop", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop", "MeanCoopOmitTies", "AnyNoncoopEqualsNoncoop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "DaleCoop", "DowningCoop", "BiagoliniCoop", "System_Jetz2011")

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
      oneColGrouped = df %>% group_by(.data[[othertempCol]]) %>% summarize(Count = n(), .groups='drop')
      oneColGrouped$Count[which(oneColGrouped[,1] == 0)]
      wide_df = as.data.frame(matrix(data = c(oneColGrouped$Count[which(oneColGrouped[,1] == 0)],NA, NA,oneColGrouped$Count[which(oneColGrouped[,1] == 1)]), nrow = 2, ncol = 2)) 
    }
    rownames(wide_df) = paste(othertempCol, c(0,1), sep = ".")
    wide_df = wide_df[,which(colnames(wide_df) != othertempCol)]
    colnames(wide_df) = paste(tempCol, c(0,1), sep = ".")
    tablesegment = rbind(tablesegment, wide_df)
  }
  tablesegment = tablesegment[which(rownames(tablesegment) != 1),]
  fulltable = cbind(fulltable, tablesegment)
}
fulltable = fulltable[,which(colnames(fulltable) != "fulltable")]
fulltable
write.csv(fulltable, "SuppTable_binary traits state intersections NumSpecies_AllSourceCoops.csv", row.names = T)



#### Simmap Overlap Source-specific Coop classifications and female song
CoopCols = c("Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "JetzCoop", "DaleCoop", "DowningCoop", "BiagoliniCoop", "RubensteinCoop")

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
dfIn = read.csv(newdata)
dfIn$AnyNoncoopEqualsNoncoop = dfIn$MeanCoopTie2Noncoop
dfIn$AnyNoncoopEqualsNoncoop[which(dfIn$SourceDiscrepancy == 1)] = 0
CoopCols = c("HighConfidence_Coop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "MeanCoopTie2Noncoop", "Griesser2017Coop", "CornwallisCoop", "CockburnCoop", "JetzCoop", "DaleCoop", "DowningCoop", "BiagoliniCoop")
newdata = dfIn
nsims_real = 500
nsims_dummy = 500

for (i in 2: length(CoopCols)) {
  tempMetric = CoopCols[i]
  print(i)
  print(tempMetric)
  dfout2 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = FALSE)
  dfDummy2 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  print(paste("starting calcHuel", tempMetric, Sys.time()))
  calcHuelout = calcHuel(dfout2, dfDummy2)
  print(paste("starting simmap overlap plots", tempMetric, Sys.time()))
  require(gridExtra)
  plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts-ExpObsCompare.pdf")
  nPlots = 6# length(calcHuelout)-5
  m3 <- marrangeGrob(calcHuelout[c("p1", "p2", "p3", "p4", "p5", "p6")], ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 4.2*nPlots, units = "in")
}
