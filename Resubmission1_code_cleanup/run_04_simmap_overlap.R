### Overlapping stochastic character maps to assess co-occurrence of discrete trait states in evolutionary history ----

# Figure 3A; Supplemental Table 9 - Co-occurrance of cooperative breeding/familial living and female song  ----
source("test_trait_overlap_simmaps.R")
nsims_real = 10
nsims_dummy = 10

targetMetrics = c("HighConfidence_Coop", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "JetzCoopInclCockburn", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnInferred")
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

# Supplemental Table 7 - Co-occurrance jackknife  ----
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
templabel = ""
columns = c(trait1, trait2)

nsims_real = 5
nsims_dummy = 10

source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = treefile)
qrates= Qout$qrates

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 65)]
familyvec = c("None", familyvec)

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


# Extended Data Figure 5A, Extended Data Table 1 and Supplemental Table 10 - Co-occurrance of female song (or cooperative breeding) with sociality traits ----
source("test_trait_overlap_simmaps.R")
nsims_real = 20
nsims_dummy = 20

socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Territory_12vs3", "TerritorialityWeakVsStrong")

temptrait2 = "FemaleSong_Agg01"
#temptrait2 = "HighConfidence_Coop"  # uncomment to run analyses found in Supplemental Table 10

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
