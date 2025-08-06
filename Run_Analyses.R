## Run analyses to generate figures in manuscript 
## Kate T Snyder
## Edited 6/19/2024
## Edited 11/27/2024 - added new Cooperative Breeding classification column based on HighConfidence_Coop but defaulting to Cockburn "Inferred", also used CockburnInferred
## Edited 12/9/2024 - added Tobias Territory
## Edited 2/3/2025 - added tests of character overlap within just "weakly territorial" and "year-round territorial" species subsets
## Edited June 8, 2025 - updating to be the overall runner script 
# Have added: Three-state simmap transitions (Female Song, Cooperative Breeding, Territoriality); Bias tests and downsampling calculations; phylopath and phylopath with downsampling
# Edited 6/9/2025 - removed calculation of plumage dimorphism (now done in merge_data_allcolumns.R)

#  Need to add: phyloglm, brownie resampled, territory data table, brownie with territory multistate, 


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
source("TransitionCounts_3trait_flexTerr_fxns.R")



# Set data, tree, and variable vectors ----
newdata = "Data_R_2025-06-09.csv"
df = read.csv(newdata)
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex" # to perform test using the alternative consensus tree for any given analysis, replace this with "ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex"
# treefile can also be set to any individual tree extracted from a BirdTree multiphylo object in order to test across many individual trees
songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")
socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "grouping_Griesser2023", "social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "social_bonds_Griesser2023", "Territory", "Territory_12vs3", "TerritorialityWeakVsStrong")
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"

# New CockburnInferred-using Cooperative Breeding variable
df = read.csv(newdata)
df$HighConf_Coop_DefaultToCockburnInferred <- df$HighConfidence_Coop
df$HighConf_Coop_DefaultToCockburnInferred[which(!is.na(df$CockburnInferred))] <- df$CockburnInferred[which(!is.na(df$CockburnInferred))]
df %>% group_by(HighConfidence_Coop, HighConf_Coop_DefaultToCockburnInferred) %>% count # Changes the classifications of 17 species, adds classifications to 169 species
df %>% group_by(FemaleSong_Agg01, HighConfidence_Coop, HighConf_Coop_DefaultToCockburnInferred) %>% count # changes CB classification for 6 species for which there is FS data, gain CB classification for 35 species that had FS data but no CB data

# Compare Female Song x Coop Breed across the old vs Cockburn-oriented cooperative breeding classification schemes
df_high_conf <- df %>% 
  group_by(HighConfidence_Coop, FemaleSong_Agg01) %>% 
  count() %>% 
  rename(HighConfidence_Coop_n = n)

df_high_conf_default <- df %>% 
  group_by(HighConf_Coop_DefaultToCockburnInferred, FemaleSong_Agg01) %>% 
  count() %>% 
  rename(HighConf_Coop_DefaultToCockburnInferred_n = n)

df_cockburn <- df %>% 
  group_by(CockburnInferred, FemaleSong_Agg01) %>% 
  count() %>% 
  rename(CockburnInferred_n = n)

combined_df <- full_join(
  df_high_conf, 
  df_high_conf_default, 
  by = c("HighConfidence_Coop" = "HighConf_Coop_DefaultToCockburnInferred", "FemaleSong_Agg01")
) %>%
  full_join(
    df_cockburn, 
    by = c("HighConfidence_Coop" = "CockburnInferred", "FemaleSong_Agg01")
  )
combined_df <- combined_df %>%
  rename(Coop = HighConfidence_Coop)
print(combined_df, n = Inf)
#write.csv(combined_df, "species counts comparison - FS vs CockburnInferred Coop Classifications.csv", row.names = F)


# Add abs_Latitude, Migration_num, and logMass_normalized variables
df <- read.csv("Data_R_2025-06-09.csv")
df$Migration_num <- as.numeric(df$Migration_AVONET)-2
df$abs_Latitude <- abs(df$Centroid.Latitude_AVONET)
df$abs_Latitude_normalized <- scale(df$abs_Latitude, center = TRUE, scale = TRUE)
df$logMass_normalized <- scale(df$logMass_AVONET, center = TRUE, scale = TRUE)
df$PercentAbsLogWingDimorphism_normalized <- scale(df$PercentAbsLogWingDimorphism, center = TRUE, scale = TRUE)
df$logMaleFemalePlumageDiffAbs_normalized <- scale(df$logMaleFemalePlumageDiffAbs, center = TRUE, scale = TRUE)
#write.csv(df, "Data_R_2025-07-21.csv", row.names = FALSE)

newdata = "Data_R_2025-07-21.csv"
# Add DefaultToCockburnInferred variable
df = read.csv(newdata)
df$HighConf_Coop_DefaultToCockburnInferred <- df$HighConfidence_Coop
df$HighConf_Coop_DefaultToCockburnInferred[which(!is.na(df$CockburnInferred))] <- df$CockburnInferred[which(!is.na(df$CockburnInferred))]
#write.csv(df, "Data_R_2025-07-23.csv", row.names = FALSE)

newdata = df

# Supplemental Table 3 - phylANOVA ----
phynovaDF = set.seed(10)
for (j in 1:length(songtraits)) {
  for (k in 1:length(socialityMetrics)) {
  songtrait = songtraits[j]
  tempgrouptrait = socialityMetrics[k]
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
      phylANOVAout = phylANOVA(subsettree, x = discvec, y = contvec, nsim = 50000, posthoc = TRUE)
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
}
write.csv(phynovaDF, file = paste(Sys.Date(), "phylANOVA outputs Songs.csv"), row.names = F)






# Figure 1, Supplemental Figure 1, Supplemental Table 4 - Brownie ----
source("brownie relative rates.R")
discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 500
currentlabel = "newBinaryTerrs"
CBcolumn = "HighConfidence_Coop"
songtraits = "Song.rep.final"
socialtraits = c("TerritorialityWeakVsStrong", "Territory_12vs3")

for (l in 1:length(socialtraits)) {
for (k in 1:length(songtraits)) {
  print(Sys.time())
  feature <- songtraits[k]
  CBcolumn <- socialtraits[l]
  discreteCatLabels = getLabels(CBcolumn)
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
} # end for l
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
source("BrownieMultistate.R")
BrownieMultistate(DiscreteTrait = "grouping_Griesser2023", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = 5, plotsimmaps = F, plotResults = T, otherlabel = "test")
BrownieMultistate(DiscreteTrait = "social_system_incl_nk_coop_Griesser2017", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = 10, plotsimmaps = F, plotResults = T, otherlabel = "test")

df_brown = read.csv(newdata)
df_brown$TerritorialityWeakVsStrong[which(df_brown$TerritorialityWeakVsStrong == 0)] <- "Weak"
df_brown$TerritorialityWeakVsStrong[which(df_brown$TerritorialityWeakVsStrong == 1)] <- "Strong"
df_brown$HighConfidence_Coop[which(df_brown$HighConfidence_Coop == 0)] <- "Noncooperative"
df_brown$HighConfidence_Coop[which(df_brown$HighConfidence_Coop == 1)] <- "Cooperative"
df_brown$TerrWeakStrongXHighConfCoop <- paste(df_brown$HighConfidence_Coop, df_brown$TerritorialityWeakVsStrong, sep = "_")
require(stringr)
df_brown$TerrWeakStrongXHighConfCoop[which(str_detect(df_brown$TerrWeakStrongXHighConfCoop, "NA"))] <- NA

tree = read.nexus(treefile)
TerrQout = findQrates(columns = "TerritorialityWeakVsStrong", plot = FALSE, newtree = treefile, newdata = df_brown)
CoopQout = findQrates(columns = "HighConfidence_Coop", plot = FALSE, newtree = treefile, newdata = df_brown)
TerrQ = TerrQout$qrates
CoopQ = CoopQout$qrates

nsims = 1000
subsetdf = df_brown[complete.cases(df_brown[,c("TerritorialityWeakVsStrong", "HighConfidence_Coop")]),]
# Find tips to drop (those not in the dataframe)
tips_to_drop <- setdiff(tree$tip.label, subsetdf$species)
# Drop tips using phytools function
pruned_tree <- drop.tip(tree, tips_to_drop)

Terrtraitvec = subsetdf[,"TerritorialityWeakVsStrong"]
names(Terrtraitvec) = subsetdf$species

Cooptraitvec = subsetdf[,"HighConfidence_Coop"]
names(Cooptraitvec) = subsetdf$species

simmapTerritoryWeakStrong = make.simmap(pruned_tree, Terrtraitvec, nsim = nsims, Q= TerrQ, type = "discrete") 

simmapHCCoop = make.simmap(pruned_tree, Cooptraitvec, nsim = nsims, Q= CoopQ, type = "discrete")

TerrCoopsims = list()
for (i in 1:nsims) {
  TerrCoopsims[[i]] <- merge_simmaps(simmap1 = simmapTerritoryWeakStrong[[i]], simmap2 = simmapHCCoop[[i]])
}
class(TerrCoopsims) <- c("multiSimmap", "multiPhylo")

BrownieMultistate(DiscreteTrait = "TerrWeakStrongXHighConfCoop", ContinuousTrait = "Song.rep.final", newdata = df_brown, treefile = treefile, nsim = 1000, plotsimmaps = T, plotResults = T, otherlabel = "MergedSimmaps", importSimmaps = TerrCoopsims)



# Figure 3A & B; Supplemental Table 9 - Co-occurrance of cooperative breeding/familial living and female song  ----
source("test_trait_overlap_simmaps.R")
nsims_real = 10
nsims_dummy = 10

targetMetrics = c("HighConfidence_Coop", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "JetzCoopInclCockburn", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop")
targetMetrics = c("HighConf_Coop_DefaultToCockburnInferred", "CockburnInferred")
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


# New simmap overlap test - FS x CB in territory == 2 vs territory == 3 adapting jackknife
newdata = "Data_R_Passerine_withTobias.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
trait1 = tempMetric = "HighConfidence_Coop"
trait2 = temptrait2 = "FemaleSong_Agg01" 
columns = c(tempMetric, temptrait2)
dfIn = read.csv(newdata)

nsims_real = 100
nsims_dummy = 100

targetMetrics = c(2,3)
for (i in 1: length(targetMetrics)) {
  subsetvalue = targetMetrics[i]
  territoryLabel = paste0("Territory",subsetvalue)
  templabel = paste0(columns[1], " ", columns[2], " ", territoryLabel)
  tempdfIn = dfIn[which(dfIn$Territory == subsetvalue),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  print(paste("Simmap Overlap", trait1, trait2, "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", setQratesTree = treefile, setQratesData = newdata, plotSampleSimmaps = T)
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", setQratesTree = treefile, setQratesData = newdata, plotSampleSimmaps = T)
  
  calcHuelOut = calcHuel(dfout = dfout4, dfDummy = dfDummy4, newplot = FALSE, otherlabel = templabel)
  
  require(gridExtra)
  plotname = file.path("Simmap Overlap Outputs", paste(templabel, nsims_real, nsims_dummy, "simmap overlap withTransCounts.pdf"))
  
  calcHuelout2 = calcHuelOut[c("p1", "p2", "p3", "p4","p5","p6")]
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
socialityMetrics = c("HighConfidence_Coop", "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "TerritorialityWeakVsStrong")
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
SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0VsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "Final.polygyny", "HighConfidence_Coop", "FemaleSong_Agg01", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Territory_12vs3", "TerritorialityWeakVsStrong")
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
multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "grouping_Griesser2023", "Territory")
multistateTraits = c("Territory", "Social.bond", "Social.bond")
othertraits = c("FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
nsims = 500
for (i in 1:length(multistateTraits)) {
  trait1 = multistateTrait = multistateTraits[i]
  trait2 = othertrait = othertraits[i]
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
  write.csv(overlapdf, file = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.csv"), row.names = F)
  
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.pdf")
  require(gridExtra)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:3], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 12, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = paste("simmap overlap FractionDummyLessThanMedianReal", multistateTrait, othertrait, nsims,"sims.csv"))
}



#### Run 3-trait transition counts and plot ----
columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong")
plot_transition_counts_3trait(Qdata = "Data_R_2025-06-09.csv", Qtree = '2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex', columns = columns, nsims = 10)


#### Phylopath analyses ----
source("run_phylopath_fxns.R")

# Prepare data for phylopath
df_phylo <- read.csv("Data_R_2025-06-09.csv")
tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Filter to tree species
dfIn_phylo <- df_phylo[df_phylo$species %in% tree_phylo$tip.label, ]
rownames(dfIn_phylo) <- dfIn_phylo$species

# Output directory
phylopath_output_dir <- "Outputs/PhylopathPlots"
if (!dir.exists(phylopath_output_dir)) dir.create(phylopath_output_dir, recursive = TRUE)

# 1. Phylopath with body mass
cat("\nRunning phylopath with body mass...\n")
result_phylopath <- run_CB_FS_Terr_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  female_song_var = "FemaleSong_Agg01",
  coop_breeding_var = "HighConfidence_Coop",
  territoriality_var = "TerritorialityWeakVsStrong",
  mass_var = "logMass_AVONET",
  plots2pdf = TRUE
)

plots_bodymass <- create_all_phylopath_plots(
  analysis_type = "nondownsampled",
  phylopath_output = result_phylopath,
  output_dir = phylopath_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

# 2. Phylopath with sexual dichromatism
if ("logMaleFemalePlumageDiffAbs" %in% colnames(dfIn_phylo)) {
  cat("\nRunning phylopath with sexual dichromatism...\n")
  result_phylopath <- run_CB_FS_Terr_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    female_song_var = "FemaleSong_Agg01",
    coop_breeding_var = "HighConfidence_Coop",
    territoriality_var = "TerritorialityWeakVsStrong",
    mass_var = "logMaleFemalePlumageDiffAbs",
    plots2pdf = TRUE
  )
  
  plots_dichrom <- create_all_phylopath_plots(
    analysis_type = "nondownsampled",
    phylopath_output = result_phylopath,
    output_dir = phylopath_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

# 3. Phylopath with sexual dimorphism
if ("PercentAbsLogWingDimorphism" %in% colnames(dfIn_phylo)) {
  cat("\nRunning phylopath with sexual dimorphism...\n")
  result_phylopath <- run_CB_FS_Terr_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    female_song_var = "FemaleSong_Agg01",
    coop_breeding_var = "HighConfidence_Coop",
    territoriality_var = "TerritorialityWeakVsStrong",
    mass_var = "PercentAbsLogWingDimorphism",
    plots2pdf = TRUE
  )
  
  plots_dimorph <- create_all_phylopath_plots(
    analysis_type = "nondownsampled",
    phylopath_output = result_phylopath,
    output_dir = phylopath_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

# 4. Phylopath with alternative cooperative breeding classifications
# Redone in RunAnalyses_nondownsampled_Phylopaths.R 7/23/2025
altCoops <- c("MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop")
for (tempCoop in altCoops) {
  if ("logMass_AVONET" %in% colnames(dfIn_phylo)) {
    cat("\nRunning phylopath with sexual dimorphism...\n")
    result_phylopath <- run_CB_FS_Terr_phylopath(
      dfIn = dfIn_phylo,
      tree = tree_phylo,
      female_song_var = "FemaleSong_Agg01",
      coop_breeding_var = tempCoop,
      territoriality_var = "TerritorialityWeakVsStrong",
      mass_var = "logMass_AVONET",
      plots2pdf = TRUE
    )
    
    plots_altCoops <- create_all_phylopath_plots(
      analysis_type = "nondownsampled",
      phylopath_output = result_phylopath,
      output_dir = phylopath_output_dir,
      save_png = TRUE,
      save_pdf = TRUE
    )
  }
}

# 5. Phylopath with binarized territory from Tobias et al (2016) body mass
cat("\nRunning phylopath with Terr12vs3 from Tobias et al (2016)...\n")
result_phylopath <- run_CB_FS_Terr_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  female_song_var = "FemaleSong_Agg01",
  coop_breeding_var = "HighConfidence_Coop",
  territoriality_var = "Territory_12vs3",
  mass_var = "logMass_AVONET",
  plots2pdf = TRUE
)

plots_Terr12v3 <- create_all_phylopath_plots(
  analysis_type = "nondownsampled",
  phylopath_output = result_phylopath,
  output_dir = phylopath_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)




#### Run brownie resampled min/max values ----
source("browniefunction.R")
source("findQrates.R")
source("plotbrownie.R")

dfIn <- read.csv("Data_R_2025-06-09.csv")
tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 2 # due to the way each simmap is called in browniefunction() (as one element in a list of simmaps), we can't do just 1 sim
currentlabel = "ResampleMinMedMax"
nSeeds = 500

# get Q rates out here so we can use the pared-down dataframe in the loop
Qoutput <- findQrates(columns = "HighConfidence_Coop", plot=F, newtree = tree, newdata = dfIn)
qrates <- Qoutput$qrates

# Resample song features and run brownie
for (song_base in c("Song.rep.", "Syllable.rep.")) {
  
  songcols <- grep(paste0("^",song_base), colnames(dfIn), value = TRUE)
  songcols
  
  # subset dataframe
  complete_vars <- c("HighConfidence_Coop", songcols)
  df <- dfIn[complete.cases(dfIn[,complete_vars]),]
  df = df[,c("species", complete_vars)]
  
  dfloop <- df
  outdf = set.seed(10)
  for (i in 1:nSeeds) {
    set.seed(i)
    sampleColName = paste0(song_base,"Sample",i)
    dfloop[[sampleColName]] <- apply(dfloop[, songcols], 1, sample, size = 1)
    brownieout = browniefunction(columns = c("HighConfidence_Coop", sampleColName), newdata = dfloop, newtree = treefile, nsim = nsim, islog = sampleColName, plotsimmaps = FALSE, setQrates = qrates)
    outdf = rbind(outdf, brownieout)
  }
  
  write.csv(outdf, paste0(Sys.Date(),"HighConfidence_Coop_",song_base, currentlabel, "_brownie",nsim,"sim", nSeeds, "resamples.csv"), row.names = F)
  write.csv(dfloop, paste0(Sys.Date()," HighConfidence_Coop_",song_base, currentlabel, " ", nSeeds, "resampled columns.csv"), row.names = F)
  
  plotbrownie(data = outdf, columns = c("HighConfidence_Coop",song_base), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nSeeds, islog = TRUE)
  
}

# Aggregate results
source("brownie relative rates.R")
# compile results into one table
#BrownieRelativeRates(BrownieOutputFolder = ".", otherlabel = currentlabel)
BrownieRelativeRates(BrownieOutputFolder = "OutputFiles", otherlabel = currentlabel)


#### Run bias analyses and calculate downsampling ----
source("generate_bias_report.R")
source("run_bias_tests.R")
source("run_dimorphism_bias_tests.R")
source("calculate_stratified_downsampling_with_territoriality.R")

df_bias = read.csv("Data_R_2025-06-09.csv")
tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# 1. Generate comprehensive bias report (all 6 tests)
cat("\nGenerating comprehensive bias report...\n")
generate_bias_report(
  df = df_bias,
  tree = tree,
  output_file = "Outputs/Bias_Test_Results.md",
  output_dir = "Outputs/"
)

# 2. Run standard bias tests for downsampling recommendations
cat("\nRunning bias tests for downsampling...\n")
bias_results <- run_bias_tests(
  df = df_bias,
  tree = tree,
  output_dir = "Outputs/BiasTests",
  geographic_col = "GeographicRegion_Jetz",
  territoriality_cols = c("TerritorialityWeakVsStrong", "Territory_12vs3"),
  save_plots = TRUE
)

# 3. Run dimorphism bias tests (optional - creates violin plots)
cat("\nRunning dimorphism bias tests...\n")
dimorphism_results <- run_dimorphism_bias_tests(
  df = df_bias,
  tree = tree,
  output_dir = file.path("Outputs","DimorphismBias")
)

# 4. Calculate stratified downsampling (for detailed bias correction)
cat("\nCalculating stratified downsampling for bias correction...\n")
downsampling_results <- calculate_stratified_downsampling(
  df = df_bias,
  stratify_vars = c("GeographicRegion_Jetz", "HighConfidence_Coop"),
  data_col = "FemaleSong_Agg01",
  output_file = file.path("Outputs", "PhylopathDownsampled", "Stratified_Downsampling_Calculations3.md")
)

cat("Downsampling calculations saved to: Outputs/Stratified_Downsampling_Calculations.md\n")
cat("\nSummary of species to remove:\n")
if (!is.null(downsampling_results$downsampling$holarctic_noncoop)) {
  cat("- Holarctic non-cooperative:", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "species\n")
} # should be 83
if (!is.null(downsampling_results$downsampling$tropical_coop)) {
  cat("- Tropical cooperative:", downsampling_results$downsampling$tropical_coop$n_to_remove, "species\n")
} # should be 24
if (!is.null(downsampling_results$downsampling$global_coop)) {
  cat("- Global cooperative:", downsampling_results$downsampling$global_coop$n_to_remove, "species\n")
} # should be 15
if (!is.null(downsampling_results$downsampling$territoriality_bias)) {
  cat("- Territoriality bias:", downsampling_results$downsampling$territoriality_bias$calculation, "\n")
} # should be 266
if (!is.null(downsampling_results$downsampling$territory_12vs3_bias)) {
  cat("- Territory_12vs3 bias:", downsampling_results$downsampling$territory_12vs3_bias$calculation, "\n")
} # should be 155


#### Run phylopath analyses with downsampling - binary traits ----

# Must first run calculate_stratified_downsampling() in above section

source("run_phylopath_fxns.R")
#source("calculate_downsampling_function.R") # calculate_downsampling_function.R probably obsolete because of calculate_stratified_downsampling(), but maybe it would be better to keep/use this one instead? 

dfIn_phylo = read.csv("Data_R_2025-06-09.csv")
tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")

# Set number of iterations (use 500 for publication, 50 for testing)
n_iterations <- 500  # Change to 50 for testing

female_song_var = "FemaleSong_Agg01"
coop_breeding_var = "HighConfidence_Coop"
territoriality_var = "Territory_12vs3" # "TerritorialityWeakVsStrong"
mass_var = "logMass_AVONET"

all_traits_phylopath_label = paste(female_song_var, coop_breeding_var, territoriality_var, mass_var)

trait_set_output_dir = file.path("Outputs", "PhylopathDownsampled", paste0(all_traits_phylopath_label, " models"))

if (!dir.exists(trait_set_output_dir)) {
  dir.create(trait_set_output_dir, recursive = T)
}

# Add geographic regions if needed
if (!"GeographicRegion_Jetz" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$GeographicRegion_Jetz <- NA
  dfIn_phylo$GeographicRegion_Jetz[which(dfIn_phylo$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
  dfIn_phylo$GeographicRegion_Jetz[which(dfIn_phylo$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
}

# Add data availability indicators
if (!"HaveFSData" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$HaveFSData <- !is.na(dfIn_phylo$FemaleSong_Agg01)
}
if (!"HaveFSCBData" %in% colnames(dfIn_phylo)) {
  dfIn_phylo$HaveFSCBdata <- !is.na(dfIn_phylo$HighConfidence_Coop) & !is.na(dfIn_phylo$FemaleSong_Agg01)
}

################################# 1. Geographic bias correction - HOLARCTIC NONCOOPERATIVE ---
cat("\n\nRunning geographic bias correction - Holarctic non-cooperative...\n")
cat("Removing", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "species\n")
result_geo_holarctic <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
  downsample_values = c("Holarctic", 0, TRUE),
  numToRemove = downsampling_results$downsampling$holarctic_noncoop$n_to_remove,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "HolarcticNoncoop")
)

# Create plots with proper file naming
geo_holarctic_info <- list(
  proportions = data.frame(
    territoriality = c("Holarctic Non-coop", "Tropical Non-coop"),
    species_with_data = c(209, 428),
    total_species = c(369, 977),
    proportion = c(0.566, 0.438)
  ),
  target_proportion = 0.438,
  n_to_remove = downsampling_results$downsampling$holarctic_noncoop$n_to_remove, #83,
  downsample_info = list(
    group_to_downsample = "Holarctic non-cooperative",
    territory_value = "0"
  )
)

prefix_geo <- paste(all_traits_phylopath_label, "Remove83HolarcticNoncoop")
saveRDS(result_geo_holarctic,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_geo, "_n", n_iterations, "_", Sys.Date(), ".rds")))

plots_geo <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_geo_holarctic,
  downsampling_info = geo_holarctic_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_prefix = prefix_geo,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 2. Geographic bias correction - TROPICAL COOPERATIVE ---
cat("\n\nRunning geographic bias correction - Tropical cooperative...\n")
result_geo_tropical <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("GeographicRegion_Jetz", "HighConfidence_Coop", "HaveFSData"),
  downsample_values = c("Tropical", 1, TRUE),
  numToRemove = downsampling_results$downsampling$tropical_coop$n_to_remove, #24,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove24TropicalCoop")
)

# Create plots
trop_coop_info <- list(
  proportions = data.frame(
    group = c("Tropical Coop", "Other"),
    species_with_data = c(68, 807),
    total_species = c(159, 2005),
    proportion = c(0.428, 0.402)
  ),
  target_proportion = 0.402,
  n_to_remove = downsampling_results$downsampling$tropical_coop$n_to_remove,
  downsample_info = list(
    group_to_downsample = "Tropical cooperative",
    territory_value = "1"
  )
)

prefix_trop <- paste(all_traits_phylopath_label, "Remove24TropicalCoop") 
saveRDS(result_geo_tropical,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_trop, "_n", n_iterations, "_", Sys.Date(), ".rds")))

plots_trop <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_geo_tropical,
  downsampling_info = trop_coop_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_prefix = prefix_trop,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 3. GLOBAL COOPERATIVE bias correction ---
cat("\n\nRunning global cooperative bias correction...\n")
result_global_coop <- run_multiple_phylopath(
  dfIn = dfIn_phylo,
  tree = tree_phylo,
  downsample_columns = c("HighConfidence_Coop", "HaveFSData"),
  downsample_values = c(1, TRUE),
  numToRemove = downsampling_results$downsampling$global_coop$n_to_remove, #15,
  n_iterations = n_iterations,
  female_song_var = female_song_var,
  coop_breeding_var = coop_breeding_var,
  territoriality_var = territoriality_var,
  mass_var = mass_var,
  save_conditional_plots = TRUE,
  plotlabel = paste0("Remove15GlobalCoop")
)

# Create plots - info
global_coop_info <- list(
  proportions = data.frame(
    group = c("Cooperative", "Non-cooperative"),
    species_with_data = c(80, 795),
    total_species = c(226, 1938),
    proportion = c(0.354, 0.410)
  ),
  target_proportion = 0.069,
  n_to_remove = downsampling_results$downsampling$global_coop$n_to_remove, #15,
  downsample_info = list(
    group_to_downsample = "Cooperative",
    territory_value = "1"
  )
)

# actual plot creation
prefix_global <- paste(all_traits_phylopath_label, "Remove15GlobalCoop")
saveRDS(result_global_coop,
        file.path(trait_set_output_dir,
                  paste0("result_", prefix_global, "_n", n_iterations, "_", Sys.Date(), ".rds")))
plots_global <- create_all_phylopath_plots(
  analysis_type = "downsampled",
  phylopath_output = result_global_coop,
  downsampling_info = global_coop_info,
  full_dataset = "Data_R_2025-06-09.csv",
  tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
  output_prefix = prefix_global,
  output_dir = trait_set_output_dir,
  save_png = TRUE,
  save_pdf = TRUE
)

################################# 4. TERRITORIALITY WEAK/STRONG bias correction ---
cat("\n\nRunning territoriality bias correction...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
dfIn_phylo_char$TerritorialityWeakVsStrong <- as.character(dfIn_phylo_char$TerritorialityWeakVsStrong)
dfIn_phylo_char$Territory_12vs3 <- as.character(dfIn_phylo_char$Territory_12vs3)

if (downsampling_results$downsampling$territoriality_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territoriality_bias$n_to_remove, "species from StrongTerr\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("TerritorialityWeakVsStrong", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territoriality_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territoriality_bias$n_to_remove, "StrongTerr")
  )
  
  # Create plots
  prefix_terr <- paste0(all_traits_phylopath_label, "Remove",
 downsampling_results$downsampling$territoriality_bias$n_to_remove, "StrongTerr")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = downsampling_results$downsampling$territoriality_bias$n_to_remove,
    full_dataset = "Data_R_2025-06-09.csv",
    tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

################################# 5. TERRITORIALITY 12 VS 3 bias correction using Territory_12vs3 as TERR variable ---
cat("\n\nRunning Territory_12vs3 bias correction...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
dfIn_phylo_char$Territory_12vs3 <- as.character(dfIn_phylo_char$Territory_12vs3)

if (downsampling_results$downsampling$territory_12vs3_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "species from", "Terr3\n")
     # terr_downsample_info$report$downsample_info$group_to_downsample, "\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("Territory_12vs3", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territory_12vs3_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = "Territory_12vs3",
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  )
  
  # Create plots
  prefix_terr <- paste0("FemaleSong_Agg01 HighConfidence_Coop Territory_12vs3 logMass_AVONET Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = NULL, #terr_downsample_info$report,
    full_dataset = "Data_R_2025-06-09.csv",
    tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}


################################# 6. TERRITORIALITY 12 VS 3 bias correction using TerritorialityWeakVsStrong as TERR variable ---
cat("\n\nRunning Territory_12vs3 bias correction - uses TerritorialityWeakVsStrong as TERR var...\n")

# Convert territoriality to character for calculation
dfIn_phylo_char <- dfIn_phylo
#dfIn_phylo_char$Terr <- as.character(dfIn_phylo_char$Territory_12vs3)
dfIn_phylo_char$TerritorialityWeakVsStrong <- as.character(dfIn_phylo_char$TerritorialityWeakVsStrong)

if (downsampling_results$downsampling$territory_12vs3_bias$n_to_remove > 0) {
  cat("Need to remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "species from Terr3", "\n")
  
  result_terr <- run_multiple_phylopath(
    dfIn = dfIn_phylo,
    tree = tree_phylo,
    downsample_columns = c("Territory_12vs3", "HaveFSCBdata"),
    downsample_values = c(1, TRUE),
    numToRemove = downsampling_results$downsampling$territory_12vs3_bias$n_to_remove,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = "TerritorialityWeakVsStrong",
    mass_var = mass_var,
    save_conditional_plots = TRUE,
    plotlabel = paste0("Remove", downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3-uses-TerritorialityWeakVsStrong-as-TERR")
  )
  
  # Create plots
  prefix_terr <- paste0("FemaleSong_Agg01 HighConfidence_Coop TerritorialityWeakVsStrong logMass_AVONET Remove",
                        downsampling_results$downsampling$territory_12vs3_bias$n_to_remove, "Terr3")
  plots_terr <- create_all_phylopath_plots(
    analysis_type = "downsampled",
    phylopath_output = result_terr,
    downsampling_info = NULL, #terr_downsample_info$report,
    full_dataset = "Data_R_2025-06-09.csv",
    tree = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex",
    output_prefix = prefix_terr,
    output_dir = trait_set_output_dir,
    save_png = TRUE,
    save_pdf = TRUE
  )
}

# 5. Process all detailed model files with enhanced plots
cat("\n\nProcessing detailed model files for enhanced plots...\n")
detailed_files <- list.files(path = trait_set_output_dir, pattern = paste0("detailed_models_.*", ".*\\.csv$"), 
                             full.names = TRUE)

## Obsolete due to no create_enhanced_phylopath_plots() anymore, but perhaps worth using this framework for plotting with create_all_phylopath_plots()
# for (file in detailed_files) {
#   cat("Processing:", basename(file), "\n")
#   
#   scenario <- gsub("detailed_models_", "", basename(file))
#   scenario <- gsub(paste0("", ".*\\.csv$"), "", scenario)
#   
#   enhanced_plots <- create_enhanced_phylopath_plots( #OBSOLETE
#     csv_file = file,
#     output_prefix = file.path(trait_set_output_dir, paste0("enhanced_", scenario)),
#     save_png = TRUE,
#     save_pdf = FALSE
#   )
# }

cat("\n\nAll phylopath analyses complete!\n")
cat("Results saved to:", trait_set_output_dir, "\n")


#### Run phylopath with dimorphism bias correction/downsampling ----
## Run early AM 6/10/2025 - seems to be working? Except phylopath plots 
## ALMOST READY TO DELETE - USE JUST "ALTERNATIVE" SECTION BELOW
# source("run_phylopath_fxns.R") # this section uses downsample_dimorphism_bias() found in the main phylopath fxns script
# 
# dfIn_phylo = read.csv("Data_R_2025-06-09.csv")
# tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
# phylopath_output_dir = "Outputs/PhylopathDownsampled"
# n_iterations = 5
# 
# # Define dimorphism variables to test
# dimorphism_vars <- list(
#   plumage = list(
#     col = "logMaleFemalePlumageDiffAbs",
#     label = "PlumageDimorphism",
#     description = "log Plumage Dimorphism (Absolute Value)"
#   ),
#   wing = list(
#     col = "PercentAbsLogWingDimorphism",
#     label = "WingDimorphism", 
#     description = "Percent Absolute-Value Log Wing Dimorphism"
#   )
# )
# 
# # Loop through each dimorphism variable
# for (dim_type in names(dimorphism_vars)) {
#   dim_info <- dimorphism_vars[[dim_type]]
#   
#   # Create output directory specific to this dimorphism type
#   dim_output_dir <- file.path(phylopath_output_dir, dim_info$label)
#   if (!dir.exists(dim_output_dir)) {
#     dir.create(dim_output_dir, recursive = TRUE)
#   }
#   
#   cat("\n\n========================================\n")
#   cat("Running", dim_info$description, "bias correction...\n")
#   cat("========================================\n")
#   
#   # Check if the column exists
#   if (!dim_info$col %in% colnames(dfIn_phylo)) {
#     cat("Warning: Column", dim_info$col, "not found. Skipping...\n")
#     next
#   }
#   
#   # Run the dimorphism downsampling analysis
#   dimorphism_downsample <- downsample_dimorphism_bias( 
#     df = dfIn_phylo,
#     dimorphism_col = dim_info$col,
#     data_col = "FemaleSong_Agg01",
#     n_iterations = n_iterations
#   )
#   
#   dist_plot_result <- create_downsampled_dimorphism_distribution_plot( # Working 6/10/2025, but not useful?
#     dfIn_phylo = dfIn_phylo,
#     dimorphism_downsample = dimorphism_downsample,
#     dim_info = dim_info,
#     output_dir = dim_output_dir,
#     prefix = prefix_dimorphism,
#     save_plot = TRUE
#   )
#   
#   # Create the violin plot (multiple iterations view)
#   violin_plot_result <- create_dimorphism_downsamples_violin_plot(
#     dfIn_phylo = dfIn_phylo,
#     dimorphism_downsample = dimorphism_downsample,
#     dim_info = dim_info,
#     output_dir = dim_output_dir,
#     prefix = prefix_dimorphism,
#     save_plot = TRUE,
#     show_iterations = 20
#   )
#   
#   cat("Need to remove", dimorphism_downsample$n_to_remove, "species to correct", dim_info$description, "bias\n")
#   cat("Current mean:", round(dimorphism_downsample$current_stats$current_mean, 3), "\n")
#   cat("Target mean:", round(dimorphism_downsample$target_stats$target_mean, 3), "\n")
#   
#   # Run phylopath on each iteration
#   all_results <- list()
#   all_model_summaries <- list()
#   detailed_models_list <- list()
#   
#   for (i in 1:n_iterations) {
#     # Get the downsampled dataset for this iteration
#     df_downsampled <- dimorphism_downsample$iterations[[i]]$remaining_data
#     
#     # Filter to species in tree
#     df_downsampled <- df_downsampled[df_downsampled$species %in% tree_phylo$tip.label, ]
#     rownames(df_downsampled) <- df_downsampled$species
#     
#     # Check if we have enough species
#     if (nrow(df_downsampled) < 50) {
#       warning(paste("Iteration", i, "has only", nrow(df_downsampled), "species. Skipping..."))
#       next
#     }
#     
#     # Run phylopath on this iteration
#     result_i <- run_CB_FS_Terr_phylopath(
#       dfIn = df_downsampled,
#       tree = tree_phylo,
#       female_song_var = "FemaleSong_Agg01",
#       coop_breeding_var = "HighConfidence_Coop",
#       territoriality_var = "TerritorialityWeakVsStrong",
#       mass_var = "logMass_AVONET",
#       plots2pdf = FALSE
#     )
#     
#     all_results[[i]] <- result_i
#     
#     # Extract model summaries from the phylopath result
#     if (!is.null(result_i$result)) {
#       # Get the summary of the phylopath result
#       summary_i <- summary(result_i$result)
#       
#       # Convert to data frame and add iteration number
#       if (!is.null(summary_i)) {
#         model_summary_i <- as.data.frame(summary_i)
#         model_summary_i$iteration <- i
#         model_summary_i$model_name <- rownames(model_summary_i)
#         all_model_summaries[[i]] <- model_summary_i
#       }
#     } # end if !is.null(result_i$result)
#   } # end for i in 1:n_iterations
#   
#   # Aggregate results across iterations
#   aggregated_summaries <- do.call(rbind, all_model_summaries)
#   
#   # Calculate average model performance
#   avg_model_performance <- aggregated_summaries %>%
#     group_by(model_name) %>%
#     summarise(
#       mean_CICc = mean(CICc, na.rm = TRUE),
#       sd_CICc = sd(CICc, na.rm = TRUE),
#       times_best = sum(delta_CICc == 0, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     arrange(mean_CICc)
#   
#   # Create a summary phylopath output structure
#   result_dimorphism <- list(
#     model_summaries = aggregated_summaries,
#     avg_performance = avg_model_performance,
#     n_iterations = n_iterations,
#     downsampling_info = dimorphism_downsample,
#     all_results = all_results,
#     dimorphism_type = dim_info$label
#   )
#   
#   # Create plots
#   cat("\nCreating", dim_info$description, "bias correction plots...\n")
#   
#   dimorphism_plot_info <- list(
#     proportions = data.frame(
#       group = c("No FS Data", "Has FS Data"),
#       mean_dimorphism = c(dimorphism_downsample$target_stats$target_mean,
#                           dimorphism_downsample$current_stats$current_mean),
#       total_species = c(sum(is.na(dfIn_phylo$FemaleSong_Agg01)),
#                         sum(!is.na(dfIn_phylo$FemaleSong_Agg01))),
#       species_with_data = c(sum(is.na(dfIn_phylo$FemaleSong_Agg01)),
#                             sum(!is.na(dfIn_phylo$FemaleSong_Agg01)))  # Same as total for this use case
#     ),
#     target_mean = dimorphism_downsample$target_stats$target_mean,
#     n_to_remove = dimorphism_downsample$n_to_remove,
#     downsample_info = list(
#       group_to_downsample = paste("High", dim_info$description, "species"),
#       method = "propensity-weighted"
#     )
#   )
#   
#   # Create output prefix with dimorphism type
#   prefix_dimorphism <- paste0("Remove", 
#                               dimorphism_downsample$n_to_remove, 
#                               "High", 
#                               dim_info$label)
#   
#   plots_dimorphism <- create_all_phylopath_plots( # not working 6/10/2025
#     analysis_type = "downsampled",
#     phylopath_output = result_dimorphism,
#     downsampling_info = NULL,
#     output_prefix = prefix_dimorphism,
#     output_dir = dim_output_dir,
#     save_png = TRUE
#   )
#   
#   # Save detailed results with dimorphism type in filename
#   write.csv(
#     aggregated_summaries,
#     file = file.path(dim_output_dir, 
#                      paste0("detailed_models_", dim_info$label, "", ".csv")),
#     row.names = FALSE
#   )
#   
#   # Save downsampling summary
#   downsample_summary <- data.frame(
#     dimorphism_type = dim_info$label,
#     dimorphism_variable = dim_info$col,
#     n_removed = dimorphism_downsample$n_to_remove,
#     target_mean = dimorphism_downsample$target_stats$target_mean,
#     original_mean = dimorphism_downsample$current_stats$current_mean,
#     corrected_mean = dimorphism_downsample$summary$mean_dimorphism_after,
#     convergence = dimorphism_downsample$summary$convergence
#   )
#   
#   write.csv(
#     downsample_summary,
#     file = file.path(dim_output_dir, 
#                      paste0("downsampling_summary_", dim_info$label, ".csv")),
#     row.names = FALSE
#   )
#   
#   cat("\n", dim_info$description, "bias correction complete!\n")
#   cat("Mean after correction:", 
#       round(dimorphism_downsample$summary$mean_dimorphism_after, 3), "\n")
#   cat("Results saved to:", dim_output_dir, "\n")
# }
# 
# cat("\n\nAll dimorphism bias corrections complete!\n")


#### ALTERNATIVE RUN DIMORPH DOWNSAMPLING ----
# Source the improvements file
source("run_phylopath_fxns.R") # this section uses downsample_dimorphism_bias() found in the main phylopath fxns script

dfIn_phylo = read.csv("Data_R_2025-06-09.csv")
tree <- tree_phylo <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
phylopath_output_dir = file.path("Outputs","PhylopathDownsampled", paste(all_traits_phylopath_label, "models"))
n_iterations = 500

# Define dimorphism variables to test
dimorphism_vars <- list(
  plumage = list(
    col = "logMaleFemalePlumageDiffAbs",
    label = "PlumageDimorphism",
    description = "log Plumage Dimorphism (Absolute Value)"
  ),
  wing = list(
    col = "PercentAbsLogWingDimorphism",
    label = "WingDimorphism", 
    description = "Percent Absolute-Value Log Wing Dimorphism"
  )
)

# Loop through each dimorphism variable
for (dim_type in names(dimorphism_vars)) {
  dim_info <- dimorphism_vars[[dim_type]]
  
  # Run the improved phylopath dimorphism correction
  result_dimorphism <- run_phylopath_dimorphism_correction(
    dfIn_phylo = dfIn_phylo,
    tree = tree_phylo,
    dim_info = dim_info,
    n_iterations = n_iterations,
    female_song_var = female_song_var,
    coop_breeding_var = coop_breeding_var,
    territoriality_var = territoriality_var,
    mass_var = mass_var,
    phylopath_output_dir = phylopath_output_dir,
    save_outputs = TRUE
  )
  
  # If results were obtained, create plots
  if (!is.null(result_dimorphism)) {
    # Create output prefix with dimorphism type
    prefix_dimorphism <- paste0("Remove", 
                                result_dimorphism$downsampling_info$n_to_remove, 
                                "High", 
                                dim_info$label,"ByPropensity")
    saveRDS(result_dimorphism,
            file.path(trait_set_output_dir,
                      paste0("result_", prefix_dimorphism, "_n", n_iterations, "_", Sys.Date(), ".rds")))
    
    # Use the flexible plotting function
    plots_dimorphism <- create_downsampled_plots(
      downsampling_results = result_dimorphism,
      detailed_models_input = result_dimorphism$detailed_models,  # Pass dataframe directly
      downsampling_info = NULL,
      full_dataset = dfIn_phylo,
      tree = tree_phylo,
      full_data_phylopath_input = NULL,
      output_prefix = prefix_dimorphism,
      output_dir = file.path(phylopath_output_dir, dim_info$label),
      save_png = TRUE,
      save_pdf = TRUE
    )
  }
}

cat("\n\nAll dimorphism bias corrections complete!\n")


#### PhyloGLM ----
