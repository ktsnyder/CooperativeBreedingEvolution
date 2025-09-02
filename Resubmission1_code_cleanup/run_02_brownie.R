# Run Brownie analyses

source("browniefunction.R")
source("plotbrownie.R")
source("jackknifingbrownie.R")
source("brownie relative rates.R")
source("BrownieMultistate.R")
source("brownie relative rates.R")
source("findQrates.R")

require(dplyr)

newdata = "Resubmission1_code_cleanup/Data_R.csv"
treefile = "Resubmission1_code_cleanup/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

# Figure 1, Supplemental Figure 1, Supplemental Table 4 - Brownie ----
nsim = 50 # for demonstration; increase to 500 for full analyses
songtraits = c("Song.rep.final","Syllable.rep.final") # "Syll.song.final", "Duration.final", "Interval.final")
socialityMetrics = c("HighConfidence_Coop") #, "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Territory_12vs3", "TerritorialityWeakVsStrong")

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
subset = subsettreedata(columns = c("HighConfidence_Coop", "Song.rep.final"), newtree = treefile, newdata = newdata)
familycounts = subset$subsetdf %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familysubset = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 20)]

jackbrowniefunction(columns = c("HighConfidence_Coop", "Song.rep.final"), islog = T, matemodel = "ARD", matensim = 10, allcsvs = T, plotsimmaps = F, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2_AVONET", otherlabel = currentlabel, cladeJackvalues = familysubset)

# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = "BrownieJackknifeOutputs", otherlabel = currentlabel)


# Supplemental Figures 2 & 3; Supplemental Tables 6 & 7 - Brownie with multi-state categorical traits ----
source("Resubmission1_code_cleanup/Brownie_functions/BrownieMultistate.R")
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

nsims = 100

BrownieMultistate(DiscreteTrait = "TerrWeakStrongXHighConfCoop", ContinuousTrait = "Song.rep.final", newdata = df_brown, treefile = treefile, nsim = nsims, plotsimmaps = T, plotResults = T, otherlabel = "FourStateSimmaps")
