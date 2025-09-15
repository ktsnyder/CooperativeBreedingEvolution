# Run Brownie analyses


source(file.path("Brownie_functions","browniefunction.R"))
source(file.path("Brownie_functions","plotbrownie.R"))
source(file.path("Brownie_functions","jackknifingbrownie.R"))
source(file.path("Brownie_functions","brownie relative rates.R"))
source(file.path("Brownie_functions","BrownieMultistate.R"))
source(file.path("Brownie_functions","brownie relative rates.R"))
source("findQrates.R")
source("getLabels.R")

require(dplyr)

newdata = "Data_R.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
currentlabel = ""

# Figure 1, Extended Data Figure 1, Supplemental Table 4 - Brownie ----
nsim = 20 # for demonstration; increase to 500 for full analyses
songtraits = c("Song.rep.final","Syllable.rep.final") # "Syll.song.final", "Duration.final", "Interval.final")
socialtraits = c("HighConfidence_Coop") #, "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Territory_12vs3", "TerritorialityWeakVsStrong")
currentlabel = ""

dir.create(file.path("Outputs", "Brownie_outputs"), recursive = T)

for (l in 1:length(socialtraits)) {
  for (k in 1:length(songtraits)) {
    print(Sys.time())
    feature <- songtraits[k]
    CBcolumn <- socialtraits[l]
    discreteCatLabels = getLabels(CBcolumn)
    print(feature)
    browniefunction(columns = c(CBcolumn, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = FALSE, otherlabel = currentlabel)
    
    if (file.exists(file.path("Outputs", "Brownie_outputs",paste0(Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv")))) {
      print("file exists")
      plotbrownie(data = paste0(Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
    } else if (file.exists(file.path("Outputs", "Brownie_outputs", paste0(Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv")))) {
      plotbrownie(data = paste0(Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
      print("yesterday's file exists")
    } else {
      print("file does not exist")
      print(file.path("Outputs", "Brownie_outputs",paste0(Sys.Date()-1,CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv")))
    }
  } # end for k
} # end for l
# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = file.path("Outputs", "Brownie_outputs"), otherlabel = currentlabel)


# Supplemental Table 5 - Brownie Jackknife ----
subset = subsettreedata(columns = c("HighConfidence_Coop", "Song.rep.final"), newtree = treefile, newdata = newdata)
familycounts = subset$subsetdf %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familysubset = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 20)]
currentlabel = "jackknife"

jackbrowniefunction(columns = c("HighConfidence_Coop", "Song.rep.final"), islog = T, matemodel = "ARD", matensim = nsim, allcsvs = T, plotsimmaps = F, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2_AVONET", otherlabel = currentlabel, cladeJackvalues = familysubset)

# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = file.path("Outputs","BrownieJackknifeOutputs"), otherlabel = currentlabel)


# Extended Data Figure 2; Supplemental Table 6 - Brownie with multi-state categorical traits ----
source(file.path("Brownie_functions","BrownieMultistate.R"))
BrownieMultistate(DiscreteTrait = "grouping_Griesser2023", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = nsim, plotsimmaps = F, plotResults = T, otherlabel = "test")
BrownieMultistate(DiscreteTrait = "social_system_incl_nk_coop_Griesser2017", ContinuousTrait = "Song.rep.final", newdata = newdata, treefile = treefile, nsim = nsim, plotsimmaps = F, plotResults = T, otherlabel = "test")

df_brown = read.csv(newdata)
df_brown$TerritorialityWeakVsStrong[which(df_brown$TerritorialityWeakVsStrong == 0)] <- "Weak"
df_brown$TerritorialityWeakVsStrong[which(df_brown$TerritorialityWeakVsStrong == 1)] <- "Strong"
df_brown$HighConfidence_Coop[which(df_brown$HighConfidence_Coop == 0)] <- "Noncooperative"
df_brown$HighConfidence_Coop[which(df_brown$HighConfidence_Coop == 1)] <- "Cooperative"
df_brown$TerrWeakStrongXHighConfCoop <- paste(df_brown$HighConfidence_Coop, df_brown$TerritorialityWeakVsStrong, sep = "_")
require(stringr)
df_brown$TerrWeakStrongXHighConfCoop[which(str_detect(df_brown$TerrWeakStrongXHighConfCoop, "NA"))] <- NA

BrownieMultistate(DiscreteTrait = "TerrWeakStrongXHighConfCoop", ContinuousTrait = "Song.rep.final", newdata = df_brown, treefile = treefile, nsim = nsim, plotsimmaps = F, plotResults = T, otherlabel = "FourStateSimmaps")
