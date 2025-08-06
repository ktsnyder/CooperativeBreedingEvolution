realdf = read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap_Overlap_Outputs_FemaleSong_Agg01_vs_TerritorialityWeakVsStrong_20250627_101945/FemaleSong_Agg01 TerritorialityWeakVsStrong REAL simmap overlap_counts output nsim 150 Re-refactored .csv')

dummydf = read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap_Overlap_Outputs_FemaleSong_Agg01_vs_TerritorialityWeakVsStrong_20250627_101945/FemaleSong_Agg01 TerritorialityWeakVsStrong DUMMYResampledMkSimmap simmap overlap_counts output nsim 150 Re-refactored .csv')
realdf$column2


dfout = read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/TerritorialityWeakVsStrong FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 Re-refactored .csv')
dfDummy = read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap_Overlap_Outputs_TerritorialityWeakVsStrong_vs_FemaleSong_Agg01_500_20250627_183157/TerritorialityWeakVsStrong FemaleSong_Agg01 DUMMYResampledMkSimmap simmap overlap_counts output nsim 500 Re-refactored .csv')

plot(density(realdf$trait2_prop_state1))
plot(density(dummydf$trait2_prop_state1))


Overlap_0_0_FractionDummyLessThanMedianReal = sum(dfDummy$Overlap_0_0 <= median(dfout$Overlap_0_0))/length(dfDummy$treenum)
Overlap_0_1_FractionDummyLessThanMedianReal = sum(dfDummy$Overlap_0_1 <= median(dfout$Overlap_0_1))/length(dfDummy$treenum)
Overlap_1_0_FractionDummyLessThanMedianReal = sum(dfDummy$Overlap_1_0 <= median(dfout$Overlap_1_0))/length(dfDummy$treenum)
Overlap_1_1_FractionDummyLessThanMedianReal = sum(dfDummy$Overlap_1_1 <= median(dfout$Overlap_1_1))/length(dfDummy$treenum)


ObsProp0Absent_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp0Absent <= median(dfout$ObsProp0Absent))/length(dfDummy$treenum)
ObsProp0Present_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp0Present <= median(dfout$ObsProp0Present))/length(dfDummy$treenum)
ObsProp1Absent_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp1Absent <= median(dfout$ObsProp1Absent))/length(dfDummy$treenum)
ObsProp1Present_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp1Present <= median(dfout$ObsProp1Present))/length(dfDummy$treenum)

Nspecies = 919
pval = 0
# Create outputs for overall multi-comparison table - from test_trait_overlap_simmaps.R
# Create outputs for overall multi-comparison table
dfout <- dfout %>%
  mutate(across(coopQ01:ObsProp1Present, as.numeric))
medians = dfout %>% select(coopQ01:ObsProp1Present) %>% summarise_all(median, na.rm = TRUE)
names(medians)[6:13] <- paste0(names(medians)[6:13],"_MedianReal")

# quartiles <- dfout %>%
#   select(ObsProp0Absent:ObsProp1Present) %>%
#   summarise_all(function(x) list(quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)))
# quartiles <- unnest(quartiles, cols = everything())


mediansReal <- cbind(trait1, trait2, Nspecies, nsims_real, nsims_dummy, pval, medians)
dfDummy <- dfDummy %>%
  mutate(across(coopQ01:ObsProp1Present, as.numeric))
mediansDummy = dfDummy %>% select(propFSabsent:ObsProp1Present) %>% summarise_all(median, na.rm = TRUE) 
names(mediansDummy) <- paste0(names(mediansDummy), "_MedianDummy")
mediansRow = cbind(mediansReal, mediansDummy)

require(emmeans)
#require(ggpubr)
lmStates = lm(ObservedState.prop ~ Which*ObservedState, data = dfCombined)
emm <- emmeans(lmStates, pairwise ~ Which | ObservedState)
contrast <- emm$contrasts
PairwisePostHoc = summary(contrast, adjust = "tukey")
PairwisePvals = PairwisePostHoc$p.value
PairwisePvals = as.data.frame(t(as.data.frame(PairwisePvals)))
colnames(PairwisePvals) <- paste0(PairwisePostHoc$ObservedState, "_DummyVsRealPval")
mediansRow = cbind(mediansRow, PairwisePvals)
mediansRow = cbind(mediansRow, ObsProp0Absent_FractionDummyLessThanMedianReal, ObsProp0Present_FractionDummyLessThanMedianReal, ObsProp1Absent_FractionDummyLessThanMedianReal, ObsProp1Present_FractionDummyLessThanMedianReal)