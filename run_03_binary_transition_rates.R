# Test models of transition rates for binary traits ----

# Supplemental Table 19 - ARD vs ER rates for binary traits ----
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
write.csv(allQout, file.path("Outputs","binary trait Qrates_rounded.csv"))
